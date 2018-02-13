from multiplierz.mzAPI import mzFile
import matplotlib.pyplot as pyt
from collections import defaultdict, deque
import time
import multiplierz.mzReport as mzReport
import os
import cPickle as pickle
import re
from multiplierz.internalAlgorithms import ProximityIndexedSequence, inPPM, average, pts_to_bins
from multiplierz import vprint, verbose_mode

from multiplierz.mzTools.featureUtilities import save_feature_database, FeatureInterface
from multiplierz.internalAlgorithms import peak_pick_PPM
import multiprocessing

from multiplierz.mgf import standard_title_parse

try:
    from scipy.stats import kurtosis, skew
except ImportError:
    kurtosis = lambda x: 'NA'
    skew = lambda x: 'NA'

__all__ = ['feature_analysis', 'detect_features']

# These (and the spectrumDescriptionTo... functions) will be modified when invoked by the GUI.
signalToNoiseThreshold = 15
#peakFindTolerance = 0.02

featureMatchupTolerance = 0.05

whitelist_tol = 0.1
def spectrumDescriptionToMZ(description):
    if 'MultiplierzMGF' in description:
        return float(standard_title_parse(description)['mz'])
    try:
        return float(description.split('-|')[0].split('-')[-1])
    except ValueError:
        return float(description.split("|")[1])
    
def mzFromPSM(psm):
    if 'Spectrum Description' in psm:
        return spectrumDescriptionToMZ(psm['Spectrum Description'])
    else:
        return float(psm['m/z [Da]']) # Proteome Discoverer output.

def spectrumDescriptionToScanNumber(description):
    if 'MultiplierzMGF' in description:
        return int(standard_title_parse(description)['scan'])
    else:
        return int(description.split('.')[1])

def scanFromPSM(psm):
    if 'Spectrum Description' in psm:
        return spectrumDescriptionToScanNumber(psm['Spectrum Description'])
    else:
        return int(psm['First Scan']) # Proteome Discoverer output.


#allowedC12RatioDifference = 0.3
allowedC12RatioDifference = 2.0
#c12RelativeIntensity = 0.6
#c12RelativeLimit = 2
isotopicRatioTolerance = 2
noisePeaksAllowed = 1000000

#splitDueToIntensity = False
featureAbsenceTolerance = 10
dropoutTimeTolerance = 0.5

def unzip(thing): return [list(x) for x in zip(*thing)]

#curveargs = [-1.02857097, 0.000113693166, 8.53554707]

#def getC12Ratio(mz, charge):
    #mass = mz * charge # Very approximately.
    #return np.log(mass) * curveargs[0] + mass * curveargs[1] + curveargs[2]


class Feature():
    def __init__(self):
        self.regions = []
        self.spectrum = None
        self.scanrange = None
        self.allmzs = None
        self.wasSplit = False
    def add(self, datapoints, index, charge):
        self.regions.append((index, datapoints))
        self.charge = charge
    def allIndexedPoints(self):
        try:
            return sum([[(i, x, y, c, n) for (x, y, c, n) in region] for i, region in self.regions], [])
        except ValueError: # Not lscan-derived points.
            return sum([[(i, x, y) for (x, y) in region] for i, region in self.regions], [])
    def strengthAt(self, index):
        try:
            return sum([x[1] for x in self.regions[index][1]])
        except IndexError:
            return 0
    def topSignal(self):
        return max([self.strengthAt(x) for x in range(0, len(self.regions))])
    def totalIntensity(self):
        return sum([self.strengthAt(x) for x in range(0, len(self.regions))])
    def c12Intensity(self):
        return sum([min(x, key = lambda x: x[0])[1] for _, x in self.regions])
    def segment(self, start, end):
        subFeature = Feature()
        subFeature.regions = self.regions[start:end]
        return subFeature
    def scanxic(self):
        return [(r[0], r[1][0][1]) for r in sorted(self.regions)]
    # Previous version didn't have absolute-scan-index conversion on regions.
    def calculate_bounds(self, absolute_scan_lookup):
        regionIndices = zip(*self.regions)[0]
        self.scans = [absolute_scan_lookup[x] for x in regionIndices]
        self.scanrange = min(self.scans), max(self.scans)
        
        self.regions = [(absolute_scan_lookup[x], y) for x, y in self.regions]
        
        minMZs = [min(zip(*peaks)[0]) for y, peaks in self.regions]
        avgC12 = average(minMZs)
        while not all([abs(x - avgC12) < 0.05 for x in minMZs]):
            oddOneOut = max(minMZs, key = lambda x: abs(x - avgC12))
            minMZs.remove(oddOneOut)
            avgC12 = average(minMZs)
        self.mz = avgC12
        
        xic = pts_to_bins([(x[0], min(x[1], key = lambda pt: (pt[0] - avgC12))[1]) for x in self.regions],
                          100)
        ints = [0.0] + list(zip(*xic)[1]) + [0.0]
        self.skewness = skew(ints)
        self.kurtosis = kurtosis(ints)
        
        
  
    def containsPoint(self, mz, scan, charge):
        return (charge == self.charge and abs(mz - self.mz) < featureMatchupTolerance
                and self.scanrange[0] < scan < self.scanrange[1])
        
    def tidyData(self):
        self.allmzs = sum([[x[0] for x in points] for (index, points) in self.regions], [])
        
    def bordersPoint(self, mz, scan, charge):
        #if not self.scanrange:
            #self.scanrange = min(self.scans) - 1, max(self.scans) + 1
        if not self.allmzs:
            self.allmzs = sum([[x[0] for x in points] for (index, points) in self.regions], [])
        
        #assert not self.containsPoint(mz, scan, charge)
        if self.containsPoint(mz, scan, charge):
            return "Contained."
        
        if charge != self.charge:
            return ""
        
        edge = []
        
        if self.scanrange[0] == scan and abs(mz - self.mz) < 0.05:
            edge.append("Scan before feature")
        if self.scanrange[1] == scan and abs(mz - self.mz) < 0.05:
            edge.append("Scan after feature")
        if any([abs(x - mz) < 0.05 for x in self.allmzs]) and scan in self.scans:
            edge.append("Non-C12 peak")
        
        if edge:
            return '; '.join(edge)
        else:
            return ""
      
      
      
def setGlobals(constants):
    if 'mzRegex' in constants:
        global spectrumDescriptionToMZ
        
        mzRegCompiled = re.compile(constants['mzRegex'])
        
        def newParser(description):
            return float(mzRegCompiled.search(description).group())
        spectrumDescriptionToMZ = newParser
        
    if 'scanRegex' in constants:
        global spectrumDescriptionToScanNumber
        
        scanRegCompiled = re.compile(constants['scanRegex'])
        
        def newParser(description):
            return float(scanRegCompiled.search(description).group())
        spectrumDescriptionToScanNumber = newParser 
    
    #if 'featureTolerance' in constants:
        #global peakFindTolerance
        #peakFindTolerance = constants['featureTolerance']
    
    if 'signalNoiseThreshold' in constants:
        global signalToNoiseThreshold
        signalToNoiseThreshold = constants['signalNoiseThreshold']
        
    
    
    
    
def getAcqPoints(datafile, resultFile):
    data = mzFile(datafile)
    scans = data.scan_info(0, 999999)
    ms2toms1 = {}
    ms1 = scans[0][2]
    ms2s = []
    assert scans[0][3] == 'MS1'
    for scan in scans:
        if scan[3] == 'MS1':
            for ms2 in ms2s:
                ms2toms1[ms2] = ms1
            ms1 = scan[2]
            ms2s = []
        elif scan[3] == 'MS2':
            ms2s.append(scan[2])
        else:
            raise Exception, "Unidentified scan type of %s" % scan[3]
    for ms2 in ms2s:
        ms2toms1[ms2] = ms1
        
    acqPoints = []
    for result in resultFile:      
        mz = spectrumDescriptionToMZ(result['Spectrum Sescription'])
        scan = spectrumDescriptionToScanNumber(result['Spectrum Description'])
        scan = data.timeForScan(ms2toms1[scan])
        acqPoints.append((mz, scan))    
    
    return acqPoints
        
        

def binByFullFeature(datafile, featureDB, results):
    data = mzFile(datafile)
    
    scans = data.scan_info(0, 999999)
    ms2toms1 = {}
    ms1 = None
    ms2s = []
    # MS2s are dropped until the first MS1.
    for scan in scans:
        if scan[3] == 'MS1':
            for ms2 in ms2s:
                ms2toms1[ms2] = ms1
            ms1 = scan[2]
            ms2s = []
        elif scan[3] == 'MS2':
            if ms1 != None:
                ms2s.append(scan[2])
        else:
            raise Exception, "Unidentified scan type of %s" % scan[3]
    for ms2 in ms2s:
        ms2toms1[ms2] = ms1   
           
    matchesToSplits = 0
    matchesToUnsplit = 0
    featureItems = defaultdict(list)
    edgeItems = defaultdict(list)
    inexplicableItems = []
    for result in results:
        #mz = spectrumDescriptionToMZ(result['Spectrum Description'])
        #scan = spectrumDescriptionToScanNumber(result['Spectrum Description']) 
        mz = mzFromPSM(result)
        scan = scanFromPSM(result)
        charge = int(result['Charge'])
        try:
            scan = ms2toms1[scan]
        except:
            continue
        
        features = [(i, x) for i, x in featureDB.mz_range(mz - 0.01, mz + 0.01)
                    if x.containsPoint(mz, scan, charge)]
        if features:
            index, feature = min(features, key = lambda x: abs(x[1].mz - mz))
            scans = min(feature.scans), max(feature.scans)
            intensity = feature.c12Intensity()
            kurtosis = feature.kurtosis
            skew = feature.skewness
            featureItems[index].append((result, scans, intensity, kurtosis, skew))
        else:
            features = [(i, x) for i, x in featureDB.mz_range(mz - 1, mz + 1)
                        if x.bordersPoint(mz, scan, charge)]
            if features:
                index, feature = min(features, key = lambda x: abs(x[1].mz - mz))
                edge = feature.bordersPoint(mz, scan, charge)
                scans = min(feature.scans), max(feature.scans)
                intensity = feature.c12Intensity()
                kurtosis = feature.kurtosis
                skew = feature.skewness                
                edgeItems[index].append((result, edge, scans, intensity, kurtosis, skew))
            else:
                inexplicableItems.append(result)
                
        
        
    groupedResults = []
    overFitCount = 0
    for feature, results in featureItems.items():
        try:
            pep = results[0][0]['Peptide Sequence']
            if not all([x['Peptide Sequence'] == pep for x, s, i, k, sk in results]):
                overFitCount += 1            
        except KeyError:
            pep = results[0][0]['Annotated Sequence']
            if not all([x['Annotated Sequence'] == pep for x, s, i, k, sk in results]):
                overFitCount += 1               

        
        for result, scans, intensity, kurtosis, skew in results:
            result['Feature'] = feature
            result['feature error'] = '-'
            result['feature start scan'] = scans[0]
            result['feature end scan'] = scans[1]
            result['feature start time'] = data.timeForScan(scans[0])  if scans[0] else '-'
            result['feature end time'] = data.timeForScan(scans[1])  if scans[1] else '-'
            result['feature intensity'] = intensity
            result['feature kurtosis'] = kurtosis
            result['feature skewness'] = skew
            groupedResults.append(result)
    for feature, resultEdges in edgeItems.items():
        for result, edge, scans, intensity, kurtosis, skew in resultEdges:
            result['Feature'] = '-'
            result['feature error'] = str(feature) + " " + edge
            result['feature start scan'] = scans[0]
            result['feature end scan'] = scans[1]
            result['feature start time'] = data.timeForScan(scans[0]) if scans[0] else '-'
            result['feature end time'] = data.timeForScan(scans[1]) if scans[1] else '-'    
            result['feature intensity'] = intensity
            result['feature kurtosis'] = kurtosis
            result['feature skewness'] = skew
            groupedResults.append(result)
    for result in inexplicableItems:
        result['Feature'] = '-'
        result['feature error'] = 'Feature not found'
        result['feature start scan'] = '-'
        result['feature end scan'] = '-'
        result['feature start time'] = '-'
        result['feature end time'] = '-'
        result['feature intensity'] = '-'
        result['feature kurtosis'] = '-'
        result['feature skewness'] = '-'
        groupedResults.append(result)

    data.close()
    return groupedResults
        
        
  
  
def runSearch(datafile, resultFiles):
    assert datafile.lower().endswith('.raw'), "Only .raw files are currently supported."
  
    features = detect_features(datafile)  
  

def feature_analysis(datafile, resultFiles,
                     featureFile = None,
                     tolerance = None,
                     mzRegex = None, scanRegex = None,
                     **constants):
    
    """
    Performs feature-detection analysis on the given .RAW file and PSM
    reports. The output files group the given PSMs by feature, with the
    addition of source feature extent and intensity information.
    
    """


    import os

    if mzRegex:
        import re
        global spectrumDescriptionToMZ
        
        mzRegCompiled = re.compile(mzRegex)
        
        def newParser(description):
            return float(mzRegCompiled.search(description).group())
        spectrumDescriptionToMZ = newParser
    
    if scanRegex:
        import re
        global spectrumDescriptionToScanNumber
        
        scanRegCompiled = re.compile(scanRegex)
        
        def newParser(description):
            return int(scanRegCompiled.search(description).group())
        spectrumDescriptionToScanNumber = newParser
        
    #if tolerance:
        #global peakFindTolerance
        #peakFindTolerance = tolerance
    
    #if signalNoise:
        #global signalToNoiseThreshold
        #signalToNoiseThreshold = signalNoise
        
        


    assert os.path.exists(datafile), "%s not found!" % datafile
    for resultfile in resultFiles:
        assert os.path.exists(resultfile), "%s not found!" % resultfile
    assert datafile.lower().endswith('.raw'), "Only .raw files are currently supported."
    
    if featureFile:
        assert os.path.exists(featureFile), "Specified feature data file %s not found!" % featureFile
    else:
        featureFile = detect_features(datafile, tolerance = tolerance, **constants)
    features = FeatureInterface(featureFile)
    
    outputfiles = []
    if resultFiles:
        print resultFiles
        print "Categorizing search results by file."
        for resultfile in resultFiles:
            resultfile = os.path.abspath(resultfile)
            inputResults = mzReport.reader(resultfile)
            outputfile = '.'.join(resultfile.split('.')[:-1] + ['featureDetect', 'xlsx']) 
            outputfiles.append(outputfile)
            
            
            resultsByFeature = binByFullFeature(datafile, features, inputResults)
            
            output = mzReport.writer(outputfile,
                                     columns = inputResults.columns + ['Feature',
                                                                       'feature error',
                                                                       'feature start scan',
                                                                       'feature end scan',
                                                                       'feature start time',
                                                                       'feature end time',
                                                                       'feature intensity',
                                                                       'feature kurtosis',
                                                                       'feature skewness'])
            
            for result in resultsByFeature:
                output.write(result)
            
            output.close()
            
            print "Output saved to %s ." % outputfile
    else:
        print "No PSM data given; skipping annotation step."
        
    return featureFile, outputfiles



def dataReaderProc(datafile, que, scanNumbers):
    try:
        data = mzFile(datafile)
        
        for scanNum in scanNumbers:
            scan = data.scan(scanNum, centroid = True)
            que.put((scanNum, scan), block = True)
    
        que.put('done')
        data.close()
    except Exception as err:
        import traceback
        print "READ THREAD ERROR."
        traceback.print_exc()
        print '------------------'
        raise err

    


    
    
def detect_features(datafile, **constants):
    """
    Runs the feature detection algorithm on the target data file (currently,
    only Thermo .RAW is supported.)  Returns the path to the feature data
    file.
    
    Optional arguments:
    - tolerance (default 10): MZ tolerance in parts-per-million for all determinations
    of peak identity.  Should usually correspond to the mass precision of the
    source instrument.
    - force (default False): If True, feature detection is run even if a
    feature data file already exists for the target data.
    """
    
    
    if 'outputfile' in constants:
        featurefile = constants['outputfile']
    else:
        featurefile = datafile + '.features'
    
    if 'tolerance' in constants and constants['tolerance']:
        global tolerance
        tolerance = constants['tolerance']
        if tolerance < 1:
            print "\n\n\nWARNING- tolerance value for SILAC analysis should now be in PPM!\n\n\n"
    else:
        tolerance = 10
        
    if 'partial' in constants:
        # This is primarily for testing purposes only.
        scanrange = constants['partial']
    else:
        scanrange = None
        
    if 'force' in constants:
        force = constants['force']
    else:
        force = False
        
    if 'whitelist_psms' in constants:
        whitelist_mzs = constants['whitelist_psms']
        featurefile = datafile + '.partial%s.features' % (str(hash(frozenset(whitelist_mzs)))[:5])
    else:
        whitelist_mzs = None
        
    if 'peak_picking_params' in constants:
        peak_pick_params = constants['peak_picking_params']
    elif 'tolerance' in constants and constants['tolerance']:
        peak_pick_params = {'tolerance':constants['tolerance']}
    else:
        peak_pick_params = {'tolerance' : 10}
    
    if os.path.exists(featurefile) and not force:
        vprint("Feature data file already exists: %s" % featurefile)
        return featurefile
    
    setGlobals(constants)

    
    times = []
    times.append(time.clock())
    data = mzFile(datafile)
    
    times.append(time.clock())
    vprint("Opened data file; getting isotopes...")

    scaninfo = [x for x in data.scan_info(0, 99999999) if x[3] == 'MS1']
    rtLookup = dict([(x[2], x[0]) for x in scaninfo])
    scaninfo = [x[2] for x in scaninfo]
    
    if scanrange:
        scaninfo = [x for x in scaninfo if scanrange[0] < x < scanrange[1]]

    data.close()
    
    que = multiprocessing.Queue(maxsize = 20)
    reader = multiprocessing.Process(target = dataReaderProc,
                                     args = (datafile, que, scaninfo))
    reader.start()
    
    isotopeData = deque()
    thing = que.get(block = True)
    bar = 0
    while thing != 'done':
        scanNum, scan = thing
        foo = time.clock()
        isotopeData.append((scanNum, peak_pick_PPM(scan, **peak_pick_params)[0]))
        bar += time.clock() - foo
        
        thing = que.get(block = True)
        
        if verbose_mode and len(isotopeData) % 100 == 0:
            print len(isotopeData) # Shielded by explicit verbose_mode check.
    
    reader.join()
    # Could just discard the un-feature'd peaks immediately.
    vprint("Isotopic features acquired; finding features over time...")

    times.append(time.clock())

    ms1ToIndex = {}
    indexToMS1 = {}
    for index, scanNum in enumerate(scaninfo):
        ms1ToIndex[scanNum] = index
        indexToMS1[index] = scanNum

            
    isotopesByChargePoint = defaultdict(lambda: defaultdict(lambda: ProximityIndexedSequence([], lambda x: x[0][0])))
    allIsotopes = []
    for scanNum, isotopesByCharge in isotopeData:
        scanIndex = ms1ToIndex[scanNum]
        for charge, isotopes in isotopesByCharge.items():
            for isoSeq in isotopes:
                isotopesByChargePoint[charge][scanIndex].add(isoSeq)
                allIsotopes.append((isoSeq, scanIndex, charge))
    
    del isotopeData

    for scanlookup in isotopesByChargePoint.values():
        for proxseq in scanlookup.values():
            proxseq.rebalance()
            

    if whitelist_mzs:
        vprint("Screening out irrelevant MZs; starting with %s..." % len(allIsotopes))
        allIsotopes.sort(key = lambda x: x[0][0][0])
        whitelist_mzs = sorted(list(set([round(x, 2) for x in whitelist_mzs])))
        isoAcc = []
        whitemz = whitelist_mzs.pop()
        while allIsotopes:
            iso = allIsotopes.pop()
            mz = iso[0][0][0]
            while whitelist_mzs and whitemz - mz > whitelist_tol:
                whitemz = whitelist_mzs.pop()
            if abs(whitemz - mz) < whitelist_tol:
                isoAcc.append(iso)
        
        allIsotopes = isoAcc
        vprint("...%s remain." % len(allIsotopes))
    
    
    
    allIsotopes.sort(key = lambda x: x[0][0][1])
    

    
    times.append(time.clock())    
    
    seenIsotopes = set()
    # Can assume isotopic sequences are unique because floats.
    # (But it may not be a valid assumption, because detectors
    # and floating point approximations!)
    
    featureList = []
    while allIsotopes:
        highIso, highScan, highChg = allIsotopes.pop()
        if tuple(highIso) in seenIsotopes:
            continue
        
        centerIndex, (centerMZ, _) = max(enumerate(highIso), 
                                         key = lambda x: x[1][1])
        
        newFeature = [[highScan, highIso]]
        curScan = highScan
        continuing = True
        lastSeen = rtLookup[indexToMS1[curScan]]
        while continuing: # Trailing the feature backwards.
            curScan -= 1
            try:
                curRT = rtLookup[indexToMS1[curScan]]
            except KeyError:
                assert curScan < max(indexToMS1.keys())
                break
            
            
            scanSeqs = isotopesByChargePoint[highChg][curScan].returnRange(centerMZ - 2, centerMZ + 1.5)
            scanSeqs.sort(key = lambda x: x[centerIndex][1], reverse = True)
            
            found = False
            for iso in scanSeqs: # These are known to have centerMZ in common.
                # The indexes between iso and highIso may not be equivalent
                # if there's sub-C12 peak(s) in either.  For a first draft
                # this can be considered a feature, since C12s should be
                # consistent throughout features, but in some cases like
                # single-scan-dropouts of the C12 this is insufficient
                # and such discrepancies should be accounted for.
                
                if (inPPM(tolerance, iso[0][0], highIso[0][0])
                    and inPPM(tolerance, iso[1][0], highIso[1][0])
                    and tuple(iso) not in seenIsotopes):
                    newFeature.append([curScan, iso])
                    found = True
                    break # From "for iso in scanSeqs"                    
            
            if found:
                lastSeen = curRT
            elif abs(curRT - lastSeen) > dropoutTimeTolerance:
                continuing = False
        
        curScan = highScan
        continuing = True
        lastSeen = rtLookup[indexToMS1[curScan]]
        while continuing: # Trailing the feature forwards; mostly repeat code.
            curScan += 1
            try:
                curRT = rtLookup[indexToMS1[curScan]]
            except KeyError:
                assert curScan > max(indexToMS1.keys())
                break
            
            
            scanSeqs = isotopesByChargePoint[highChg][curScan].returnRange(centerMZ - 2, centerMZ + 1.5)
            scanSeqs.sort(key = lambda x: x[centerIndex][1], reverse = True)
            
            found = False
            for iso in scanSeqs: # These are known to have centerMZ in common.
                # Ditto.
                               
                if (inPPM(tolerance, iso[0][0], highIso[0][0])
                    and inPPM(tolerance, iso[1][0], highIso[1][0])
                    and tuple(iso) not in seenIsotopes):
                    newFeature.append([curScan, iso])
                    found = True
                    break # From "for iso in scanSeqs"                    

            if found:
                lastSeen = curRT
            elif abs(curRT - lastSeen) > dropoutTimeTolerance:
                continuing = False

        if len(newFeature) > 1:
            featureList.append((highChg, newFeature))
        
        for _, iso in newFeature:
            seenIsotopes.add(tuple(iso))
    times.append(time.clock())
    
    for chg, feature in featureList:
        for stage in feature:
            stage[0] = indexToMS1[stage[0]]
            
    class idLookup():
        def __getitem__(self, thing):
            return thing
    lookup = idLookup()

    if scanrange:
        featurefile = datafile + ('%s-%s.features' % scanrange)

    featureObjects = []
    for chg, feature in featureList:
        newfeature = Feature()
        for scan, envelope in feature:
            newfeature.add(envelope, scan, chg)
        
        newfeature.calculate_bounds(lookup)
        
        #newfeature.prepareBoxes(lookup)
        #newfeature.prepareBoxes() # It's entirely different, for some reason?

        #test = Feature()
        #for scan, envelope in feature:
            #test.add(envelope, scan, chg)
        #test.calculate_bounds(lookup)
        
        #assert test.mz == newfeature.mz and test.charge == newfeature.charge
        
        featureObjects.append(newfeature)
    save_feature_database(featureObjects, featurefile)
    
    vprint("Saved feature file.")
    times.append(time.clock())
    
    return featurefile



# Old-naming-convention alias, for legacy purposes.
detectFeatures = detect_features


    
    
    
# RUNNING FEATURE DETECTION BY USING THIS FILE AS __MAIN__ DOESN'T WORK.