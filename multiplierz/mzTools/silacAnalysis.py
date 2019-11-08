from multiplierz.mzTools.featureDetector import detectFeatures, Feature
from multiplierz.mzAPI import mzFile
from multiplierz.mzReport import reader, writer, mzSpreadsheet
from multiplierz import myData, vprint
from collections import defaultdict
import pickle
import matplotlib.pyplot as pyt
from multiplierz.mass_biochem import mz, unimod
#try:
    #from scipy.stats import pearsonr
#except ImportError:
    #print "Could not import scipy.stats.pearsonr!"
    #def pearsonr(x, y):
        #return float('NaN'), float('NaN')
def pearsonr(x, y):
    return float('NaN'), float('NaN')        
from numpy import average, median, floor, ceil
import os
import re
from multiplierz.internalAlgorithms import ProximityIndexedSequence, inPPM, select
from multiplierz.mzTools.featureUtilities import FeatureInterface
from multiplierz.mgf import standard_title_parse
import warnings

__all__ = ['SILAC2Plex', 'SILAC3Plex']

    


peakFindTolPPM = 10
XICTol = 0.008    
    

def spectrumDescriptionToMZ(description):
    try:
        return float(standard_title_parse(description)['mz'])
    except ValueError:
        try:
            return float(description.split('-|')[0].split('-')[-1])
        except ValueError:
            return float(description.split("|")[1])

def spectrumDescriptionToScanNumber(description):
    try:
        return int(standard_title_parse(description)['scan'])
    except ValueError:
        return int(description.split('.')[1])
        
mediumShifts = None
heavyShifts = None

protonMass = 1.00727647

# Precise shift from isotopes.
isotopeDiffs = {'C':1.003355,
                'N':0.997034}

allTags = None
heavyR = None
heavyK = None
mediumR = None
mediumK = None


def unzip(thing, default = 2): 
    if thing:
        return [list(x) for x in zip(*thing)]
    else:
        return [[]]*default





class FeatureMemo(object):
    """
    Wrappers within wrappers!
    
    Obtains bins that memoize and easily retrieve MZ-ranges of features.
    """
    def __init__(self, featurefile):
        self.features = FeatureInterface(featurefile)
        self.bins = defaultdict(list)
    
    def get_bin(self, down, up):
        key = down, up
        return self.features.mz_range(*key)
        
    def mzs_around(self, mz, width):
        low, high = mz - (width/2), mz + (width/2)
        return self.get_bin(low, high)
            
        
        
        

def parseModifications(modstring):
    labels = [x.strip() for x in modstring.split(";") if x.strip()]
    
    tags = []
    mods = []
    for label in labels:
        loc, modtype = label.split()
        # Throw out modification location?  Does that matter?
        
        if modtype in allTags:
            tags.append(modtype)
        else:
            mods.append(loc + modtype)

    return (tags, mods)

def findDoubles(data):
    pepMods = defaultdict(set)
    
    def unsilacMods(modstring):
        if modstring:
            mods = modstring.split(';')
            return frozenset([mod.strip() for mod in mods if not any([tag in mod for tag in allTags])])
        else:
            return frozenset()
    
    for psm in list(data.values()):
        charge = psm['Charge']
        sequence = psm['Peptide Sequence']
        mods = unsilacMods(psm['Variable Modifications'])
        
        # Only peptides that can be SILAC labelled (have K or R) are included in the analysis.
        if 'K' in sequence or 'R' in sequence:
            pepMods[sequence, mods, charge].add(psm['Spectrum Description'])

    doubles = []
    for cosilacSpecDescs in list(pepMods.values()):
        cosilacPSMs = [data[x] for x in cosilacSpecDescs]
        
        lights = [x for x in cosilacPSMs if not any(parseModifications(x['Variable Modifications'])[0])]
        heavies = [x for x in cosilacPSMs if any(parseModifications(x['Variable Modifications'])[0])]
        
        assert len(lights) + len(heavies) == len(cosilacPSMs)
        
        # Filtering irrelevant label mods; essentially for testing with medium-as-

        lights = [x for x in lights if 
                  not any([label in x['Variable Modifications'] for label in allTags])]
        
        # Incompletely labelled peptides are not counted in the analysis.
        heavies = [x for x in heavies if
                   len([m for m in x['Variable Modifications'].split(';') if any([p in m for p in allTags])])
                   == x['Peptide Sequence'].count('K') + x['Peptide Sequence'].count('R')]
        
        doubles.append(([psm['Spectrum Description'] for psm in lights],
                        [psm['Spectrum Description'] for psm in heavies]))
        
    return doubles


def findDoublesAdapter(data):
    doubles = findDoubles(data)
    
    newdoubles = []
    for lights, heavies in doubles:
        newdoubles.append((list(zip(lights, [None]*len(lights), [None]*len(lights))),
                          [],
                          list(zip(heavies, [None]*len(heavies), [None]*len(heavies)))))
        
    return newdoubles
        
        
def findTriples(data):
    pepMods = defaultdict(set)
    
    def unsilacMods(modstring):
        if modstring:
            mods = modstring.split(';')
            return frozenset([mod.strip() for mod in mods if not any([tag in mod for tag in allTags])])
        else:
            return frozenset()
    
    for psm in list(data.values()):
        charge = psm['Charge']
        sequence = psm['Peptide Sequence']
        mods = unsilacMods(psm['Variable Modifications'])
        
        # Filtering out unlabelable peptides.
        if 'K' in sequence or 'R' in sequence:
            pepMods[sequence, mods, charge].add(psm['Spectrum Description'])
        
    triples = []
    for cosilacSpecDescs in list(pepMods.values()):
        cosilacPSMs = [data[spectrumDescriptionToScanNumber(x)] for x in cosilacSpecDescs]
        
        lights = []
        mediums = []
        heavies = []        
        for psm in cosilacPSMs:
            if psm['Variable Modifications']:
                varMods = [x.strip() for x in psm['Variable Modifications'].split(';') if x]
            else:
                varMods = []
            medRMods = [x for x in varMods if x[0] == 'R' and mediumR in x and heavyR not in x and heavyK not in x]
            medKMods = [x for x in varMods if mediumK in x and heavyR not in x and heavyK not in x]
            heavyKMods = [x for x in varMods if x[0] == 'K' and heavyK in x]
            heavyRMods = [x for x in varMods if x[0] == 'R' and heavyR in x]
            
            if (medRMods or medKMods) and (heavyKMods or heavyRMods):
                raise Exception
            
            # Filtering out incompletely labelled peptides.
            seq = psm['Peptide Sequence']
            if (medRMods and len(medRMods) != seq.count('R') or
                medKMods and len(medKMods) != seq.count('K') or
                heavyRMods and len(heavyRMods) != seq.count('R') or
                heavyKMods and len(heavyKMods) != seq.count('K')):
                continue
            
            if medRMods or medKMods: mediums.append(psm)
            elif heavyKMods or heavyRMods: heavies.append(psm)
            else: lights.append(psm)
        
        assert len(lights) + len(mediums) + len(heavies) == len(cosilacPSMs), "Bad tag lists?"
        
        triples.append(([psm['Spectrum Description'] for psm in lights],
                        [psm['Spectrum Description'] for psm in mediums],
                        [psm['Spectrum Description'] for psm in heavies]))
    
    return triples

def findTriplesAdapter(data):
    triples = findTriples(data)
    newtriples = []
    for lights, mediums, heavies in triples:
        newtriples.append((list(zip(lights, [None] * len(lights), [None]* len(lights))),
                          list(zip(mediums, [None] * len(mediums), [None]* len(mediums))),
                          list(zip(heavies, [None] * len(heavies), [None]* len(heavies)))))
        
    return newtriples



def fetchData(peptideFile):
    #data = mzSpreadsheet.get_xl_sheet(peptideFile, sheet_name = 'Data')
    #mzSpreadsheet._stop_excel() # Just to be safe.
    data = list(reader(peptideFile, autotypecast = False))
    return data

def indexTable(peptideData, indexTitle, uniqueIndex = True):
    if isinstance(peptideData[0], dict):
        data = {}
        for psm in peptideData:
            data[psm[indexTitle]] = psm
            
        return data
    
    header = peptideData[0]
    titleCol = header.index(indexTitle)
    
    if uniqueIndex:
        data = {}
        for line in peptideData[1:]:
            entry = {}
            for index, item in enumerate(header):
                entry[item] = line[index] if line[index] else ""
            data[line[titleCol]] = entry
    else:
        # This suddenly becomes dubious business.
        class hashabledict(dict):
            def __hash__(self):
                return hash(frozenset(sorted(self.items())))
            def __eq__(self, other):
                return self.__hash__() == other.__hash__()
        
        data = defaultdict(set)
        for line in peptideData[1:]:
            entry = hashabledict()
            for index, item in enumerate(header):
                entry[item] = line[index] if line[index] else ""
            data[line[titleCol]].add(entry)

    return data


def getPSMIntensity(resultIndex, featureList, ms2toms1, specDesc, shift = None):
    #scan = ms2toms1[int(specDesc.split('.')[1])]
    scan = ms2toms1[spectrumDescriptionToScanNumber(specDesc)]
    chg = float(resultIndex[specDesc]['Charge'])
    #mz = float(specDesc.split('|')[1]) + (shift/float(chg))
    
    try:
        mz = spectrumDescriptionToMZ(specDesc) + (shift/float(chg))
    except IndexError:
        mz = float(resultIndex[specDesc]['Experimental mz']) + (shift/float(chg))
    
    try:
        match = next((x for x in featureList if x.containsPoint(mz, scan, chg)))
    except StopIteration:
        return 0
    
    totalInt = 0
    for index, region in match.regions:
        for mz, intensity, _ in region:
            if mz >= match.mz:
                totalInt += intensity
    
    return totalInt
    


def sumFeatures(features):
    summation = Feature()
    for feature in features:
        charge = feature.charge
        for index, region in feature.regions:
            summation.add(region, index, charge)
    return summation

        
def getRatio(peaks):
    if not peaks: return []
    peaks.sort(key = lambda x: x[0])
    intensities = [x[1] for x in peaks]
    highest = float(max(intensities))
    return [x / highest for x in intensities]

def featureSimilarities(firstFeatures, secondFeatures):
    if not (firstFeatures and secondFeatures):
        raise Exception
    
    first = sumFeatures(firstFeatures)
    second = sumFeatures(secondFeatures)
    
    allRegions = sorted(list(set(unzip(first.regions)[0] + unzip(second.regions)[0])))
    firstRegions = dict(first.regions)
    secondRegions = dict(second.regions)

    ratioScores = []
    coIntensities = []
    
    for index in allRegions:
        try: firstPeaks = firstRegions[index]
        except KeyError: firstPeaks = []
        try: secondPeaks = secondRegions[index]
        except KeyError: secondPeaks = []
        
        firstInt = sum([x[1] for x in firstPeaks])
        secondInt = sum([x[1] for x in secondPeaks])
        coIntensities.append((firstInt, secondInt))
        
        firstRatio = getRatio(firstPeaks)
        secondRatio = getRatio(secondPeaks)
        
        try: firstMost = firstRatio.index(max(firstRatio))
        except ValueError: firstMost = 0
        try: secondMost = secondRatio.index(max(secondRatio))
        except ValueError: secondMost = 0
        
        difference = 0
        for offset in range(-10, 10, 1):
            first = firstMost + offset
            second = secondMost + offset
            haveFirst = 0 <= first < len(firstRatio)
            haveSecond = 0 <= second < len(secondRatio)
            
            if not (haveFirst or haveSecond): continue
            elif haveFirst and not haveSecond:
                difference += firstRatio[first]**2
            elif haveSecond and not haveFirst:
                difference += secondRatio[second]**2
            else:
                difference += (firstRatio[first] - secondRatio[second])**2
        ratioScores.append((difference, firstInt + secondInt))
                
    correlation = pearsonr(*unzip(coIntensities))[0]
    finalRatioScore = average(unzip(ratioScores)[0], weights = unzip(ratioScores)[1])
    
    
    return correlation, finalRatioScore

def tripleFeatureSimilarities(firstFeatures, secondFeatures, thirdFeatures):
    if firstFeatures and secondFeatures:
        firstSecondCorrelation, firstSecondRatioMatch = featureSimilarities(firstFeatures,
                                                                            secondFeatures)
    else:
        firstSecondCorrelation, firstSecondRatioMatch = '-', '-'
    
    if firstFeatures and thirdFeatures:
        firstThirdCorrelation, firstThirdRatioMatch = featureSimilarities(firstFeatures,
                                                                          thirdFeatures)
    else:
        firstThirdCorrelation, firstThirdRatioMatch = '-', '-'
    
    return (firstSecondCorrelation, firstThirdCorrelation,
            firstSecondRatioMatch, firstThirdRatioMatch)

def overlapIntensity(data, lightFeatures, heavyFeatures, shift):
    light = sumFeatures(lightFeatures)
    heavy = sumFeatures(heavyFeatures)
    
    #allRegions = sorted(list(set(unzip(light.regions)[0] + unzip(heavy.regions)[0])))
    lightRegions = dict(light.regions)
    heavyRegions = dict(heavy.regions)
    allRegions = list(set(list(lightRegions.keys()) + list(heavyRegions.keys())))
    try:
        assert sorted(allRegions) == sorted(list(set(unzip(light.regions)[0] + unzip(heavy.regions)[0])))
    except IndexError:
        pass
    lightIntensity = 0
    heavyIntensity = 0

    
    if lightRegions and heavyRegions:
        note = ''
        for index in allRegions:
            try: lightPeaks = lightRegions[index]
            except KeyError: continue
            try: heavyPeaks = heavyRegions[index]
            except KeyError: continue
            # If there's not both features in the scan, its irrelevant
            # to the overlap.
            
            for lightPeak in lightPeaks:
                try:
                    heavyPeak = min(heavyPeaks, 
                                    key = lambda x: abs((lightPeak[0]+shift) - x[0]))
                    
                    if abs((lightPeak[0]+shift) - heavyPeak[0]) > 0.1: continue
                    else:
                        lightIntensity += lightPeak[1]
                        heavyIntensity += heavyPeak[1]
                except ValueError: continue
    
    elif lightRegions and not heavyRegions:
        note = 'Imputed missing heavy feature intensity from noise'
        
        for index in lightRegions:
            scan = data.lscan(index)
            for lightPeak in lightRegions[index]:
                assert lightPeak in [x[:2] for x in scan]
                
                lightIntensity += lightPeak[1]
                
                heavymz = lightPeak[0] + shift
                nearHeavy = min(scan, key = lambda x: abs(x[0] - heavymz))
                #if abs(nearHeavy[0] - heavymz) < peakFindTolerance:
                if inPPM(peakFindTolPPM, nearHeavy[0], heavymz):
                    heavyIntensity += nearHeavy[1]
                else:
                    heavyIntensity += nearHeavy[2] # Noise.
    
    elif heavyRegions and not lightRegions:
        note = 'Imputed missing light feature intensity from noise'
        
        for index in heavyRegions:
            scan = data.lscan(index)
            for heavyPeak in heavyRegions[index]:
                assert heavyPeak in [x[:2] for x in scan]
                
                heavyIntensity += heavyPeak[1]
                
                lightmz = heavyPeak[0] - shift
                nearLight = min(scan, key = lambda x: abs(x[0] - lightmz))
                #if abs(nearLight[0] - lightmz) < peakFindTolerance:
                if inPPM(peakFindTolPPM, nearLight[0], lightmz):
                    lightIntensity += nearLight[1]
                else:
                    lightIntensity += nearLight[2] # Noise.
                
    
    return lightIntensity, heavyIntensity, note
            
            
def interpolatePoints(firstPoints, secondPoints):
    """
    Its a good idea to look for complimentary SILAC features at each other's
    times, even in cases where the PSMs were time-offset, since the PSMs are
    more likely to occur outside of their features' times (perversely!).
    """
    
    firstScans = select(1, firstPoints)
    secondScans = select(1, secondPoints)
    
    onlyFirstScans = set(firstScans) - set(secondScans)
    onlySecondScans = set(secondScans) - set(firstScans)
    
    firstMZ = median(select(0, firstPoints))
    secondMZ = median(select(0, secondPoints))
    firstChg = median(select(2, firstPoints))
    secondChg = median(select(2, secondPoints))
    
    firstPoints += [(firstMZ, scan, firstChg) for scan in onlySecondScans]
    secondPoints += [(secondMZ, scan, secondChg) for scan in onlyFirstScans]
    
    return firstPoints, secondPoints
    
    
    
def interpolateTriplePoints(firstPoints, secondPoints, thirdPoints):
    firstScans = select(1, firstPoints)
    secondScans = select(1, secondPoints)
    thirdScans = select(1, thirdPoints)
    
    firstMZ = median(select(0, firstPoints))
    secondMZ = median(select(0, secondPoints))
    thirdMZ = median(select(0, thirdPoints))
    firstChg = median(select(2, firstPoints))
    secondChg = median(select(2, secondPoints))
    thirdChg = median(select(2, thirdPoints))
    
    # There's presumably some better form of this set arithmetic.
    
    onlyFirstScans = set(firstScans) - (set(thirdScans) | set(secondScans))
    onlySecondScans = set(secondScans) - (set(firstScans) | set(thirdScans))
    onlyThirdScans = set(thirdScans) - (set(firstScans) | set(secondScans))
    
    firstPoints += [(firstMZ, scan, firstChg) for scan in onlySecondScans | onlyThirdScans]
    secondPoints += [(secondMZ, scan, secondChg) for scan in onlyFirstScans | onlyThirdScans]
    thirdPoints += [(secondMZ, scan, secondChg) for scan in onlyFirstScans | onlySecondScans] 
    
    return firstPoints, secondPoints, thirdPoints
            

def featureIntensities(points, pointTol, features):
    foundFeatureIndices = []
    foundFeatures = set([])
    totalIntensity = 0
    
    #assert isinstance(features, ProximityIndexedSequenceAgain)
    #assert isinstance(features, FeatureInterface)
    
    for pt in points:
        try:
            #allMatches = [(i, x) for (i, x) in enumerate(features) if x.containsPoint(*pt)]
            #index, match = max(allMatches, key = lambda x: x[1].totalIntensity())
            #allMatches = features.returnRange(pt[0] - 1, pt[0] + 1)
            #allMatches = features.mz_range(pt[0] - 1, pt[0] + 1)
            
            daTol = (pt[0]/1000000) * pointTol
            
            allMatches = features.mzs_around(pt[0], width = daTol) # FeatureMemo method.
            allMatches = [(i, x) for (i, x) in allMatches if x.containsPoint(*pt)]
            index, match = max(allMatches, key = lambda x: x[1].totalIntensity())
                
        except ValueError:
            continue
        #if match in foundFeatures: continue
        if str(index) in foundFeatureIndices:
            continue
        
        foundFeatures.add(match)
        foundFeatureIndices.append(str(index))
        for index, region in match.regions:
            totalIntensity += sum([x[1] for x in region if abs(x[0] - pt[0]) < 0.1])
        #break
        
            
    return totalIntensity, foundFeatures, foundFeatureIndices


def overlapIntensityTriple(data,
                           lights, mediums, heavies,
                           lightMediumShift, lightHeavyShift):
    light = sumFeatures(lights)
    medium = sumFeatures(mediums)
    heavy = sumFeatures(heavies)

    allRegions = sorted(list(set(unzip(light.regions)[0] +
                                 unzip(medium.regions)[0] +
                                 unzip(heavy.regions)[0])))
    lightRegions = dict(light.regions)
    mediumRegions = dict(medium.regions)
    heavyRegions = dict(heavy.regions)

    def imputeAtPoint(scanNum, peaks):
        scan = data.lscan(scanNum)
        imputedInt = 0
        for pt in peaks:
            mz = pt[0]
            nearPt = min(scan, key = lambda x: abs(x[0] - mz))
            if inPPM(peakFindTolPPM, nearPt[0], mz):
                imputedInt += nearPt[1]
            else:
                imputedInt += nearPt[2] # Noise
        return imputedInt

    # Double dispatch here just to make the map expressions below neater.
    # Currying-enabled languages are so nice!  But Python is not one of them.
    def shiftPeaks(peaks, shift):
        return [(x[0] + shift, x[1]) for x in peaks]
        
    def commonPeaks(first, second, shift):
        firstcommon = []
        secondcommon = []
        for peak in first:
            matches = [x for x in second if inPPM(peakFindTolPPM*5,
                                                  (peak[0] + shift),
                                                  x[0])]
            if matches:
                firstcommon.append(peak)
                secondcommon.append(matches[0])
        return firstcommon, secondcommon

    lightIntensity = 0
    mediumIntensity = 0
    heavyIntensity = 0
    for index in allRegions:
        try: lightPeaks = lightRegions[index]
        except KeyError: lightPeaks = []

        try: mediumPeaks = mediumRegions[index]
        except KeyError: mediumPeaks = []

        try: heavyPeaks = heavyRegions[index]
        except KeyError: heavyPeaks = []
        
        
        # A bunch of crazy business to get a common set of peaks and/or places
        # to look for imputing peaks.
        # If all labelled states are present, the overlap is all peaks shared between
        # the three;
        # If two are present, the overlap is all peaks shared between the two, and the
        # third corresponding set of points is imputed
        # If one is present, the "overlap" is all peaks in the one that's present,
        # and these are all imputed in the others.
        if lightPeaks and mediumPeaks and heavyPeaks:
            lightCommon, mediumCommon = commonPeaks(lightPeaks, mediumPeaks, lightMediumShift)
            lightCommon, heavyCommon = commonPeaks(lightCommon, heavyPeaks, lightHeavyShift)
            mediumCommon, heavyCommon = commonPeaks(mediumCommon, heavyCommon,
                                                    lightHeavyShift - lightMediumShift)
            note = ''
        elif lightPeaks and mediumPeaks:
            lightCommon, mediumCommon = commonPeaks(lightPeaks, mediumPeaks, lightMediumShift)
            heavyCommon = shiftPeaks(lightCommon, lightHeavyShift)
            note = 'Imputed missing heavy feature intensity from noise'
        elif lightPeaks and heavyPeaks:
            lightCommon, heavyCommon = commonPeaks(lightPeaks, heavyPeaks, lightHeavyShift)
            mediumCommon = shiftPeaks(lightCommon, lightMediumShift)
            note = 'Imputed missing medium feature intensity from noise'
        elif mediumPeaks and heavyPeaks:
            mediumCommon, heavyCommon = commonPeaks(mediumPeaks, heavyPeaks,
                                                    lightHeavyShift - lightMediumShift)
            lightCommon = shiftPeaks(heavyCommon, -1*lightHeavyShift)
            note = 'Imputed missing light feature intensity from noise'
        elif lightPeaks:
            lightCommon = lightPeaks
            mediumCommon = shiftPeaks(lightPeaks, lightMediumShift)
            heavyCommon = shiftPeaks(lightPeaks, lightHeavyShift)
            note = 'Imputed missing medium and heavy feature intensities from noise'
        elif mediumPeaks:
            lightCommon = shiftPeaks(mediumPeaks, -lightMediumShift)
            mediumCommon = mediumPeaks
            heavyCommon = shiftPeaks(mediumPeaks, lightHeavyShift - lightMediumShift)
            note = 'Imputed missing light and heavy feature intensities from noise'
        elif heavyPeaks:
            lightCommon = shiftPeaks(heavyPeaks, -lightHeavyShift)
            mediumCommon = shiftPeaks(heavyPeaks, lightMediumShift - lightHeavyShift)
            heavyCommon = heavyPeaks
            note = 'Imputed missing light and medium feature intensities from noise'
        else:
            raise Exception
            
        scan = data.scan(index)
        
        assert len(lightCommon) == len(mediumCommon) == len(heavyCommon)
        
        if lightPeaks:
            lightIntensity += sum([x[1] for x in lightCommon])
        else:
            lightIntensity += imputeAtPoint(index, lightCommon)
        if mediumPeaks:
            mediumIntensity += sum([x[1] for x in mediumCommon])
        else:
            mediumIntensity += imputeAtPoint(index, mediumCommon)
        if heavyPeaks:
            heavyIntensity += sum([x[1] for x in heavyCommon])
        else:
            heavyIntensity += imputeAtPoint(index, heavyCommon)
            

    return lightIntensity, mediumIntensity, heavyIntensity, note



def collidablePSM(psm):
    sequence = psm['Peptide Sequence']
    Rs, Ks = sequence.count('R'), sequence.count('K')
    mediumShift = mediumShifts['R'] * Rs + mediumShifts['K'] * Ks
    heavyShift = heavyShifts['R'] * Rs + heavyShifts['K'] * Ks
    
    expMass = psm['Predicted mr']
    charge = psm['Charge']
    
    
    expectedLights = [(mass / charge) - (charge * protonMass) for mass
                      in [expMass + d for d in range(0, 5)]]
    expectedMediums = [(mass / charge) - (charge * protonMass) for mass
                      in [expMass + mediumShift + d for d in range(0, 5)]]
    expectedHeavies = [(mass / charge) - (charge * protonMass) for mass
                      in [expMass + heavyShift + d for d in range(0, 5)]]
    
    return (any([any([inPPM(peakFindTolPPM, lightX, mediumX) for mediumX in expectedMediums])
                for lightX in expectedLights]) or
            any([any([inPPM(peakFindTolPPM, heavyX, mediumX) for mediumX in expectedMediums])
                 for heavyX in expectedHeavies]))
                 


def reconstructTagFromCollision(featureList, psm, tagKind,
                                lowFeatures, tagFeatures,
                                lowPoints, tagPoints):
    expectedMass = psm['Predicted mr']
    charge = psm['Charge']
    expectedMZ = (expectedMass / charge) + protonMass
    
    sequence = psm['Peptide Sequence']
    Rs, Ks = sequence.count('R'), sequence.count('K')   
    
    assert lowFeatures
    
    c13TagPoints = [(mz + 1/charge, scan, chg) for (mz, scan, chg) in tagPoints]
    
    c13TagFeatures = [x for x in featureList if 
                         any([x.containsPoint(*pt) for pt in c13TagPoints])]    
    
    if not c13TagFeatures:
        raise NotImplementedError
    
    for feature in c13TagFeatures:
        feature.index = str(featureList.index(feature))
        
    scanLookups = []
    for feature in lowFeatures:
        scanLookups.append(dict(feature.regions))
    
    allrecs = []
    #for c13Feature in c13TagFeatures:
    for j in range(0, len(c13TagFeatures)):
        c13Feature = c13TagFeatures[j]
        if tagKind == 'medium':
            expectedC12 = expectedMZ + (mediumShifts['K'] * Ks)/charge + (mediumShifts['R'] * Rs)/charge
        elif tagKind == 'heavy':
            expectedC12 = expectedMZ + (heavyShifts['K'] * Ks)/charge + (heavyShifts['R'] * Rs)/charge
        else:
            raise Exception
        
        #for scan, points in c13Feature.regions:
        for i in range(0, len(c13Feature.regions)):
            scan, points = c13Feature.regions[i]
            
            recoveredC12 = None
            for lowPoints in [x[scan] for x in scanLookups if scan in x]:
                try:
                    recoveredC12 = (x for x in lowPoints if
                                    #abs(x[0] - expectedC12) < peakFindTolerance).next()
                                    inPPM(peakFindTolPPM, expectedC12, x[0]))
                except StopIteration:
                    continue
            
                if recoveredC12:
                    #points.append(recoveredC12)
                    allrecs.append((scan, recoveredC12))
                    c13TagFeatures[j].regions[i] = scan, points + [recoveredC12]
    
    
    for feature in c13TagFeatures: feature.prepareBoxes()
    
    recoveredTagFeatures = [x for x in c13TagFeatures if
                            any([x.containsPoint(*pt) for pt in tagPoints])]
    recoveredIntensity = sum([x.totalIntensity() for x in recoveredTagFeatures])
    recoveredIndices = [x.index for x in recoveredTagFeatures]
    
    return recoveredIntensity, recoveredTagFeatures, recoveredIndices
    
def getTripleRatios(data, resultIndex, featureList, ms2toms1, tagTuples):
    def getShifts(specDesc):
        sequence = resultIndex[specDesc]['Peptide Sequence']
        charge = resultIndex[specDesc]['Charge']
        Rs = sequence.count('R')
        Ks = sequence.count('K')

        lightMediumShift = Rs * mediumShifts['R'] + Ks * mediumShifts['K']
        lightHeavyShift = Rs * heavyShifts['R'] + Ks * heavyShifts['K']
        return lightMediumShift / charge, lightHeavyShift / charge


    def getAcqPoints(psms, predictFrom = None, predictTo = None):
        assert (not predictFrom) or predictFrom in ['light', 'medium', 'heavy']
        assert (not predictTo) or predictTo in ['light', 'medium', 'heavy']
        
        points = []
        for specDesc, _, _ in psms:
            scanNum = spectrumDescriptionToScanNumber(specDesc)
            scan = ms2toms1[scanNum]
            chg = float(resultIndex[scanNum]['Charge'])
            try:
                mz = spectrumDescriptionToMZ(specDesc)
            except IndexError:
                mz = float(resultIndex[scanNum]['Experimental mz'])
            
            if predictFrom and predictTo:
                seq = resultIndex[scanNum]['Peptide Sequence']
                Rs = seq.count('R')
                Ks = seq.count('K')
                
                shift = 0
                if predictFrom == 'heavy':
                    shift -= ((Rs * heavyShifts['R'])/chg + (Ks * heavyShifts['K'])/chg)
                elif predictFrom == 'medium':
                    shift -= ((Rs * mediumShifts['R'])/chg + (Ks * mediumShifts['K'])/chg)
                
                if predictTo == 'heavy':
                    shift += ((Rs * heavyShifts['R'])/chg + (Ks * heavyShifts['K'])/chg)
                elif predictTo == 'medium':
                    shift += ((Rs * mediumShifts['R'])/chg + (Ks * mediumShifts['K'])/chg)
                    
                mz += shift
            
            points.append((mz, scan, chg))
            
        return points
        
    def getDeltaTolerance(psmIndexes):
        scanNums = [spectrumDescriptionToScanNumber(x[0]) for x in psmIndexes]
        avgDelta = average([abs(resultIndex[x]['Delta']) for x in scanNums])
        avgMZ = average([abs(resultIndex[x]['Experimental mz']) for x in scanNums])
        return ((avgDelta*1000000)/avgMZ) + (peakFindTolPPM/2)               

    ratios = []
    for lights, mediums, heavies in tagTuples:
        if not (lights or mediums or heavies): continue
        
        notes = []
        if not heavies:
            notes.append("Extrapolated heavy-tag mz")
        if not mediums:
            notes.append("Extrapolated medium-tag mz")
        if not lights:
            notes.append("Extrapolated untagged mz")        

        lightPoints = []
        mediumPoints = []
        heavyPoints = []

        if lights: lightPoints += getAcqPoints(lights) 
        if mediums: lightPoints += getAcqPoints(mediums, predictTo = 'light',
                                                predictFrom = 'medium')
        if heavies: lightPoints += getAcqPoints(heavies, predictTo = 'light',
                                                predictFrom = 'heavy')
        
        if mediums: mediumPoints += getAcqPoints(mediums)
        if lights: mediumPoints += getAcqPoints(lights, predictTo = 'medium', 
                                                predictFrom = 'light')
        if heavies: mediumPoints += getAcqPoints(heavies, predictTo = 'medium',
                                                 predictFrom = 'heavy')
        
        if heavies: heavyPoints += getAcqPoints(heavies)
        if lights: heavyPoints += getAcqPoints(lights, predictFrom = 'light',
                                               predictTo = 'heavy')
        if mediums: heavyPoints += getAcqPoints(mediums, predictFrom = 'medium',
                                                predictTo = 'heavy')
        
        # Why wasn't this added earlier?
        lightPoints, mediumPoints, heavyPoints = interpolateTriplePoints(lightPoints, 
                                                                         mediumPoints,
                                                                         heavyPoints)
        
        pointTolerance = getDeltaTolerance(lights + mediums + heavies)
        
        (lightIntensity, lightFeatures,
         lightIndices) = featureIntensities(lightPoints, pointTolerance, featureList)
        (mediumIntensity, mediumFeatures,
         mediumIndices) = featureIntensities(mediumPoints, pointTolerance, featureList)
        (heavyIntensity, heavyFeatures,
         heavyIndices) = featureIntensities(heavyPoints, pointTolerance, featureList)

        
        exampleScannum = spectrumDescriptionToScanNumber((lights+mediums+heavies)[0][0])
        examplePSM = resultIndex[exampleScannum]
        if lightFeatures and (not (mediumFeatures and heavyFeatures)) and collidablePSM(examplePSM):
            notes.append("Potential case of overlapping SILAC features")
            # There used to be code to try to account for this; see old versions in Hg.
            
        if lightFeatures or mediumFeatures or heavyFeatures:
            lightMediumShift, lightHeavyShift = getShifts(exampleScannum)
            (overlapLight, overlapMedium, overlapHeavy, 
             overlapNote) = overlapIntensityTriple(data,
                                                   lightFeatures, mediumFeatures, heavyFeatures,
                                                   lightMediumShift, lightHeavyShift)
            
            if overlapNote:
                notes.append(overlapNote)
                
            if lightFeatures:
                strongestLight = max(lightFeatures, key = lambda x: x.totalIntensity())
            else:
                strongestLight = None
            if mediumFeatures:
                strongestMedium = max(mediumFeatures, key = lambda x: x.totalIntensity())
            else:
                strongestMedium = None
            if heavyFeatures:
                strongestHeavy = max(heavyFeatures, key = lambda x: x.totalIntensity())
            else:
                strongestHeavy = None                
            
            def isStronger(first, second, third):
                return (((not second) or first.totalIntensity() > second.totalIntensity()) and
                        ((not third) or first.totalIntensity() > third.totalIntensity()))
                
                
            if strongestLight and isStronger(strongestLight, strongestMedium, strongestHeavy):
                (topSignalOverlapLight,
                 topSignalOverlapMedium,
                 topSignalOverlapHeavy,
                 topSignalNote) = overlapIntensityTriple(data,
                                                         [strongestLight],
                                                         mediumFeatures,
                                                         heavyFeatures,
                                                         lightMediumShift,
                                                         lightHeavyShift)
            elif strongestMedium and isStronger(strongestMedium, strongestLight, strongestHeavy):
                (topSignalOverlapLight,
                 topSignalOverlapMedium,
                 topSignalOverlapHeavy,
                 topSignalNote) = overlapIntensityTriple(data,
                                                         lightFeatures,
                                                         [strongestMedium],
                                                         heavyFeatures,
                                                         lightMediumShift,
                                                         lightHeavyShift)
            elif strongestHeavy and isStronger(strongestHeavy, strongestLight, strongestMedium):
                (topSignalOverlapLight,
                 topSignalOverlapMedium,
                 topSignalOverlapHeavy,
                 topSignalNote) = overlapIntensityTriple(data,
                                                         lightFeatures,
                                                         mediumFeatures,
                                                         [strongestHeavy],
                                                         lightMediumShift,
                                                         lightHeavyShift)
                 
                
            if topSignalNote:
                notes.append('TS- ' + topSignalNote)
                
                
            xicParameters = getXICParametersTriple(lightFeatures, mediumFeatures, 
                                                   heavyFeatures, lightMediumShift,
                                                   lightHeavyShift)
            (xicByFeature,
             (totalLightXIC, totalMediumXIC, totalHeavyXIC)) = getRatioXICsTriple(data,
                                                                                  xicParameters)
            
            
        else:
            overlapLight, overlapMedium, overlapHeavy = '-', '-', '-'
            xicParameters = []
            (topSignalOverlapLight,
             topSignalOverlapMedium,
             topSignalOverlapHeavy,
             topSignalNote) = ('-', '-', '-', '')
            (totalLightXIC, totalMediumXIC, totalHeavyXIC) = 0, 0, 0
            xicByFeature = []            
            
        
        (lightMediumCorrelation,
         lightHeavyCorrelation,
         lightMediumRatioScore,
         lightHeavyRatioScore) = tripleFeatureSimilarities(lightFeatures,
                                                           mediumFeatures,
                                                           heavyFeatures)
        
        if not lightFeatures: 
            notes.append("No untagged feature")
        if not mediumFeatures:
            notes.append("No medium-tagged feature")
        if not heavyFeatures: 
            notes.append("No heavy-tagged feature")      
            
        if len(lightFeatures) > 1:
            notes.append("Untagged intensity over %s features" % len(lightFeatures))
        if len(mediumFeatures) > 1:
            notes.append("Medium-tag intensity over %s features" % len(mediumFeatures))
        if len(heavyFeatures) > 1:
            notes.append("Heavy-tag intensity over %s features" % len(heavyFeatures))
            
        notes = '; '.join(notes)
        identifier = (resultIndex[spectrumDescriptionToScanNumber(lights[0][0])]['Spectrum Description']
                      if lights else '-',
                      resultIndex[spectrumDescriptionToScanNumber(mediums[0][0])]['Spectrum Description']
                      if mediums else '-',
                      resultIndex[spectrumDescriptionToScanNumber(heavies[0][0])]['Spectrum Description']
                      if heavies else '-')
        

        lightIndices = ';'.join(lightIndices)
        mediumIndices = ';'.join(mediumIndices)
        heavyIndices = ';'.join(heavyIndices)
        
        ratios.append((identifier,
                       lightIntensity, mediumIntensity, heavyIntensity,
                       notes,
                       lightMediumCorrelation, lightHeavyCorrelation, 
                       lightMediumRatioScore, lightHeavyRatioScore, 
                       overlapLight, overlapMedium, overlapHeavy,
                       lightIndices, mediumIndices, heavyIndices,
                       topSignalOverlapLight, topSignalOverlapMedium,
                       topSignalOverlapHeavy,
                       xicByFeature,
                       totalLightXIC, totalMediumXIC, totalHeavyXIC))
        
    return ratios

def getDoubleRatios(data, resultIndex, featureList, ms2toms1, tagTuples):
    featureIndex = featureList # Use of FeatureInterface.
    #featureIndex = ProximityIndexedSequenceAgain(enumerate(featureList),
                                                 #indexer = lambda x: x[1].mz)
    
    def getShift(psm):
        specDesc, _, _ = psm
        chg = float(resultIndex[specDesc]['Charge'])
        seq = resultIndex[specDesc]['Peptide Sequence']
        Rs = seq.count('R')
        Ks = seq.count('K')
        return ((Rs * heavyShifts['R'])/chg + (Ks * heavyShifts['K'])/chg) 
    
    def getAcqPoints(psms, predict = None, recover = False):
        points = []
        for specDesc, _, _ in psms:
            try:
                scan = ms2toms1[spectrumDescriptionToScanNumber(specDesc)]
            except KeyError as err:
                print("Scan number not found in target file.")
                print(specDesc)
                print(spectrumDescriptionToScanNumber(specDesc))
                print(list(ms2toms1.keys())[:5])
                raise err
            
            chg = float(resultIndex[specDesc]['Charge'])
    
            if predict == 'light': sign = -1
            elif predict == 'heavy': sign = 1
            elif predict: raise Exception("Invalid value of predict.")
            
            shift = 0
            if predict:
                seq = resultIndex[specDesc]['Peptide Sequence']
                Rs = seq.count('R')
                Ks = seq.count('K')
                shift = ((Rs * heavyShifts['R'])/chg + (Ks * heavyShifts['K'])/chg) * sign
            
            #mz = float(specDesc.split('|')[1])
            try:
                mz = spectrumDescriptionToMZ(specDesc)
            except IndexError:
                mz = float(resultIndex[specDesc]['Experimental mz'])
            mz += shift
            
            points.append((mz, scan, chg))
            
        return points
    
    def getDeltaTolerance(psmIndexes):
        avgDelta = average([abs(resultIndex[x[0]]['Delta']) for x in psmIndexes])
        avgMZ = average([abs(resultIndex[x[0]]['Experimental mz']) for x in psmIndexes])
        return ((avgDelta*1000000)/avgMZ) + (peakFindTolPPM/2)
        
        
        
    
    
    ratios = []
    for lights, _, heavies in tagTuples:
        if not (lights or heavies): continue
        
        notes = []
        if not heavies:
            notes.append("Extrapolated tagged mz")
        if not lights:
            notes.append("Extrapolated untagged mz")        
        
        lightIntensity = 0
        heavyIntensity = 0
        
        if lights: 
            lightPoints = getAcqPoints(lights) 
        else:
            lightPoints = getAcqPoints(heavies, predict = 'light')
        
        if heavies:
            heavyPoints = getAcqPoints(heavies)
        else:
            heavyPoints = getAcqPoints(lights, predict = 'heavy')
            
        pointTolerance = getDeltaTolerance(lights + heavies)
        
        lightPoints, heavyPoints = interpolatePoints(lightPoints, heavyPoints)
        
        lightIntensity, lightFeatures, lightIndices = featureIntensities(lightPoints, pointTolerance, featureIndex)
        heavyIntensity, heavyFeatures, heavyIndices = featureIntensities(heavyPoints, pointTolerance, featureIndex)
        

        shift = getShift((lights + heavies)[0])
        if lightFeatures or heavyFeatures:
            overlapLight, overlapHeavy, overlapNote = overlapIntensity(data,
                                                                       lightFeatures, 
                                                                       heavyFeatures,
                                                                       shift)
            if overlapNote:
                notes.append(overlapNote)
            
            if lightFeatures:
                strongestLight = max(lightFeatures, key = lambda x: x.totalIntensity())
            else:
                strongestLight = None
            if heavyFeatures:
                strongestHeavy = max(heavyFeatures, key = lambda x: x.totalIntensity())
            else:
                strongestHeavy = None
            if strongestLight and ((not strongestHeavy)
                                   or strongestLight.totalIntensity() > strongestHeavy.totalIntensity()):
                (topSignalOverlapLight,
                 topSignalOverlapHeavy,
                 topSignalNote) = overlapIntensity(data,
                                                   [strongestLight],
                                                   heavyFeatures,
                                                   shift)
            else:
                (topSignalOverlapLight,
                 topSignalOverlapHeavy,
                 topSignalNote) = overlapIntensity(data,
                                                   lightFeatures,
                                                   [strongestHeavy],
                                                   shift)
            
            if topSignalNote:
                notes.append('TS- ' + topSignalNote)
                
            
        
            xicParameters = getXICParametersDouble(lightFeatures, heavyFeatures, shift)            
            xicByFeature, (totalLightXIC, totalHeavyXIC) = getRatioXICsDouble(data, xicParameters)
        else:
            overlapLight, overlapHeavy = '-', '-'
            xicParameters = []
            (topSignalOverlapLight,
             topSignalOverlapHeavy,
             topSignalNote) = ('-', '-', '')
            (totalLightXIC, totalHeavyXIC) = 0, 0
            xicByFeature = []
            

            

        if not lightFeatures: 
            notes.append("No untagged feature")

        if not heavyFeatures: 
            notes.append("No tagged feature")
            
        if len(lightFeatures) > 1:
            notes.append("Untagged intensity over %s features" % len(lightFeatures))
        if len(heavyFeatures) > 1:
            notes.append("Tagged intensity over %s features" % len(heavyFeatures))
            
        if lightFeatures and heavyFeatures:
            correlationScore, ratioScore = featureSimilarities(lightFeatures, heavyFeatures)
        else: correlationScore, ratioScore = '-', '-'
        
        notes = '; '.join(notes)
        
        identifier = (resultIndex[lights[0][0]]['Spectrum Description'] if lights else '-', 
                      resultIndex[heavies[0][0]]['Spectrum Description'] if heavies else '-')        
        
        lightIndices = ';'.join(lightIndices)
        heavyIndices = ';'.join(heavyIndices)        

        ratios.append((identifier, 
                       lightIntensity, heavyIntensity, 
                       notes, 
                       correlationScore,
                       ratioScore, 
                       overlapLight, overlapHeavy,
                       lightIndices, heavyIndices,
                       topSignalOverlapLight, topSignalOverlapHeavy,
                       xicByFeature, totalLightXIC, totalHeavyXIC))
        
    return ratios
            
            
         
         
def getXICParametersDouble(firsts, seconds, shift):
    firstPts = set([(x[1], x[0]) for x in sum([f.allIndexedPoints() for f in firsts], [])])
    secondPts = set([(x[1], x[0]) for x in sum([f.allIndexedPoints() for f in seconds], [])])
    
    try:
        firstMin = min([x[0] for x in firstPts])
        firstMZ = average([x[0] for x in firstPts if abs(x[0] - firstMin) < 0.1])
    except ValueError:
        firstMZ = None
    
    try:
        secondMin = min([x[0] for x in secondPts]) 
        secondMZ = average([x[0] for x in secondPts if abs(x[0] - secondMin) < 0.1]) 
    except ValueError:
        secondMZ = None
    
    if firstMZ and secondMZ:
        assert abs((firstMZ + shift) - secondMZ) < 0.1, "Non-C12 feature MZ? %s" % abs((firstMZ + shift) - secondMZ)
    else:
        if not firstMZ:
            firstMZ = secondMZ - shift
        else:
            secondMZ = firstMZ + shift
            
    
    # If this were to strictly adhere to the old system, it would only
    # take scans where the C12 pt was present.
    scans = set([x[1] for x in firstPts] + [x[1] for x in secondPts])
    
    startscan, endscan = min(scans), max(scans)
    
    return [(0, startscan, endscan, firstMZ, secondMZ)]
    
            
    
def getXICParametersTriple(lights, mediums, heavies, lightMediumShift, lightHeavyShift):
    def featurePoints(features):
        return set([(x[1], x[0]) for x in sum([f.allIndexedPoints() for f in features], [])])
    lPts = featurePoints(lights)
    mPts = featurePoints(mediums)
    hPts = featurePoints(heavies)
    
    try:
        lMin = min(select(0, lPts))
        lMZ = average([x[0] for x in lPts if abs(x[0] - lMin) < 0.1])
    except ValueError:
        lMZ = None
    try:
        mMin = min(select(0, mPts))
        mMZ = average([x[0] for x in mPts if abs(x[0] - mMin) < 0.1])
    except ValueError:
        mMZ = None
    try:
        hMin = min(select(0, hPts))
        hMZ = average([x[0] for x in hPts if abs(x[0] - hMin) < 0.1])
    except ValueError:
        hMZ = None
    
    
    if not lMZ:
        if mMZ:
            lMZ = mMZ - lightMediumShift
        elif hMin:
            lMZ = hMZ - lightHeavyShift
    if not mMZ:
        if lMZ:
            mMZ = lMZ + lightMediumShift
        elif hMZ:
            mMZ = hMZ - (lightHeavyShift - lightMediumShift)
    if not hMZ:
        if lMZ:
            hMZ = lMZ + lightHeavyShift
        elif mMZ:
            hMZ = mMZ + (lightHeavyShift - lightMediumShift)
    
    assert all([lMZ, mMZ, hMZ])
    
    scans = set([x[1] for x in lPts] + [x[1] for x in mPts] + [x[1] for x in hPts])
    startscan, endscan = min(scans), max(scans)
    
    return [(0, startscan, endscan, lMZ, mMZ, hMZ)]
    
        
    
    
        
def getMS1Lookup(datafile):
    data = mzFile(datafile)
    scans = data.scan_info(0, 9999999)
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
            raise Exception("Unidentified scan type of %s" % scan[3])
    for ms2 in ms2s:
        ms2toms1[ms2] = ms1
        
    # This should be gotten rid of when the MS1-relative scan numbering
    # in feature detector output is fixed.
    absMS1Lookup = dict(enumerate([x[2] for x in scans if x[3] == 'MS1']))

        
    return ms2toms1, absMS1Lookup


def showRatioHistogram(ratioList):
    ratios = [sorted([x[1], x[2]]) for x in ratioList if x[1] or x[2]]
    ratios = [x[0]/x[1] for x in ratios if x[1] and x[0] != '-' and x[1] != '-']
    
    print(len([x for x in ratios if x > 0.6]))    
    
    pyt.hist(ratios, bins=100)
    pyt.show()
    
def showOverlapHistogram(ratioList):
    ratios = [sorted([x[6], x[7]]) for x in ratioList if x[6] or x[7]]
    ratios = [x[0]/x[1] for x in ratios if x[1] and x[0] != '-' and x[1] != '-']
    
    print(len([x for x in ratios if x > 0.6]))
    
    pyt.hist(ratios, bins=100)
    pyt.show()
    
    
    
def renderC12XIC(xics):
    if xics:
        assert xics[0][3] == min([x[3] for x in xics])
    
        ratio = xics[0][-2] / xics[0][-1] if xics[0][-1] else float('nan')
        return str(round(ratio, 2))
    else:
        return '-'
    
def renderC12XICTriple(xics, mode):
    if not xics:
        return '-'
    else:
        xic = xics[0]
        assert xic[3] == min([x[3] for x in xics]) # Dunno why?
        
        if mode == 'ml':
            ratio = xic[-2]/xic[-3] if xic[-3] else float('nan')
        elif mode == 'hl':
            ratio = xic[-1]/xic[-3] if xic[-3] else float('nan')
        return str(round(ratio, 2))
    

def renderXICs(xics):
    output = []
    isotopes = ['C12', 'C13'] + ['%d-C13' % x for x in range(2, 20)]
    for index, startScan, endScan, lightMZ, heavyMZ, lightAOC, heavyAOC in xics:
        ratio = lightAOC / heavyAOC if heavyAOC else float('nan') # Reverse this?
        token = '%s : %.2f' % (isotopes[index], ratio)
        output.append(token)
    return '|'.join(output)
    

    
def renderTripleXICs(xics, mode):
    output = []
    isotopes = ['C12', 'C13'] + ['%d-C13' % x for x in range(2, 20)]
    
    for index, startScan, endScan, lightMZ, mediumMZ, heavyMZ, lightAOC, mediumAOC, heavyAOC in xics:
        if mode == 'hl':
            ratio = heavyAOC / lightAOC if lightAOC else float('nan')
        elif mode == 'ml':
            ratio = mediumAOC / lightAOC if lightAOC else float('nan')
        else:
            raise Exception
        token = "%s : %.2f" % (isotopes[index], ratio)
        output.append(token)
    
    return '|'.join(output)
    
def writeDuplexSILAC(resultIndex, ratios, columns, outputFile):    
    output = writer(outputFile,
                    columns = columns +
                    ['Light Intensity', 'Heavy Intensity', 'SILAC notes', 
                     'Elution Correlation', 'Ratio Variance',
                     'Light Features', 'Heavy Features',
                     'Light Overlap Intensity', 'Heavy Overlap Intensity',
                     'Light Top-Signal Intensity', 'Heavy Top-Signal Intensity',
                     'Light XIC Intensity', 'Heavy XIC Intensity',
                     'Direct Ratio (H/L)', 'Overlap Ratio (H/L)', 'XIC Ratio (H/L)'])

    molecIndex = {}
    #for identifier, light, heavy, notes, correlation, ratioMatch, over1, over2, xics in ratios:
    for ratio in ratios:
        identifier = ratio[0]
        if identifier[0] != '-':
            lightPSM = resultIndex[identifier[0]]
            lightMol = lightPSM['Peptide Sequence'], lightPSM['Variable Modifications'], lightPSM['Charge']
            molecIndex[lightMol] = ratio[1:]
        
        if identifier[1] != '-':
            heavyPSM = resultIndex[identifier[1]]
            heavyMol = heavyPSM['Peptide Sequence'], heavyPSM['Variable Modifications'], heavyPSM['Charge']
            molecIndex[heavyMol] = ratio[1:]
    
    for line in list(resultIndex.values()):
        molecule = line['Peptide Sequence'], line['Variable Modifications'], line['Charge']
        try:
            (light, heavy, notes, correlation, ratioMatch,
             overlap1, overlap2, lightFeatures, heavyFeatures,
             lightTopSignal, heavyTopSignal,
             xics, totalLightAOC, totalHeavyAOC) = molecIndex[molecule]
        except KeyError:
            (light, heavy, notes, correlation, ratioMatch,
             overlap1, overlap2, lightFeatures, heavyFeatures,
             lightTopSignal, heavyTopSignal,
             xics, totalLightAOC, totalHeavyAOC) = ("-", "-", "No SILAC pair.", "-", "-",
                                                    "-", "-", "-", "-", '-', '-',
                                                    [], '', '')
        
        line['Light Intensity'] = light
        line['Heavy Intensity'] = heavy
        line['SILAC notes'] = notes
        line['Elution Correlation'] = correlation
        line['Ratio Variance'] = ratioMatch
        line['Light Overlap Intensity'] = overlap1
        line['Heavy Overlap Intensity'] = overlap2
        line['Light Top-Signal Intensity'] = lightTopSignal
        line['Heavy Top-Signal Intensity'] = heavyTopSignal        
        #line['C12 XIC Ratio'] = renderC12XIC(xics)
        #line['XIC Ratios'] = renderXICs(xics)
        line['Light Features'] = lightFeatures
        line['Heavy Features'] = heavyFeatures
        
        if isinstance(overlap1, float):
            line['Direct Ratio (H/L)'] = heavy / light if light else 'NaN'
            line['Overlap Ratio (H/L)'] = overlap2 / overlap1 if overlap1 else 'NaN'
            line['Light XIC Intensity'] = totalLightAOC
            line['Heavy XIC Intensity'] = totalHeavyAOC            
            line['XIC Ratio (H/L)'] = totalHeavyAOC / totalLightAOC if totalLightAOC else 'NaN'
        else:
            line['Light XIC Intensity'] = '-'
            line['Heavy XIC Intensity'] = '-'
            line['Direct Ratio (H/L)'] = '-'
            line['Overlap Ratio (H/L)'] = '-'            
            line['XIC Ratio (H/L)'] = '-'

        output.write(line)
        
    output.close()
    
def writeTriplexSILAC(resultIndex, ratios, columns, outputFile):
    output = writer(outputFile,
                    columns = columns +
                    ['Light Intensity', 'Medium Intensity', 'Heavy Intensity', 'SILAC notes', 
                     'Light-Medium Elution Correlation', 'Light-Heavy Elution Correlation', 
                     'Light-Medium Ratio Variance', 'Light-Heavy Ratio Variance',
                     'Light Features', 'Medium Features', 'Heavy Features',
                     'Light Overlap Intensity', 'Medium Overlap Intensity',
                     'Heavy Overlap Intensity',
                     'C12 XIC Ratio (H/L)', 'XIC Ratios (H/L)',
                     'C12 XIC Ratio (M/L)', 'XIC Ratios (M/L)',
                     'Direct Ratio (H/L)', 'Overlap Ratio (H/L)', 'XIC Ratio (H/L)',
                     'Direct Ratio (M/L)', 'Overlap Ratio (M/L)', 'XIC Ratio (M/L)'
                     ])
    
    molecIndex = {}
    for ratioData in ratios:
        identifier = ratioData[0]
        
        if identifier[0] != '-':
            lightPSM = resultIndex[identifier[0]]
            lightMol = lightPSM['Peptide Sequence'], lightPSM['Variable Modifications'], lightPSM['Charge']
            molecIndex[lightMol] = ratioData
        
        if identifier[1] != '-':
            mediumPSM = resultIndex[identifier[1]]
            mediumMol = mediumPSM['Peptide Sequence'], mediumPSM['Variable Modifications'], mediumPSM['Charge']
            molecIndex[mediumMol] = ratioData
            
        if identifier[2] != '-':
            heavyPSM = resultIndex[identifier[2]]
            heavyMol = heavyPSM['Peptide Sequence'], heavyPSM['Variable Modifications'], heavyPSM['Charge']
            molecIndex[heavyMol] = ratioData
            
    for line in list(resultIndex.values()):
        molecule = line['Peptide Sequence'], line['Variable Modifications'], line['Charge']
        try:
            (identifier, light, medium, heavy, notes, 
             correlationLM, correlationLH, ratioMatchLM, ratioMatchLH,
             overlap1, overlap2, overlap3,
             lightFe, mediumFe, heavyFe,
             xics, totalLightAOC, totalMediumAOC,
             totalHeavyAOC) = molecIndex[molecule]
        except KeyError:
            (identifier, light, medium, heavy, notes, 
             correlationLM, correlationLH, ratioMatchLM, ratioMatchLH,
             overlap1, overlap2, overlap3,
             lightFe, mediumFe, heavyFe,
             xics, totalLightAOC, totalMediumAOC,
             totalHeavyAOC) = ("-", "-", "-", "-",
                               "No SILAC pair.", 
                               "-", "-", "-", "-", "-", "",
                               "", "", "", "", "", "")
            
        line['Light Intensity'] = light
        line['Medium Intensity'] =  medium
        line['Heavy Intensity'] = heavy
        line['SILAC notes'] = notes
        line['Light-Medium Elution Correlation'] = correlationLM
        line['Light-Heavy Elution Correlation'] = correlationLH
        line['Light-Medium Ratio Variance'] = ratioMatchLM
        line['Light-Heavy Ratio Variance'] = ratioMatchLH
        line['Light Overlap Intensity'] = overlap1
        line['Medium Overlap Intensity'] = overlap2
        line['Heavy Overlap Intensity'] = overlap3
        line['C12 XIC Ratio (M/L)'] = renderC12XICTriple(xics, 'ml')
        line['C12 XIC Ratio (H/L)'] = renderC12XICTriple(xics, 'hl')
        line['XIC Ratios (H/L)'] = renderTripleXICs(xics, 'hl')
        line['XIC Ratios (M/L)'] = renderTripleXICs(xics, 'ml')
        line['Light Features'] = lightFe
        line['Medium Features'] = mediumFe
        line['Heavy Features'] = heavyFe
        
        if isinstance(overlap1, float):
            line['Direct Ratio (H/L)'] = heavy / light if light else 'NaN'
            line['Overlap Ratio (H/L)'] = overlap3 / overlap1 if overlap1 else 'NaN'
            line['XIC Ratio (H/L)'] = totalHeavyAOC / totalLightAOC if totalLightAOC else 'NaN'
            line['Direct Ratio (M/L)'] = medium / light if light else 'NaN'
            line['Overlap Ratio (M/L)'] = overlap2 / overlap1 if overlap1 else 'NaN'
            line['XIC Ratio (M/L)'] = totalMediumAOC / totalLightAOC if totalLightAOC else 'NaN'
        else:
            line['Direct Ratio (H/L)'] = '-'
            line['Overlap Ratio (H/L)'] = '-'
            line['XIC Ratio (H/L)'] = '-'
            line['Direct Ratio (M/L)'] = '-'
            line['Overlap Ratio (M/L)'] = '-'
            line['XIC Ratio (M/L)'] = '-'
        
        
        output.write(line)
    
    output.close()
        
        
    
def combineTags(Rs, Ks):
    sR = set([tuple(map(tuple, x)) for x in Rs])
    sK = set([tuple(map(tuple, x)) for x in Ks])
    
    sBoth = sR | sK
    
    return [tuple(map(list, x)) for x in sBoth]




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
            return int(scanRegCompiled.search(description).group())
        spectrumDescriptionToScanNumber = newParser
        
    if 'tolerance' in constants:
        #peakFindTolerance = constants['tolerance']
        #XICTol = peakFindTolerance/2
        global peakFindTolPPM
        global XICTol
        peakFindTolPPM = constants['tolerance']
        if peakFindTolPPM < 1:
            print("\n\n\nWARNING- tolerance value for SILAC analysis should be in PPM!\n\n\n")
        XICTol = 0.0008 * peakFindTolPPM # "Typical" Da tolernace for PPM.
        
        

def getRatioXICsDouble(data, xicparameters):
    def AOC(xic):
        aoc = 0
        for i in range(0, len(xic)-1):
            cur = xic[i]
            next = xic[i+1]
            width = next[0] - cur[0]
            height = cur[1]
            aoc += width * height
        return aoc    
    
    totalLightAOC = 0
    totalHeavyAOC = 0
    #for j in range(0, len(xicparameters)):
    xicForFeature = []
    for index, startScan, endScan, lightMZ, heavyMZ in xicparameters:
        try:
            startRT = data.timeForScan(startScan)
        except IOError:
            warnings.warn('Feature scan out of range: %s %s' % (datafile, startScan))
            startRT = data.timeForScan(startScan+1)
        try:
            stopRT = data.timeForScan(endScan)
        except IOError:
            warnings.warn('Feature scan out of range: %s %s' % (datafile, endScan))
            stopRT = data.timeForScan(endScan-1)
            
        try:
            lightXIC = data.xic(startRT, stopRT, lightMZ - XICTol, lightMZ + XICTol)
            heavyXIC = data.xic(startRT, stopRT, heavyMZ - XICTol, heavyMZ + XICTol)
        except AssertionError:
            lightXIC = []
            heavyXIC = []
        
        lightAOC = AOC(lightXIC)
        heavyAOC = AOC(heavyXIC)
        
        totalLightAOC += lightAOC
        totalHeavyAOC += heavyAOC
        
        xicForFeature.append((index, startScan, endScan, lightMZ, heavyMZ, lightAOC, heavyAOC))
    
    return xicForFeature, (totalLightAOC, totalHeavyAOC)
    
    
                
                
def getRatioXICsTriple(data, ratios):
    def AOC(xic):
        aoc = 0
        for i in range(0, len(xic)-1):
            cur = xic[i]
            next = xic[i+1]
            width = next[0] - cur[0]
            height = cur[1]
            aoc += width * height
        return aoc
    
    lightAOC = 0
    mediumAOC = 0
    heavyAOC = 0
    
    xicForFeature = []    

    totalLightAOC = 0
    totalMediumAOC = 0
    totalHeavyAOC = 0
    for index, startScan, endScan, lightMZ, mediumMZ, heavyMZ in ratios:
        startRT = data.timeForScan(startScan)
        stopRT = data.timeForScan(endScan)
        lightXIC = data.xic(startRT, stopRT, lightMZ - XICTol, lightMZ + XICTol)
        mediumXIC = data.xic(startRT, stopRT, mediumMZ - XICTol, mediumMZ + XICTol)
        heavyXIC = data.xic(startRT, stopRT, heavyMZ - XICTol, heavyMZ + XICTol)
        
        lightAOC = AOC(lightXIC)
        mediumAOC = AOC(mediumXIC)
        heavyAOC = AOC(heavyXIC)
        totalLightAOC += lightAOC
        totalMediumAOC += mediumAOC
        totalHeavyAOC += heavyAOC
        
        xicForFeature.append((index, startScan, endScan, lightMZ, mediumMZ, heavyMZ, 
                              lightAOC, mediumAOC, heavyAOC))
    
    return xicForFeature, (totalLightAOC, totalMediumAOC, totalHeavyAOC)



def label_shift(label):
    isotopes = re.findall('([A-Z]\([0-9]+\))', label)
    shift = 0
    for isostr in isotopes:
        # Should be in the form "C(8)" or similar.
        el = isostr[0]
        count = int(isostr[2:-1])
        shift += isotopeDiffs[el] * count
    return shift
        
    
    
def silac_associated_mzs(psms, heavyK, heavyR, mediumK = None, mediumR = None, round_and_condense = False):
    if mediumK and mediumR:
        plex = 3
        medium_shifts = {'K':label_shift(mediumK), 'R':label_shift(mediumR)}
    else:
        assert not (mediumK or mediumR), "Did you mean to do 3-plex?"
        plex = 2
        medium_shifts = None
    heavy_shifts = {'K':label_shift(heavyK), 'R':label_shift(heavyR)}
    
    mzs = []
    for psm in psms:
        mz = float(psm['Experimental mz'])
        mzs.append(mz)
        
        varmods = psm['Variable Modifications']
        pepseq = psm['Peptide Sequence']
        
        if 'K' in pepseq or 'R' in pepseq:
            kcount = pepseq.count('K')
            rcount = pepseq.count('R')
            heavy = heavy_shifts['K']*kcount + heavy_shifts['R']*rcount
            if medium_shifts:
                medium = medium_shifts['R']*rcount + medium_shifts['K']*kcount
            if any([x in varmods for x in [heavyK, heavyR]]):
                mzs.append(mz - heavy)
                if medium_shifts:
                    mzs.append((mz - heavy) + medium)
            elif medium_shifts and any([x in varmods for x in [mediumK, mediumR]]):
                mzs.append(mz - medium)
                mzs.append((mz - medium) + heavy)
            else:
                if medium_shifts:
                    mzs.append(mz + medium)
                mzs.append(mz + heavy)
    
    if round_and_condense:
        mzs = list(set([round(x, 2) for x in mzs]))
    return mzs
                
        
        #for aa, label, shift in medium_shifts + heavy_shifts:
            #sign = -1 if label in varmods else 1
            
            
            
        




    
    
def SILAC2Plex(datafiles, resultfiles, heavyTags, **constants):
    global heavyR, heavyK, allTags
    global heavyShifts
    
    if isinstance(datafiles, str):
        datafiles = [datafiles]
    if isinstance(resultfiles, str):
        resultfiles = [resultfiles]
    
    if isinstance(heavyTags, list):
        heavyK, heavyR = heavyTags
    elif isinstance(heavyTags, dict):
        assert set(heavyTags.keys()) == set(['K', 'R']), 'Only K and R labelling supported.'
        heavyK = heavyTags['K']
        heavyR = heavyTags['R']
    allTags = [heavyR, heavyK]
    heavyShifts = {'K' : unimod.get_mod_delta(heavyK),
                   'R' : unimod.get_mod_delta(heavyR)}     
    
    if constants:
        setGlobals(constants)
    if 'whitelist_psms' in constants and constants['whitelist_psms']:
        psms = sum([list(reader(x)) for x in resultfiles], [])
        mzs = silac_associated_mzs(psms, heavyK = heavyK, heavyR = heavyR)
        constants['whitelist_psms'] = mzs
     
    
    ms2toms1_byfile = {}
    absMS1Lookup_byfile = {}
    featureMemo_byfile = {}
    for datafile in datafiles:
        datafilebase = os.path.basename(datafile)
        (ms2toms1_byfile[datafilebase],
         absMS1Lookup_byfile[datafilebase]) = getMS1Lookup(datafile)
    
        featureMemo_byfile[datafilebase] = FeatureMemo(detectFeatures(datafile,
                                                                      **constants))

    # This is fairly odd, but just for the sake of least possible
    # assumptions: all result files are combined together, and then parcelled
    # out by their source raw file for analysis. Analyzed results are then
    # brought together into a single annotated result file.
    combinedResults = []
    resultToData = {}
    for resultFile in resultfiles:
        rdr = reader(resultFile)
        #assert 'Source' in rdr.columns
        psms = list(rdr)
        
        if 'Source' in psms[0]:
            partialResultFiles = set([x['Source'] for x in psms])
        else:
            partialResultFiles = set([resultFile])
            rdr.columns.insert(0, 'Source')
            for psm in psms:
                psm['Source'] = resultFile
            
        for partialResultFile in partialResultFiles:
            #datafileList = [x for x in datafiles if os.path.basename(x).split('.')[0] in partialResultFile]
            #assert len(datafileList) == 1
            if partialResultFile == 'None':
                assert len(datafiles) == 1
                datafile = datafiles[0]
            else:
                # Will fail if datafile base name and PSM file base name don't match.
                # Used to use 'in' instead of equality, but this lead to problems
                # for, e.g., 'some_fractionated_run_1' vs 'some_fractionated_run_12'.
                datafile = [x for x in datafiles if
                            os.path.basename(partialResultFile).split('.')[0] == os.path.basename(x).split('.')[0]][0] 
                
            resultToData[partialResultFile] = datafile
        combinedResults += psms
    
    assert len(set(resultToData.values())) == len(datafiles), str(list(resultToData.values()))
    
    
    results_byfile = defaultdict(list)
    for psm in combinedResults:
        results_byfile[os.path.basename(resultToData[psm['Source']])].append(psm)    
    results_byfile = dict(results_byfile)
        
    
    ratios_byfile = {}
    resultIndex_byfile = {}
    for datafilebase, psms in list(results_byfile.items()):
        assert datafilebase == os.path.basename(resultToData[psms[0]['Source']])
        datafile = resultToData[psms[0]['Source']]
        
        #resultIndex = {x['Spectrum Description']:x for x in psms}   
        resultIndex = {}
        data = mzFile(datafile)
        for psm in psms:
            desc = psm['Spectrum Description']
            psm['Retention Time'] = data.timeForScan(spectrumDescriptionToScanNumber(desc))
            resultIndex[desc] = psm
        
        taggedPSMs = findDoublesAdapter(resultIndex)
        vprint("Got %s doubles." % len(taggedPSMs))
        
        ratios = getDoubleRatios(data, resultIndex, featureMemo_byfile[datafilebase],
                                 ms2toms1_byfile[datafilebase], taggedPSMs)
        #ratios = getRatioXICsDouble(datafile, ratios)
        
        ratios_byfile[datafilebase] = ratios
        resultIndex_byfile[datafilebase] = resultIndex
    
    if len(resultfiles) == 1:
        outputfile = '.'.join(resultfiles[0].split('.')  + ['SILAC_annotated', 'xlsx'])
    else:
        outputfile = os.path.join(os.path.dirname(resultfiles[0]),
                                  'combined_results.SILAC_annotated.xlsx')    
    
    return writeCombinedDuplexSILAC(resultIndex_byfile, results_byfile, ratios_byfile, rdr.columns, outputfile)
    
    
    
def SILAC3Plex(datafiles, resultfiles, mediumTags, heavyTags, **constants):
    global heavyR, heavyK, mediumR, mediumK, allTags
    global mediumShifts, heavyShifts
    
    if isinstance(datafiles, str):
        datafiles = [datafiles]
    if isinstance(resultfiles, str):
        resultfiles = [resultfiles]    
    
    if isinstance(heavyTags, list):
        heavyK, heavyR = heavyTags
        mediumK, mediumR = mediumTags
    elif isinstance(heavyTags, dict):
        assert set(heavyTags.keys()) == set(['K', 'R']), 'Only K and R labelling supported.'
        assert set(mediumTags.keys()) == set(['K', 'R']), 'Only K and R labelling supported.'
        heavyK = heavyTags['K']
        heavyR = heavyTags['R']
        mediumK = mediumTags['K']
        mediumR = mediumTags['R']
    allTags = [heavyR, heavyK, mediumR, mediumK]
    
    mediumShifts = {'K' : unimod.get_mod_delta(mediumK),
                    'R' : unimod.get_mod_delta(mediumR)}
    heavyShifts = {'K' : unimod.get_mod_delta(heavyK),
                   'R' : unimod.get_mod_delta(heavyR)}  
    
    if constants:
        setGlobals(constants)
        

    ms2toms1_byfile = {}
    absMS1Lookup_byfile = {}
    featureMemo_byfile = {}
    for datafile in datafiles:
        datafilebase = os.path.basename(datafile)
        (ms2toms1_byfile[datafilebase],
         absMS1Lookup_byfile[datafilebase]) = getMS1Lookup(datafile)
    
        featureMemo_byfile[datafilebase] = FeatureMemo(detectFeatures(datafile,
                                                                      **constants))

    # This is fairly odd, but just for the sake of least possible
    # assumptions: all result files are combined together, and then parcelled
    # out by their source raw file for analysis. Analyzed results are then
    # brought together into a single annotated result file.
    combinedResults = []
    resultToData = {}
    for resultFile in resultfiles:
        rdr = reader(resultFile)
        #assert 'Source' in rdr.columns
        psms = list(rdr)
        partialResultFiles = set([x.get('Source', 'None') for x in psms])
        for partialResultFile in partialResultFiles:
            if partialResultFile == 'None':
                assert len(datafiles) == 1
                datafile = datafiles[0]
            else:
                datafile = [x for x in datafiles if
                            os.path.basename(partialResultFile).split('.')[0] in x][0]
            resultToData[partialResultFile] = datafile
        combinedResults += psms
    
    results_byfile = defaultdict(list)
    for psm in combinedResults:
        results_byfile[os.path.basename(resultToData[psm.get('Source', 'None')])].append(psm)    
    results_byfile = dict(results_byfile)
    
    
    ratios_byfile = {}
    resultIndex_byfile = {}
    for datafilebase, psms in list(results_byfile.items()):
        assert datafilebase == os.path.basename(resultToData[psms[0].get('Source', 'None')])
        
        #resultIndex = {x['Spectrum Description']:x for x in psms}
        resultIndex = {}
        data = mzFile(datafile)
        for psm in psms:
            desc = spectrumDescriptionToScanNumber(psm['Spectrum Description'])
            psm['Retention Time'] = data.timeForScan(desc)
            resultIndex[desc] = psm
        
        taggedPSMs = findTriplesAdapter(resultIndex)
        vprint("Got %s triples." % len(taggedPSMs))
        
        ratios = getTripleRatios(data, resultIndex, featureMemo_byfile[datafilebase],
                                 ms2toms1_byfile[datafilebase], taggedPSMs)
        #ratios = getRatioXICsTriple(datafile, ratios)
        
        ratios_byfile[datafilebase] = ratios
        resultIndex_byfile[datafilebase] = resultIndex
    
    if len(resultfiles) == 1:
        outputfile = '.'.join(resultfiles[0].split('.')  + ['SILAC_annotated', 'xlsx'])
    else:
        outputfile = os.path.join(os.path.dirname(resultfiles[0]),
                                  'combined_results.SILAC_annotated.xlsx')
        
    return writeCombinedTriplexSILAC(resultIndex_byfile, results_byfile, ratios_byfile, rdr.columns, outputfile)
        
    
        
        
    
def writeCombinedTriplexSILAC(resultIndex_byfile, results_byfile, ratios_byfile, columns, outputfile):
    output = writer(outputfile,
                    columns = columns +    
                    ['Light Intensity', 'Medium Intensity', 'Heavy Intensity', 'SILAC notes', 
                     'Light-Medium Elution Correlation', 'Light-Heavy Elution Correlation', 
                     'Light-Medium Ratio Variance', 'Light-Heavy Ratio Variance',
                     'Light Features', 'Medium Features', 'Heavy Features',
                     'Retention Time',
                     'Light Overlap Intensity', 'Medium Overlap Intensity',
                     'Heavy Overlap Intensity',
                     'Light Top-Signal Intensity', 'Medium Top-Signal Intensity',
                     'Heavy Top-Signal Intensity',
                     'Direct Ratio (H/L)', 'Overlap Ratio (H/L)',
                     'XIC Ratio (H/L)', 'Top-Signal Ratio (H/L)',
                     'Direct Ratio (M/L)', 'Overlap Ratio (M/L)',
                     'XIC Ratio (M/L)', 'Top-Signal Ratio (M/L)',
                     ])
    
    for datafile, psms in list(resultIndex_byfile.items()):
        molecIndex = {}
        for ratioData in ratios_byfile[datafile]:
            identifier = ratioData[0]
            
            if identifier[0] != '-':
                lightPSM = psms[spectrumDescriptionToScanNumber(identifier[0])]
                lightMol = lightPSM['Peptide Sequence'], lightPSM['Variable Modifications'], lightPSM['Charge']
                molecIndex[lightMol] = ratioData
            
            if identifier[1] != '-':
                mediumPSM = psms[spectrumDescriptionToScanNumber(identifier[1])]
                mediumMol = mediumPSM['Peptide Sequence'], mediumPSM['Variable Modifications'], mediumPSM['Charge']
                molecIndex[mediumMol] = ratioData
                
            if identifier[2] != '-':
                heavyPSM = psms[spectrumDescriptionToScanNumber(identifier[2])]
                heavyMol = heavyPSM['Peptide Sequence'], heavyPSM['Variable Modifications'], heavyPSM['Charge']
                molecIndex[heavyMol] = ratioData            
        
        psmList = results_byfile[datafile]
        for line in psmList:
            molecule = line['Peptide Sequence'], line['Variable Modifications'], line['Charge']
            try:
                (identifier, light, medium, heavy, notes, 
                 correlationLM, correlationLH, ratioMatchLM, ratioMatchLH,
                 overlap1, overlap2, overlap3,
                 lightFe, mediumFe, heavyFe,
                 tsOverLight, tsOverMedium, tsOverHeavy,
                 xics, totalLightAOC, totalMediumAOC,
                 totalHeavyAOC) = molecIndex[molecule]
            except KeyError:
                (identifier, light, medium, heavy, notes, 
                 correlationLM, correlationLH, ratioMatchLM, ratioMatchLH,
                 overlap1, overlap2, overlap3,
                 lightFe, mediumFe, heavyFe,
                 tsOverLight, tsOverMedium, tsOverHeavy,
                 xics, totalLightAOC, totalMediumAOC,
                 totalHeavyAOC) = ("-", "-", "-", "-",
                                   "No SILAC pair.", 
                                   "-", "-", "-", "-", "-", "",
                                   '-','-','-',
                                   "", "", "", "", "", "", '', '')
                
            line['Light Intensity'] = light
            line['Medium Intensity'] =  medium
            line['Heavy Intensity'] = heavy
            line['SILAC notes'] = notes
            line['Light-Medium Elution Correlation'] = correlationLM
            line['Light-Heavy Elution Correlation'] = correlationLH
            line['Light-Medium Ratio Variance'] = ratioMatchLM
            line['Light-Heavy Ratio Variance'] = ratioMatchLH
            line['Light Overlap Intensity'] = overlap1
            line['Medium Overlap Intensity'] = overlap2
            line['Heavy Overlap Intensity'] = overlap3
            line['XIC Ratio (M/L)'] = renderC12XICTriple(xics, 'ml')
            line['XIC Ratio (H/L)'] = renderC12XICTriple(xics, 'hl')
            line['Light Features'] = lightFe
            line['Medium Features'] = mediumFe
            line['Heavy Features'] = heavyFe
            line['Light Top-Signal Intensity'] = tsOverLight
            line['Medium Top-Signal Intensity'] = tsOverMedium
            line['Heavy Top-Signal Intensity'] = tsOverHeavy
            
            if isinstance(overlap1, float):
                line['Direct Ratio (H/L)'] = heavy / light if light else 'NaN'
                line['Overlap Ratio (H/L)'] = overlap3 / overlap1 if overlap1 else 'NaN'
                line['XIC Ratio (H/L)'] = totalHeavyAOC / totalLightAOC if totalLightAOC else 'NaN'
                line['Direct Ratio (M/L)'] = medium / light if light else 'NaN'
                line['Overlap Ratio (M/L)'] = overlap2 / overlap1 if overlap1 else 'NaN'
                line['XIC Ratio (M/L)'] = totalMediumAOC / totalLightAOC if totalLightAOC else 'NaN'
                line['Top-Signal Ratio (H/L)'] = tsOverHeavy / tsOverLight if tsOverLight and tsOverLight != '-' else 'NaN'
                line['Top-Signal Ratio (M/L)'] = tsOverMedium / tsOverLight if tsOverLight and tsOverLight != '-' else 'NaN'
            else:
                line['Direct Ratio (H/L)'] = '-'
                line['Overlap Ratio (H/L)'] = '-'
                line['XIC Ratio (H/L)'] = '-'
                line['Direct Ratio (M/L)'] = '-'
                line['Overlap Ratio (M/L)'] = '-'
                line['XIC Ratio (M/L)'] = '-'
                line['Top-Signal Ratio (H/L)'] = '-'
                line['Top-Signal Ratio (M/L)'] = '-'
            
            output.write(line)            
    
    
    output.close()
    
    return outputfile

def writeCombinedDuplexSILAC(resultIndex_byfile, results_byfile, ratios_byfile, columns, outputFile):
    addedColumns = ['Retention Time',
                    'Light Intensity', 'Heavy Intensity', 'SILAC notes', 
                    'Elution Correlation', 'Ratio Variance',
                    'Light Features', 'Heavy Features',
                    'Light Overlap Intensity', 'Heavy Overlap Intensity',
                    'Light Top-Signal Intensity', 'Heavy Top-Signal Intensity',
                    'Light XIC Intensity', 'Heavy XIC Intensity',
                    'Direct Ratio (H/L)', 'Overlap Ratio (H/L)', 'XIC Ratio (H/L)',
                    'Top-Signal Ratio (H/L)']
    oldColumns = [x for x in columns if x not in addedColumns]
    
    output = writer(outputFile,
                    columns = oldColumns + addedColumns)  
    
    
    for datafile, psms in list(resultIndex_byfile.items()):
        molecIndex = {}
        for ratioData in ratios_byfile[datafile]:
            identifier = ratioData[0]
            
            if identifier[0] != '-':
                lightPSM = psms[identifier[0]]
                lightMol = lightPSM['Peptide Sequence'], lightPSM['Variable Modifications'], lightPSM['Charge']
                molecIndex[lightMol] = ratioData[1:]
                
            if identifier[1] != '-':
                heavyPSM = psms[identifier[1]]
                heavyMol = heavyPSM['Peptide Sequence'], heavyPSM['Variable Modifications'], heavyPSM['Charge']
                molecIndex[heavyMol] = ratioData[1:]      
        
        psmList = results_byfile[datafile] # Full list of PSMs rather than exemplar-per-peptide.
        for line in psmList:
            molecule = line['Peptide Sequence'], line['Variable Modifications'], line['Charge']
            try:
                (light, heavy, notes, correlation, ratioMatch,
                 overlap1, overlap2, lightFeatures, heavyFeatures,
                 lightTopSignal, heavyTopSignal,
                 xics, totalLightAOC, totalHeavyAOC) = molecIndex[molecule]
            except KeyError:
                (light, heavy, notes, correlation, ratioMatch,
                 overlap1, overlap2, lightFeatures, heavyFeatures,
                 lightTopSignal, heavyTopSignal,
                 xics, totalLightAOC, totalHeavyAOC) = ("-", "-", "No SILAC pair.", "-", "-",
                                                        "-", "-", "-", "-", '-', '-',
                                                        [], '', '')
            
            line['Light Intensity'] = light
            line['Heavy Intensity'] = heavy
            line['SILAC notes'] = notes
            line['Elution Correlation'] = correlation
            line['Ratio Variance'] = ratioMatch
            line['Light Overlap Intensity'] = overlap1
            line['Heavy Overlap Intensity'] = overlap2
            line['Light Top-Signal Intensity'] = lightTopSignal
            line['Heavy Top-Signal Intensity'] = heavyTopSignal
            #line['C12 XIC Ratio'] = renderC12XIC(xics)
            #line['XIC Ratios'] = renderXICs(xics)
            line['Light Features'] = lightFeatures
            line['Heavy Features'] = heavyFeatures
            
            if isinstance(overlap1, float):
                line['Direct Ratio (H/L)'] = heavy / light if light else 'NaN'
                line['Overlap Ratio (H/L)'] = overlap2 / overlap1 if overlap1 else 'NaN'
                line['Light XIC Intensity'] = totalLightAOC
                line['Heavy XIC Intensity'] = totalHeavyAOC
                line['XIC Ratio (H/L)'] = totalHeavyAOC / totalLightAOC if totalLightAOC else 'NaN'
                line['Top-Signal Ratio (H/L)'] = heavyTopSignal / lightTopSignal if lightTopSignal and lightTopSignal != '-' else 'Nan"'
            else:
                line['Light XIC Intensity'] = '-'
                line['Heavy XIC Intensity'] = '-'
                line['Direct Ratio (H/L)'] = '-'
                line['Overlap Ratio (H/L)'] = '-'
                line['XIC Ratio (H/L)'] = '-'
                line['Top-Signal Ratio (H/L)'] = '-'
    
            output.write(line)
            
    output.close()        
    
    return outputFile
    