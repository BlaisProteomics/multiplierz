from multiplierz.mzAPI import mzFile
from multiplierz.mzReport import reader
from collections import defaultdict
from multiplierz.internalAlgorithms import ProximityIndexedSequence, collectByCriterion
from multiplierz import vprint

import pickle
import sqlite3
import os
import base64











class idLookup():
    def __getitem__(self, thing):
        return thing


# Returns (scan number -> feature number, scan number, feature number -> peak list) dict tuple.
def getScanFeatureDicts(datafile, featureData, absScanFeatures = False):
    scans = datafile.scan_info(0, 9999999)
    
    if absScanFeatures:
        ms1ScanLookup = idLookup()
    else:
        ms1ScanLookup = dict(enumerate([x[2] for x in scans if x[3] == 'MS1']))
        
    ms2toms1s = {}
    currentMS1 = None
    for scan in scans:
        if scan[3] == 'MS1':
            currentMS1 = scan[2]
        else:
            ms2toms1s[scan[2]] = currentMS1
    
    scanToFeature = defaultdict(list)
    scanFeatureToPeaks = {}
    featureToMS1s = {}
    for featureNumber, feature in featureData:
        featureToMS1s[featureNumber] = feature.scans
        
        for relativeMs1, peaks in feature.regions:
            absMS1 = ms1ScanLookup[relativeMs1]
            scanToFeature[absMS1].append(featureNumber)
            scanFeatureToPeaks[absMS1, featureNumber] = peaks
        


    return dict(scanToFeature), scanFeatureToPeaks, featureToMS1s

def getMS2FeatureDict(datafile, featureData, absScanFeatures = False):
    scans = datafile.scan_info(0, 9999999)
    ms2toms1s = {}
    currentMS1 = None
    for scan in scans:
        if scan[3] == 'MS1':
            currentMS1 = scan[2]
        else:
            ms2toms1s[scan[2]] = currentMS1    
    
    if absScanFeatures:
        ms1ScanLookup = idLookup()
    else:
        ms1ScanLookup = dict(enumerate([x[2] for x in scans if x[3] == 'MS1']))    
    
    scanToFeature = defaultdict(list)    
    for featureNumber, feature in enumerate(featureData):
        for relativeMs1, peaks in feature.regions:
            absMS1 = ms1ScanLookup[relativeMs1]
            scanToFeature[absMS1].append(featureNumber)
            scanFeatureToPeaks[absMS1, featureNumber] = peaks
            
    featureToMS2s = defaultdict(list)
    for index, ms2 in [(x[2], x) for x in scans if x[3] == 'MS2']:
        ms1 = ms2toms1s[index]
        features = scanToFeature[ms1]
        for featureNum, feature in [(x, featureData[x]) for x in features]:
            for scanNum, peaks in feature.regions:
                if ms1ScanLookup[scanNum] == ms1:
                #if scanNum == ms1:
                    if any([abs(x[0] - ms2[1]) < 0.05 for x in peaks]):
                        featureToMS2s[featureNum].append(ms2[2])
                    break    
                
    return dict(featureToMS2s)

mediumR = "Label:13C(6)"
mediumK = "Label:2H(4)"
heavyK = "Label:13C(6)15N(2)"
heavyR = "Label:13C(6)15N(4)"   

def featureToPSM(resultFile, featureData, groupSILAC = False):
    results = reader(resultFile)
    if 'Feature' not in results.columns:
        raise IOError("Not a feature-annotated file!")
    
    featureToPSMs = defaultdict(list)
    if groupSILAC:
        for psm in results:
            mods = psm['Variable Modifications']
            if mods == None: mods = []
            
            isHeavy = heavyK in mods or heavyR in mods
            isMedium = (mediumK in mods or mediumR in mods) and not isHeavy
            isLight = not (isHeavy or isMedium)
            
            if isLight:
                if not psm['Light Features']: continue
                features = str(psm['Light Features']).split(';')
            elif isMedium:
                if not psm['Medium Features']: continue
                features = str(psm['Medium Features']).split(';')
            else:
                if not psm['Heavy Features']: continue
                features = str(psm['Heavy Features']).split(';')
            
            for feature in features:
                feature = int(float(feature))
                featureToPSMs[feature].append(psm)
    else:
        for psm in results:
            try:
                featureToPSMs[int(float(psm['Feature']))].append(psm)
            except ValueError:
                pass
                
    
    return dict(featureToPSMs)
    
    


def save_feature_database(features, outputfile, overwrite = None):
    """
    Saves a SQLite-mode feature database. Result file will have the
    extension '.features' .
    """
    
    if os.path.exists(outputfile):
        if overwrite != True:
            os.remove(outputfile)
        else:
            raise IOError("Target file %s already exists!" % outputfile)
    
    conn = sqlite3.connect(outputfile)
    cur = conn.cursor()
    
    createTable = "CREATE TABLE features(ind int, mz real, startscan int, endscan int, data text)"
    cur.execute(createTable)

    vprint("Created table.")
    for index, feature in enumerate(features):
        mz = feature.mz
        startscan, endscan = feature.scanrange
        featureData = base64.b64encode(pickle.dumps(feature, protocol = 2))
        addFeature = ('INSERT INTO features VALUES (%s, %s, %s, %s, "%s")' 
                      % (index, mz, startscan, endscan, featureData.decode()))
        cur.execute(addFeature)
    
        if index % 100 == 0:
            conn.commit()
            
    #print("TEST MODE")
    #sidechannel = open(outputfile + "SIDECHANNEL.pickle", 'wb')
    #pickle.dump(list(enumerate(features)), sidechannel)
    #sidechannel.close()
    
    vprint("Indexing...")
    createIndex = "CREATE INDEX mzindex ON features(mz, startscan)"
    cur.execute(createIndex)  
    vprint("Analyzing...")
    cur.execute("ANALYZE")
    
    vprint("Final SQLite commit...")
    conn.commit()
    
    conn.close()
    
    
    
    

def newMarshal(x):
    return pickle.loads(base64.b64decode(x))
def oldMarshal(x):
    return pickle.loads(x)

    
# Doing feature data manipulation proper-like, now.
# For now just a wrapper over loading the whole standard pickle
# object, naively; later a SQLite interface?
class FeatureInterface(object):    
    """
    A wrapper over either a pickled list of features or a SQLite database
    containing a feature table.  The indexing operator returns a feature
    according to its index value (e.g., order of discovery, and how they're
    named in annotated search results.)  For example
    
    > data = FeatureInterface('featurefile.features')
    > data[100]
    < Feature object at  ...>
    
    Commands mz_range and scan_range return iterators over features
    that fall within the specified ranges.
    """
    def __init__(self, filename):
        self.filename = filename
        assert os.path.exists(filename)
        if filename.lower().endswith('featurepickle'):
            vprint("Legacy mode enabled.")
            self.data = pickle.load(open(filename))
            self.mode = 'pickle'
        else:
            self.connection = sqlite3.connect(filename)
            self.data = self.connection.cursor()
            self.mode = 'sql'
            
            self.decoder = newMarshal
            #self.data.execute("SELECT data FROM features WHERE ind=1")
            #testfeature = str(self.data.fetchone()[0])
            #if '\n' in testfeature:
                ## Old non-base64 encoded feature file!
                #self.decoder = oldMarshal
            #else:
                #self.decoder = newMarshal
            
            
        
    def __getitem__(self, index):
        if self.mode == 'pickle':
            assert isinstance(index, int)
            return self.data[index]
        
        command = "SELECT data FROM features WHERE ind=%s" % index
        self.data.execute(command)
        feature = self.decoder(str(self.data.fetchone()[0]))
        return feature
        
    def __iter__(self):
        assert self.mode == 'sql'
        self.data.execute('SELECT ind, data FROM features')
        sublist = self.data.fetchmany()
        while sublist:
            for i_f in sublist:
                yield i_f[0], self.decoder(str(i_f[1]))
            sublist = self.data.fetchmany()        
    
    def mz_range(self, start_mz, end_mz):
        assert self.mode == 'sql', 'Requires SQLite mode!'
        
        command = "SELECT ind, data FROM features WHERE mz >= %s AND mz <= %s" % (start_mz, end_mz)
        self.data.execute(command)
 
        sublist = self.data.fetchmany()
        while sublist:
            for i_f in sublist:
                yield i_f[0], self.decoder(str(i_f[1]))
            sublist = self.data.fetchmany()
    
    def scan_range(self, start_scan, end_scan):
        # Returns features fully enclosed in the specified scan range.
        command = "SELECT ind, data FROM features where start_scan >= %s AND end_scan <= %s" % (start_scan, end_scan)
        sublist = self.data.fetchmany()
        while sublist:
            for i_f in sublist:
                yield i_f[0], self.decoder(str(i_f[1]))
            sublist = self.data.fetchmany()                
        
    def close(self):
        if self.mode == 'sql':
            self.conection.close()
    


class FeatureInterface_preload(object):
    def __init__(self, filename):
        source = FeatureInterface(filename)
        self.features = source.mz_range(0, 5000)
        
        
    
            

