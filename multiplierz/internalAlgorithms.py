from math import floor, ceil, sqrt
from collections import defaultdict, deque
import gzip

import sys

try:
    from numpy import average
except ImportError:
    # Someday, hopefully, numpy will no longer be a multiplierz dependency.        
    def average(xs):
        return sum(xs) / float(len(xs))




def floatrange(start, stop, step = 1):
    cur = start
    while cur < stop:
        yield cur
        cur += step
        
        
        
        

def print_progress(step, total = None):
    """
    Shows a single-line changing string, to indicate progress without
    scrolling the terminal too much.
    
    Fiddly, though!  If total is used, step should always be numeric, and
    smaller than total.
    
    Furthermore, only use after printing a blank line (for it to overwrite)
    and print a newline after completing (to not append to the line used
    by this).
    
    Also won't work in the Wing terminal, I guess?
    """
    if total:
        sys.stdout.write('\r%d%% complete...   ' % int(step/float(total)))
    else:
        sys.stdout.write('\r%s...   ' % step)
    sys.stdout.flush()
        
    
        
    
def gzOptOpen(filename, mode = 'r'):
    '''
    Utility function for opening files that may be compressed
    via gzip ('fasta.gz' files.)
    '''
    
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, mode = mode)
    else:
        return open(filename, mode = mode)        
        
def unzip(thing): return [list(t) for t in zip(*thing)]

def typeInDir(directory, ext, recursive = False):
    import os
    if not recursive:
        return [os.path.join(directory, x) for x in
                os.listdir(directory) if x.lower().endswith(ext.lower())]
    else:
        output = []
        for path, subdirs, files in os.walk(directory):
            for filename in files:
                if filename.lower().endswith(ext.lower()):
                    output.append(os.path.join(path, filename))
        return output



def retryTimes(function, times, exceptions):
    for _ in xrange(times-1):
        try:
            result = function()
            return result
        except exceptions:
            continue
    # So that it properly throws an exception
    # the last time.
    return function()

def multisplit(thing, delimiters):
    fullout = []
    partout = []
    for char in (x for x in thing):
        if char in delimiters:
            if partout:            
                fullout.append(''.join(partout))
            partout = []
        else:
            partout.append(char)
    fullout.append(''.join(partout))
    return fullout

def inPPM(ppm, mA, mB):
    return abs(mA - mB) < (max((mA, mB))/1000000.0)*ppm
    
    
def parseClassicSpectrumDesc(desc):
    words = desc.split('.')
    numbers = words[-1].split('|')

    if len(numbers) > 20:
        plex = 10
    elif len(numbers) > 16:
        plex = 8
    elif len(numbers) > 12:
        plex = 6
    elif len(numbers) > 8:
        plex = 4
    else:
        plex = None

    uncorrected = numbers[-1*(plex+3):-3]
    corrected = numbers[:plex]

    datafile = words[0] # Usually.
    scanNum = words[1] # Again, usually.

    mz = desc.split('-')[1]

    return {'datafile':datafile,
            'scan':scanNum,
            'mz':mz,
            'corrected':corrected,
            'uncorrected':uncorrected}



class ProximityIndexedSequence(object):
    def __init__(self, sequence, indexer = lambda x: x, 
                 binCount = None, binSequence = None,
                 dynamic = True, **kwargs):
     
        
        if dynamic:
            self.indexer = indexer
        
        if sequence:
            indSeq = [(indexer(x), x) for x in sequence]
            self.bins = {(max([x[1] for x in indSeq]), min([x[1] for x in indSeq])):indSeq}
        else:
            self.bins = {}
            
        self.dynamic = dynamic
        
        self.rebalance()

    def tidyBins(self):
        print "ProximityIndexedSequence.tidyBins() is deprecated, because it sounds silly."
        self.rebalance()
        
    def rebalance(self):
        if not self.bins:
            return
        
        sequence = sum(self.bins.values(), [])
        self.bins = {}
        sequence.sort(reverse = True)
        
        binCount = floor(sqrt(len(sequence)))
        binLength = ceil(len(sequence)/binCount)
        
        newBin = []
        while sequence:
            newvals = [sequence.pop()]
            while sequence and sequence[-1] == newvals[-1]:
                newvals.append(sequence.pop())
            newBin += newvals
            
            if len(newBin) >= binLength:
                maxVal = max(newBin)[0]
                minVal = min(newBin)[0]
                self.bins[minVal, maxVal] = newBin[:]
                newBin = []
        if newBin:
            maxVal = max(newBin)[0]
            minVal = min(newBin)[0]
            self.bins[minVal, maxVal] = newBin[:]
        
        
    def seal(self):
        self.dynamic = False
        try:
            del self.indexer
        except NameError:
            pass
    
    def add(self, value):
        assert self.dynamic
        
        if not self.bins:
            ind = self.indexer(value)
            self.bins[ind, ind] = [(ind, value)]
            return
        
        el = (self.indexer(value), value)
        try:
            key = next(x for x in self.bins.keys() if x[0] <= el[0] <= x[1])
            self.bins[key].append(el)
        except StopIteration:
            # Add is often used to iteratively build the list; tidyBins()
            # *ought* to be called if things are being used correctly, but
            # its still faster if the bins are made in some reasonable size.
            maxbinlimit = max([x[1] for x in self.bins.keys()])
            minbinlimit = min([x[0] for x in self.bins.keys()])
            binwidth = average([x[1] - x[0] for x in self.bins.keys()])
            if not binwidth:
                binwidth = 1
            if el[0] > maxbinlimit:
                if maxbinlimit + binwidth > el[0]:
                    newbin = maxbinlimit, maxbinlimit + binwidth
                else:
                    while maxbinlimit + binwidth <= el[0]:
                        binwidth *= 2
                    newbin = maxbinlimit, maxbinlimit + binwidth
            elif el[0] < minbinlimit:
                if minbinlimit - binwidth < el[0]:
                    newbin = minbinlimit - binwidth, minbinlimit
                else:
                    while minbinlimit - binwidth >= el[0]:
                        binwidth *= 2
                    newbin = minbinlimit - binwidth, minbinlimit
            else:
                raise Exception, "Add new bin error."
            
            self.bins[newbin] = [el]
            
    
    def remove(self, value):
        assert self.dynamic
        
        index = self.indexer(value)
        key = next(x for x in self.bins.keys() if x[0] <= index <= x[1])
        del self.bins[key][self.bins[key].index((index, value))]
        
        if not self.bins[key]:
            del self.bins[key]
        elif index in key:
            indexes = [x[0] for x in self.bins[key]]
            newKey = min(indexes), max(indexes)
            self.bins[newKey] = self.bins[key]
            del self.bins[key]
            
        #try:
            #del self.bins[key][self.bins[key].index((index, value))]
        #except ValueError:
            #print "Failed removing %s, performing correction." % value
            #for key in self.bins.keys():
                #self.bins[key] = [x for x in self.bins[key] if x != (index, value)]    
        
        #for key in self.bins.keys():
            #if not self.bins[key]:
                #del self.bins[key]
                
    def __getitem__(self, index):
        # Which option is actually better seems complicated.  So far going with the least-code option.
        #try:
            #key = next((x for x in self.bins.keys() if x[0] <= index <= x[1]))
        #except StopIteration:
        key = min(self.bins.keys(), key = lambda x: min(abs(x[0] - index), abs(x[1] - index)))

        return min(self.bins[key], key = lambda x: abs(x[0] - index))[1]
    
    def asList(self):
        return [x[1] for x in sum(self.bins.values(), [])]
    
    def __iter__(self):
        for key in self.bins.keys():
            for index, thing in self.bins[key]:
                yield thing
                
    def returnRange(self, begin, stop):
        #raise NotImplementedError, "This has to be fixed!"
        #keys = [k for k in self.bins.keys() if begin <= k[0] <= stop or begin <= k[1] <= stop]
        keys = [k for k in self.bins.keys() 
                if (k[0] >= begin and k[1] <= stop) # Key range is within given range.
                or (k[0] <= begin and k[1] >= stop) # Given range is within key range.
                or (begin <= k[1] and stop >= k[1]) # Given range overlaps end of key range.
                or (stop >= k[0] and begin <= k[0])] # Given range overlaps beginning of key range.
        
        output = []
        for key in keys:
            if begin <= key[0] and stop >= key[1]:
                output += [x[1] for x in self.bins[key]]
            else:
                output += [x[1] for x in self.bins[key] if begin <= x[0] <= stop]
                
        return output



        
    
    
class NaiveProximitySequence(object):
    def __init__(self, sequence, indexer = lambda x: x):
        self.sequence = [(indexer(x), x) for x in sequence]
        self.indexer = indexer
    
    def __getitem__(self, thing):
        return min(self.sequence, key = lambda x: abs(thing - self.indexer(x)))[1]
    
    def remove(self, thing):
        index = self.indexer(thing)
        del self.sequence[self.sequence.index([x for x in self.sequence if x[0] == index][0])]
        
    def asList(self):
        return zip(*self.sequence)[1]
    
    
    
    
    

def collectByCriterion(data, criterion, splitby = None):
    """
    Utility function to collect a sequence of data objects based on
    a key computed by a "criterion" function.  The result is a dict
    of lists where all objects that evaluate to the same key are
    entries in the corresponding list.
    """
    output = defaultdict(list)
    for datum in data:
        key = criterion(datum)
        if isinstance(key, basestring):
            key = key.strip()
        if splitby:
            keys = key.split(splitby)
            for key in keys:
                output[key].append(datum)
        else:
            output[key].append(datum)
    return dict(output)
    




def centroid(scan):
    """
    Centroids profile-mode data given in [(M/Z, intensity)] format.
    """
    
    threshold = average(zip(*scan)[1])
    peaks = []
    peak = []
    for pt in scan:
        if pt[1] > threshold:
            peak.append(pt)
        elif peak:
            centroid = average(zip(*peak)[0], weights = zip(*peak)[1]), max(zip(*peak)[1])
            peaks.append(centroid)
            peak = []
    return peaks



    
# Tolerances widened for elegance and caution.
isotopicRatios = {(0, 1) : (0.5, 16.2),
                  (1, 2) : (1.0, 9.7),
                  (2, 3) : (1.4, 20.7),
                  (3, 4) : (1.7, 21.2),
                  (4, 5) : (2.0, 30.0),
                  (5, 6) : (2.4, 25.2),
                  (6, 7) : (2.6, 24.0),
                  (7, 8) : (2.9, 22.2),
                  (8, 9) : (3.0, 20.0),
                  (9, 10) : (3.4, 17.0),
                  (10, 11) : (3.6, 15.0)
                  }
                  
#from collections import defaultdict
#isotopicRatios = defaultdict(lambda: (0, 100))
    
    
    # Previously scan_features_opt
def peak_pick(scan, tolerance = 0.005, max_charge = 8, min_peaks = 3, correction = None,
              cleanup = False, recover_peaks = True,
              enforce_isotopic_ratios = True):
    """
    Scans a scan and gives back a by-charge dict of lists of isotopic sequences
    found in the scan, as well as a list of the unassigned peaks.
    
    tolerance - MZ range two points can be that are the "same" mass.
    max_charge - Maximum charge of peptides that are searched for.
    min_peaks - Minimum isotopic peaks required for an isotopic feature to be recorded.
    correction - Advanced feature; recalibration factor used for repeated calls on the same file.
    """

    if len(scan) > 10000:
        if all([(scan[i+1][0] - scan[i][0]) < 0.3 for i in range(0, 100)]):
            raise NotImplementedError, "Called scan_features on a profile-mode spectrum."

    if enforce_isotopic_ratios:
        global isotopicRatios
    else:
        isotopicRatios = defaultdict(lambda: (-1000000, 1000000))
    lowC13Rat, highC13Rat = isotopicRatios[0, 1]
    
    chargeFractions = [(x, 1.0/x) for x in range(1, max_charge+1)]
    chargeFracDict = dict(chargeFractions)
    scan.sort()
    
    envelopes = defaultdict(list)
    unassigned = []
    activePts = {}
    activeSets = defaultdict(deque)
    for pt in scan:
        # Search proceeds in one sweep through the scan in order of M/Z.
        pmz, pint = pt[:2]
        accounted = False

        # A given point is first tested against the set of currently open
        # isotopic sequences. If the point is beyond the highest possible next
        # MZ of a sequence, that sequence is closed (recorded to output.) If
        # it matches the expected mz of a sequence, and it's intensity matches
        # the expected ratio with the previous intensity in the sequence, it
        # is added to that sequence and the algorithm proceeds to the next
        # point.
        for chg, chargeSet in activeSets.items():
            # Chargeset is ordered by MZ, by construction, so only the first
            # (within-range) expected-point has to be checked.
            
            while chargeSet and chargeSet[0]['next'] + tolerance < pmz:
                if chargeSet:
                    finishedSet = chargeSet.popleft()
                    if len(finishedSet['envelope']) >= min_peaks:
                        envelopes[chg].append(finishedSet['envelope'])
                    else:
                        #unassigned += finishedSet['envelope']
                        # Give these peaks a second chance.
                        if chg < max_charge and recover_peaks:
                            for oldPt in finishedSet['envelope']:
                                possibleNexts = [(z, oldPt[0] + cF) for z, cF in chargeFractions if oldPt[0] + cF > pmz]
                                if possibleNexts:
                                    activePts[oldPt] = possibleNexts

            if chargeSet and chargeSet[0]['next'] - pmz < tolerance: # Too-low has already been eliminated?
                iso = chargeSet[0]
                isoLen = len(iso['envelope'])
                lowRat, highRat = isotopicRatios[isoLen-1, isoLen]
                
                if lowRat <= (iso['prevInt'] / pint) <= highRat:
                #if True:
                    iso['envelope'].append(pt)
                    nextMZ = pmz + chargeFracDict[chg]
                    iso['next'] = nextMZ
                    iso['prevInt'] = pint
    
                    # Updated envelope goes to the back of the queue for its charge.
                    chargeSet.rotate(-1)
    
                    accounted = True
                    break

            
        if accounted: # If it was put in a pre-existing isotopic sequence.
            continue


        # "Active points" are single previously-seen points that were not matched
        # to an isotopic sequence.  They are kept separate from isotopic sequences
        # since their charge is indeterminate, so there's multiple possible-next-MZs.
        # If a point matches to an active point, they're put together in a new 
        # isotopic sequence.
        for actPt in sorted(activePts.keys(), key = lambda x: x[1], reverse = True):
            if actPt[0] + 1 + tolerance < pmz:
                del activePts[actPt]
                unassigned.append(actPt)
            #elif (actPt[1]*10 > pint) and (pint*10 > actPt[1]):
            elif lowC13Rat <= (actPt[1] / pint) <= highC13Rat:
                    charge = next((charge for (charge, mz) in activePts[actPt] if
                                   abs(mz - pmz) < tolerance), None)
                    if charge:
                        nextMZ = pmz + chargeFracDict[charge]
                        newActiveSet = {'envelope':[actPt, pt],
                                        'charge':charge,
                                        'next':nextMZ,
                                        'prevInt':pint}
                        activeSets[charge].append(newActiveSet)
        
                        
                        del activePts[actPt]
                        
                        accounted = True
                        break
        
        # If the point didn't match any isotopic sequence or active point,
        # it is made into a new active point.
        if not accounted:
            # Optimized by precomputing the possible isotopic gaps.
            possibleNexts = [(chg, pmz + cF) for chg, cF in chargeFractions]
            activePts[pt] = possibleNexts
    
    
    
    newCorrection = []
    for charge, things in activeSets.items():
        for thing in things:
            if len(thing['envelope']) >= min_peaks:
                envelopes[charge].append(thing['envelope'])
            else:
                unassigned += thing['envelope']            
                
    unassigned += activePts.keys()
    


    # Cleanup phase- each envelope checks for unassigned peaks that
    # would correctly extend it.
    if cleanup:
        unProx = ProximityIndexedSequenceAgain(unassigned, indexer = lambda x: x[0],
                                               dynamic = False)        
        for charge, chgEnvelopes in envelopes.items():
            increment = 1.0/charge
            for envelope in chgEnvelopes:
                first = min(envelope)
                fmz = first[0] - increment
                
                #scale = isotopicRatios[0, 1]
                #befPt = next((x for x in unassigned if (abs(x[0] - fmz) < tolerance and
                                                        #scale[0] < (x[1]/first[1]) < scale[1])),
                             #None)
                befPt = unProx[fmz]
                #if befPt:
                scale = isotopicRatios[0, 1]
                if (abs(befPt[0] - fmz) < tolerance and
                    scale[0] < (befPt[1]/first[1]) < scale[1]):
                    envelope.insert(0, befPt)

    if correction != None:
        return dict(envelopes), unassigned, newCorrection
    else:
        return dict(envelopes), unassigned        
    
    
    
    
    
def peak_pick_PPM(scan, tolerance = 10, max_charge = 8, min_peaks = 3, correction = None,
                  cleanup = False, recover_peaks = True,
                  enforce_isotopic_ratios = True):
    """
    Scans a scan and gives back a by-charge dict of lists of isotopic sequences
    found in the scan, as well as a list of the unassigned peaks; this is the PPM-based
    version, which will be a bit slower.
    
    tolerance - MZ range two points can be that are the "same" mass.
    max_charge - Maximum charge of peptides that are searched for.
    min_peaks - Minimum isotopic peaks required for an isotopic feature to be recorded.
    correction - Advanced feature; recalibration factor used for repeated calls on the same file.
    """

    if len(scan) > 10000:
        if all([(scan[i+1][0] - scan[i][0]) < 0.3 for i in range(0, 100)]):
            raise NotImplementedError, "Called scan_features on a profile-mode spectrum."

    if enforce_isotopic_ratios:
        global isotopicRatios
    else:
        isotopicRatios = defaultdict(lambda: (-1000000, 1000000))
    lowC13Rat, highC13Rat = isotopicRatios[0, 1]
    
    chargeFractions = [(x, 1.0/x) for x in range(1, max_charge+1)]
    chargeFracDict = dict(chargeFractions)
    scan.sort()
    
    envelopes = defaultdict(list)
    unassigned = []
    activePts = {}
    activeSets = defaultdict(deque)
    for pt in scan:
        # Search proceeds in one sweep through the scan in order of M/Z.
        pmz, pint = pt[:2]
        accounted = False

        # A given point is first tested against the set of currently open
        # isotopic sequences. If the point is beyond the highest possible next
        # MZ of a sequence, that sequence is closed (recorded to output.) If
        # it matches the expected mz of a sequence, and it's intensity matches
        # the expected ratio with the previous intensity in the sequence, it
        # is added to that sequence and the algorithm proceeds to the next
        # point.
        for chg, chargeSet in activeSets.items():
            # Chargeset is ordered by MZ, by construction, so only the first
            # (within-range) expected-point has to be checked.
            
            #while chargeSet and chargeSet[0]['next'] + tolerance < pmz:
            while chargeSet and (chargeSet[0]['next'] < pmz
                                 and not inPPM(tolerance, chargeSet[0]['next'], pmz)):
                if chargeSet:
                    finishedSet = chargeSet.popleft()
                    if len(finishedSet['envelope']) >= min_peaks:
                        envelopes[chg].append(finishedSet['envelope'])
                    else:
                        #unassigned += finishedSet['envelope']
                        # Give these peaks a second chance.
                        if chg < max_charge and recover_peaks:
                            for oldPt in finishedSet['envelope']:
                                possibleNexts = [(z, oldPt[0] + cF) for z, cF in chargeFractions if oldPt[0] + cF > pmz]
                                if possibleNexts:
                                    activePts[oldPt] = possibleNexts

            #if chargeSet and chargeSet[0]['next'] - pmz < tolerance: # Too-low has already been eliminated?
            if chargeSet and inPPM(tolerance, chargeSet[0]['next'], pmz):
                iso = chargeSet[0]
                isoLen = len(iso['envelope'])
                lowRat, highRat = isotopicRatios[isoLen-1, isoLen]
                
                if lowRat <= (iso['prevInt'] / pint) <= highRat:
                #if True:
                    iso['envelope'].append(pt)
                    nextMZ = pmz + chargeFracDict[chg]
                    iso['next'] = nextMZ
                    iso['prevInt'] = pint
    
                    # Updated envelope goes to the back of the queue for its charge.
                    chargeSet.rotate(-1)
    
                    accounted = True
                    break

            
        if accounted: # If it was put in a pre-existing isotopic sequence.
            continue


        # "Active points" are single previously-seen points that were not matched
        # to an isotopic sequence.  They are kept separate from isotopic sequences
        # since their charge is indeterminate, so there's multiple possible-next-MZs.
        # If a point matches to an active point, they're put together in a new 
        # isotopic sequence.
        for actPt in sorted(activePts.keys(), key = lambda x: x[1], reverse = True):
            #if actPt[0] + 1 + tolerance < pmz:
            if actPt[0] + 1 < pmz and not inPPM(tolerance, actPt[0] + 1, pmz):
                del activePts[actPt]
                unassigned.append(actPt)
            #elif (actPt[1]*10 > pint) and (pint*10 > actPt[1]):
            elif lowC13Rat <= (actPt[1] / pint) <= highC13Rat:
                    charge = next((charge for (charge, mz) in activePts[actPt] if
                                   #abs(mz - pmz) < tolerance), None)
                                   inPPM(tolerance, mz, pmz)), None)
                    if charge:
                        nextMZ = pmz + chargeFracDict[charge]
                        newActiveSet = {'envelope':[actPt, pt],
                                        'charge':charge,
                                        'next':nextMZ,
                                        'prevInt':pint}
                        activeSets[charge].append(newActiveSet)
        
                        
                        del activePts[actPt]
                        
                        accounted = True
                        break
        
        # If the point didn't match any isotopic sequence or active point,
        # it is made into a new active point.
        if not accounted:
            # Optimized by precomputing the possible isotopic gaps.
            possibleNexts = [(chg, pmz + cF) for chg, cF in chargeFractions]
            activePts[pt] = possibleNexts
    
    
    
    newCorrection = []
    for charge, things in activeSets.items():
        for thing in things:
            if len(thing['envelope']) >= min_peaks:
                envelopes[charge].append(thing['envelope'])
            else:
                unassigned += thing['envelope']            
                
    unassigned += activePts.keys()
    


    ## Cleanup phase- each envelope checks for unassigned peaks that
    ## would correctly extend it.
    if cleanup:
        raise NotImplementedError

    if correction != None:
        return dict(envelopes), unassigned, newCorrection
    else:
        return dict(envelopes), unassigned            
    
    
    
# Just for debugging.
def deisochart(envelopes, unassigned):
    import matplotlib.pyplot as pyt
    from random import seed, shuffle
    
    seed(1)
    
    def lines(thing, charge = None, **etc):
        pyt.vlines(zip(*thing)[0], [0]*len(thing), zip(*thing)[1], **etc)
        if charge:
            for pt in thing:
                pyt.text(pt[0], pt[1], str(charge))
    
    lines(unassigned, color = 'b')
    length = len(sum(envelopes.values(), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in envelopes.items():
        for env in chgEnvs:
            label = "%s(%s)" % (i, chg)
            lines(env, charge = label, color = colors[i], linewidth = 3)
            i += 1
    
    pyt.show()
    
    

def deisocompare(envelopes1, unassigned1, envelopes2, unassigned2):
    import matplotlib.pyplot as pyt
    from random import seed, shuffle
    
    seed(1)
    
    def lines(thing, charge = None, upness = 1, **etc):
        pyt.vlines(zip(*thing)[0], [0]*len(thing), 
                   [x*upness for x in zip(*thing)[1]], **etc)
        if charge:
            for pt in thing:
                pyt.text(pt[0], upness * pt[1], str(charge))
    
    lines(unassigned1, color = 'b')
    length = len(sum(envelopes1.values(), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in envelopes1.items():
        for env in chgEnvs:
            label = "%s(%s)" % (i, chg)
            lines(env, charge = label, color = colors[i], linewidth = 3)
            i += 1
            
            
    lines(unassigned2, upness = -1, color = 'b')
    length = len(sum(envelopes2.values(), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in envelopes2.items():
        for env in chgEnvs:
            label = "%s(%s)" % (i, chg)
            lines(env, charge = label, upness = -1, color = colors[i], linewidth = 3)
            i += 1
    
    
    pyt.show()  
    
    
        

        
        
        
        
        
        