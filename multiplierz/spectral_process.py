from collections import defaultdict
from multiplierz.internalAlgorithms import average
from multiplierz import protonMass
from collections import deque, Iterator
from numpy import std

def centroid(scan, threshold = None):
    """
    Centroids profile-mode data given in [(M/Z, intensity)] format.
    """
    from internalAlgorithms import average
    
    if not scan:
        return scan
    
    if not threshold:
        assert not isinstance(scan, Iterator), ("centroid() requires explicit "
                                                "threshold value if given an "
                                                "Iterator argument.")

        threshold = average(list(zip(*scan))[1])
        
    peaks = []
    peak = []
    for pt in scan:
        if pt[1] > threshold:
            peak.append(pt)
        elif peak:
            centMZ = average(list(zip(*peak))[0], weights = list(zip(*peak))[1])
            centInt = max(list(zip(*peak))[1])
            if len(pt) == 4:
                centNoise = average(list(zip(*peak))[2], weights = list(zip(*peak))[1])
                centCharge = max(list(zip(*peak))[3])
                peaks.append((centMZ, centInt, centNoise, centCharge))
            else:
                peaks.append((centMZ, centInt))
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
                  (10, 11) : (3.6, 15.0),
                  (11, 12) : (0.0, 0.0) # To gracefully close an endless envelope.
                  }
# Tolerances widened even further, for chlorine atoms.
isotopicRatios_permissive = {(0, 1) : (0.5, 16.2),
                             (1, 2) : (1.0, 9.7),
                             (2, 3) : (1.0, 20.7),
                             (3, 4) : (1.0, 21.2),
                             (4, 5) : (1.0, 30.0),
                             (5, 6) : (1.5, 25.2),
                             (6, 7) : (1.5, 24.0),
                             (7, 8) : (2.0, 22.2),
                             (8, 9) : (3.0, 20.0),
                             (9, 10) : (3.4, 17.0),
                             (10, 11) : (3.6, 15.0),
                             (11, 12) : (0.0, 0.0)
                           }
# Dev note: "First, do no harm"- all peaks from the input scan should be
# represented somewhere in the output (and no new peaks, obviously.)
def peak_pick(scan, tolerance = 0.005, max_charge = 8, min_peaks = 3, correction = None,
              cleanup = False, recover_peaks = True,
              enforce_isotopic_ratios = True):
    """
    Scans a scan and gives back a by-charge dict of lists of isotopic sequences
    found in the scan, as well as a list of the unassigned peaks.

    The input scan must be centroided!

    tolerance - MZ range two points can be that are the "same" mass.
    max_charge - Maximum charge of peptides that are searched for.
    min_peaks - Minimum isotopic peaks required for an isotopic feature to be recorded.
    correction - Advanced feature; recalibration factor used for repeated calls on the same file.
    """

    if enforce_isotopic_ratios == True:
        global isotopicRatios
    elif enforce_isotopic_ratios == 'permissive':
        isotopicRatios = isotopicRatios_permissive
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
        for chg, chargeSet in list(activeSets.items()):
            # Chargeset is ordered by MZ, by construction, so only the first
            # (within-range) expected-point has to be checked.

            while chargeSet and chargeSet[0]['next'] + tolerance < pmz:
                if chargeSet:
                    finishedSet = chargeSet.popleft()
                    if len(finishedSet['envelope']) >= min_peaks:
                        envelopes[chg].append(finishedSet['envelope'])
                    else:
                        # The envelope hasn't managed to accumulate enough
                        # peaks to meet requirements, but some of it's
                        # constituent peaks could still go on to be members
                        # of other, lower-charge envelopes.
                        if recover_peaks:
                            if chg < max_charge:
                                for oldPt in finishedSet['envelope']:
                                    possibleNexts = [(z, oldPt[0] + cF) for z, cF in chargeFractions if oldPt[0] + cF > pmz]
                                    if possibleNexts:
                                        activePts[oldPt] = possibleNexts
                                    else:
                                        unassigned.append(oldPt)
                            else:
                                unassigned += finishedSet['envelope']

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
        for actPt in sorted(list(activePts.keys()), key = lambda x: x[1], reverse = True):
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
    for charge, things in list(activeSets.items()):
        for thing in things:
            if len(thing['envelope']) >= min_peaks:
                envelopes[charge].append(thing['envelope'])
            else:
                unassigned += thing['envelope']            

    unassigned += list(activePts.keys())



    # Cleanup phase- each envelope checks for unassigned peaks that
    # would correctly extend it.  This *usually* does not find anything new.
    if cleanup:
        unProx = ProximityIndexedSequenceAgain(unassigned, indexer = lambda x: x[0],
                                               dynamic = False)        
        for charge, chgEnvelopes in list(envelopes.items()):
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
    
    
    
    
    
def deisotope_reduce_scan(spectrum, *peak_pick_args, **peak_pick_kwargs):
    """
    Deisotopes the given mass spectrum, and replaces isotopically pure peaks
    of determined charge with imputed peaks of charge +1. This often improves
    search scores, particularly when using a database search engine which
    does not account for all possible charge states.
    
    For best results, the scan mass tolerance should be below 0.01 Daltons.
    """
    global protonMass
    
    chargeEnvelopes, peaks = peak_pick(spectrum,
                                       *peak_pick_args,
                                       **peak_pick_kwargs)
    for charge, envelopes in list(chargeEnvelopes.items()):
        for envelope in envelopes:
            mz, ints = envelope[0]
            if charge > 1:
                mass = (mz * charge) - (protonMass * (charge-1))
            else:
                mass = mz
            peaks.append((mass, ints))
    return peaks

def deisotope_scan(spectrum, *args, **kwargs):
    """
    Deisotopes the given mass spectrum, but leaves isotopically pure peaks
    in place.  This will often improve search score (due to unscored
    isotopic peaks being absent.)
    
    For best results, the scan mass tolerance should be below 0.01 Daltons.
    """
     
    chargeEnvelopes, peaks = peak_pick(spectrum, *args, **kwargs)
    for charge, envelopes in list(chargeEnvelopes.items()):
        for envelope in envelopes:
            peaks.append(envelope[0])
    return peaks           

def get_monoisotopic(spectrum, *args, **kwargs):
    """
    Performs peak-picking and returns only known-monoisotopic peaks.
    """
    
    charge_envelopes, _ = peak_pick(spectrum, *args, **kwargs)
    peaks = []
    for charge, envelopes in list(charge_envelopes.items()):
        for envelope in envelopes:
            peaks.append((charge, envelope[0]))
    return peaks
        



def top_n_peaks(spectrum, N):
    """
    Takes the most-intense N peaks of the given scan, and discards the rest.
    """
    N = int(N) if N else 0
    return sorted(spectrum, key = lambda x: x[1], reverse = True)[:int(N)]

def exclusion_radius(spectrum, exclusion):
    """
    For every peak in the spectrum, in order of most-to-least intense, takes
    that peak into the output while discarding all other peaks within a given
    exclusion "radius" (in Daltons.)  Continues until all peaks have been taken
    or discarded.
    """
    if not exclusion or not float(exclusion):
        return spectrum
    exclusion = float(exclusion)
    
    acc = []
    while spectrum:
        peak = max(spectrum, key = lambda x: x[1])
        acc.append(peak)
        spectrum = [x for x in spectrum if not abs(x[0] - peak[0]) < exclusion]
    
    return acc    

def signal_noise(spectrum, minSN):
    """
    Filters the spectrum by a signal-to-noise threshold, by finding the subset
    of lowest-intensity peaks that match a specified minimum signal-to-noise
    ratio and discarding those.
    """
    
    # Ought to do a binary search!
    
    if not minSN or not float(minSN):
        return spectrum
    
    minSN = float(minSN)
    spectrum.sort(key = lambda x: x[1])
    for i in range(0, len(spectrum)):
        ints = [x[1] for x in spectrum[i:]]
        SN = average(ints) / std(ints)
        if SN > minSN:
            return spectrum[i:]
    
    return spectrum

def intensity_threshold(spectrum, threshold):
    """
    Discards all peaks in the spectrum that have an intensity below
    the specified threshold.
    """
    threshold = float(threshold)
    return [x for x in spectrum if float(x[1]) > threshold]

def mz_range(spectrum, range):
    """
    Discards all peaks in the spectrum that are beyond the specified
    MZ range.
    """
    
    if isinstance(range, str):
        start, stop = [int(x) for x in range.split('-')]
    else:
        start, stop = range
    return [x for x in spectrum if start <= x[0] <= stop]


def recalibrate(spectra, constant = 0, slope = 1):
    """
    For each spectrum in spectra, shift the MZ of every point to:
    
    MZ' = (MZ * slope) + constant
    """
    
    for spectrum in spectra:
        yield [((x[0]*slope)+constant, x[1]) for x in spectrum]
    