from multiplierz.mzAPI import mzFile
import matplotlib.pyplot as pyt
import os
import numpy as np
import win32com.client
import matplotlib as mpl
import wx
#from multiplierz.internalAlgorithms import ProximityIndexedSequence
from collections import defaultdict

import time

import iminuit

__all__ = ['chargeTransform']

try:
    import MSO
except ImportError:
    print ("WARNING: Microsoft Office Object Library not found."
           "This is necessary for mzTransform's Powerpoint-slide"
           "output feature.")
    MSO = None


thresholdCoefficient = 1
minimumThreshold = 3

from itertools import dropwhile

from fractions import gcd


initialPeakRemoval = str(3)
initialRequiredPeaks = str(8)
initialPeakIterations = str(2)
initialMZRange = "300-2000"


try:
    from mzGUI import file_chooser
    mzgui = True
except ImportError:
    mzgui = False

initialTolerance = 1
finalTolerance = 1
binSize = 30
maxPossibleCharge = 200
minPossibleCharge = 2
requireLength = 8

from math import isnan

proton = 1.00727647



fileTypes = ['.txt', '.ppt', '.pptx', '.pgf', '.svgz', '.tiff', '.jpg', '.raw',
             '.jpeg', '.png', '.ps', '.svg', '.eps', '.rgba', '.pdf', '.tif']


def unzip(thing): return [list(t) for t in zip(*thing)]

def numberSequence(start, end, step):
    assert (start < end and step > 0) or (start > end and step < 0), "Invalid number sequence."

    foo = [start]
    while (foo[-1] + step) < end:
        foo.append(foo[-1] + step)

    return foo

def pairplot(thing):
    pyt.plot(list(zip(*thing))[0], list(zip(*thing))[1])

class EmptyMass(Exception):
    def __init__(self, value = None):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class SpectralMiss(Exception):
    def __init__(self, value = None):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class EmptySpectrum(Exception):
    def __init__(self, value = None):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def stripNone(sequence):
    for i in range(len(sequence)-1, 0, -1):
        if not sequence[i][0]: del sequence[i]
        else: break
    return list(dropwhile(lambda x: x[0] == None, sequence))

def line(xs, a, b):
    return np.array([x*a + b for x in xs])
def curve(xs, a, b, c, d):
    return np.array([(x**3)*a + (x**2)*b + x*c + d for x in xs])

def iterativeFitBreak(data, fitFunc = curve, depth = 2):
    """
    At several points we want to find "peaks" that are actually just protuberances
    above some underlying noise factor; to do this, we fit a curve to the data
    (which winds up mostly being fit to the noise) and split off the regions which
    lie above this line.
    
    Also, peaks often have sub-peaks (and those have sub-sub-peaks, etc), so we
    may apply this recursively.
    
    The output is a set of regions, the maximum points of which are the nominal
    "peak points."
    """
    
    
    
    if min(list(zip(*data))[1]) + 3 > max(list(zip(*data))[1]):
        return [data]
    elif len(data) < 4:
        return [data]
        

    
    mzs, intensities = unzip(data)
    avgInt = np.average(intensities)
    normInts = [x - avgInt for x in intensities]
    #pars, _ = curve_fit(fitFunc, mzs, normInts) # Will randomly and without warning return garbage!  Get rid of!!
    
    #filtered = [x if x > 0 else 0 for x in (np.array(normInts) - fitFunc(mzs, *pars))]
    #pyt.plot(mzs, fitFunc(mzs, *pars), color = 'b')
    
    def mFitFunc(a, b, c, d):
        return sum([abs((((mz**3)*a) + b*(mz**2) + c*mz + d) - intensity)
                    for mz, intensity in data])
    
    m = iminuit.Minuit(mFitFunc, a = 0.0, error_a = 0.1,
                       b = 0.0, error_b = 0.1,
                       c = 0.0, error_c = 0.1,
                       d = 0.0, error_d = 0.1, 
                       print_level = 0,
                       errordef = 1)
    
    m.migrad(1000)
    
    a, b, c, d = m.args
    fitCurve = [intensity - (((mz**3)*a) + b*(mz**2) + c*mz + d) for mz, intensity in data]
    filtered = [mz if mz > 0 else 0 for mz in fitCurve]
    
    
    if not all(filtered):
        for i in range(0, len(filtered)):
            if not filtered[i]: break
            else: filtered[i] = 0
    #else:
        ## Curve-fitting has obviously failed.
        #if fitFunc == line:
            #return []
        #else:
            #return iterativeFitBreak(data, fitFunc = line, depth = depth)
    
    peaks = []
    newPeak = []
    for corInt, dataPt in zip(filtered, data):
        if corInt == 0:
            if newPeak:
                peaks.append(newPeak)
                newPeak = []
        else:
            newPeak.append(dataPt)
    
    if not peaks:
        return [data]
    elif depth > 0 and peaks[0] != data:
        peaks = sum([iterativeFitBreak(xs, fitFunc, depth-1) for xs in peaks], [])
        return peaks
    else:
        return peaks
    
    


def obtainPeaksForMass(data, mass, tolerance, recorrect = False, requireConsistent = True):
    peaks = []
    charges = []
    lowerBound, upperBound = data[0][0], data[-1][0]
    for chg in range(minPossibleCharge, maxPossibleCharge):
        mz = (mass + (chg * proton)) / chg
        if not lowerBound < mz < upperBound: continue
        
        try:
            #peak = (x for x in data if abs(x[0] - mz) < tolerance and x not in peaks).next()
            peak = max([x for x in data if abs(x[0] - mz) < tolerance and x not in peaks],
                       key = lambda x: x[1])
        #except StopIteration:
        except ValueError:
            peaks.append(None)
            charges.append(None)
            continue

        peaks.append(peak)
        charges.append(chg)
        

    
    if not recorrect:
        if not peaks: raise SpectralMiss(str(mass))
        return peaks, charges
    else:
        return obtainPeaksForMass(data, refineMass(list(zip(peaks, charges)), mass),
                                  tolerance, recorrect = False)

def refineMass(chargePeaks, mass):
    hypotheticalMasses = []
    charges = [x for x in unzip(chargePeaks)[1] if x]
    minChg, maxChg = min(charges), max(charges)
    weightConst = (maxChg - minChg) / 100.0
    for peak, chg in chargePeaks:
    #for index in range(1, len(chargePeaks)):
        #peak, chg = chargePeaks[index]
        #prevPeak, prevChg = chargePeaks[index-1]
        if not peak: continue
        
        #before = 1.0 / prevChg
        #current = ((peak[0] / prevPeak[0]) - 1) / chg
        #error = abs(before - current)
        mz = peak[0]
        errSum = 0
        for otherPeak, otherChg in chargePeaks:
            if peak == otherPeak or not otherPeak: continue
            otherMz = otherPeak[0]
            otherChg = float(otherChg)
            errSum += (abs((chg / otherChg) - (otherMz / mz)) / (chg / otherChg)) ** 2
        
        weight = errSum + 0.0000001
        
        
        #chargeEstimate = mass / peak[0]
        chargeEstimate = mass / (peak[0] - proton)
        if not round(chargeEstimate) == chg:
            continue
        if abs(round(chargeEstimate) - chg) > 0.1:
            print("Distant.")
        hypotheticalMasses.append(((peak[0] * chg) - (chg * proton), weight))
    
    #assert hypotheticalMasses
    if not hypotheticalMasses: raise SpectralMiss
    #return sum(hypotheticalMasses) / len(hypotheticalMasses)
    return np.average(unzip(hypotheticalMasses)[0], weights = unzip(hypotheticalMasses)[1])
    
        
        
        
        
        
def validatePeaks(peakPoints):
    stripPeaks = stripNone(peakPoints)
    
    # Widget to account for invalid large gaps caused by the acquisition of some spurious
    # peak far away from the actual peak sequence.
    noneCount = 0
    acc = []
    for peak in stripPeaks:
        if peak[0]:
            noneCount = 0
            acc.append(peak)
        elif not peak[0] and noneCount < 5: # Random number; maximum gap penalty?
            noneCount += 1
            acc.append(peak)
    stripPeaks = acc      
        
    # Widget to account for the acquisition of peaks in between the real peak sequence;
    # using the fact that intensity of peaks almost always remains changes smoothly from
    # charge to charge.
    for index, _ in sorted(enumerate(stripPeaks), key = lambda x: x[1][1]):
        if index == 0 or index == len(stripPeaks)-1: 
            continue
        before, current, after = stripPeaks[index-1], stripPeaks[index], stripPeaks[index+1]
        if not (before[0] and current[0] and after[0]):
            continue
        
        if before[0][1] > current[0][1]*2 or after[0][1] > current[0][1]*2:
            stripPeaks[index] = (None, None)
    
    return stripPeaks
    
        
# Takes a list of mz-space peaks and a guess for the chargeless mass that
# would account for a subset of those peaks.  Necessarily heuristic in nature,
# so there's way too many magic numbers and fairly arbitrary parameters scattered
# about.
#
# There are two important ways the mass guess could be wrong; it could be some factor
# smaller than the true mass, or some factor larger.  If its a factor larger, then
# the predicted series of mz peaks will have lots of gaps were there aren't valid 
# true mz peaks to fill the spot; this is measured in the noneGapFraction.  In this
# case, the greatest common denomenator between the charges of the greatest two mz
# peaks that are found should indicate the guess' multiple of the true mass.
#
# If the guess is a factor smaller than the true mass, it is (in general, approximately)
# smaller by some factor which is a product of reasonably small primes.  Such
# primes are each checked to see if the guess multiplied by that prime yields a
# more successful result, and if one is found, it is recursed upon.
#
# Note that non-multiplicative inaccuracies in the guess are generally solved
# by refineMass(), so nothing here really addresses that very well.
def assignMassCharge(peakPoints, mass):
    peaks = list(zip(*obtainPeaksForMass(peakPoints, mass, finalTolerance)))
    
    stripPeaks = validatePeaks(peaks) # Removes not-found peaks from beginning and end.
    

    if not stripPeaks:
        raise EmptyMass(str(mass))
    
    noneGapFraction = len([x for x in stripPeaks if x[0] == None]) / float(len(stripPeaks))
    if noneGapFraction > 0.3:
        truePeaks = [x for x in stripPeaks if x[0] != None]
        #truePeaks.sort(key = lambda x: x[0][1], reverse = True)
        #first, second = truePeaks[:2]
        first = max([x for x in stripPeaks if x[0]], key = lambda x: x[0][1])
        firstIndex = truePeaks.index(first)
        second = min(enumerate(truePeaks), key = lambda x: abs(x[0] - firstIndex) if x[0] != firstIndex else 999999)[1]
        
        
        chargeDen = gcd(first[1], second[1])
        #counts = defaultdict(int)
        #for index in range(0, len(stripPeaks)-1):
            #current, after = stripPeaks[index], stripPeaks[index+1]
            #if current[1] and after[1]:
                #counts[gcd(current[1], after[1])] += 1
            
        
        if chargeDen == 1:
            if noneGapFraction > 0.5:
                raise SpectralMiss(str(noneGapFraction))
                    
            finalMass = refineMass(peaks, mass)
            finalPeaks = list(zip(*obtainPeaksForMass(peakPoints, finalMass, finalTolerance)))
            return [x for x in finalPeaks if x[0]], finalMass
        else:    
            return assignMassCharge(peakPoints, mass / chargeDen)
    
    else:
        # Check mass multiples of the next few primes to see if its skipping over peaks.
        for multiple in [2,3,5,7,11]:
            try:
                nextPeaks = list(zip(*obtainPeaksForMass(peakPoints, mass*multiple, finalTolerance)))
                nextPeaks = validatePeaks(nextPeaks)
                if (not nextPeaks) or len(nextPeaks)*2 < len(stripPeaks): continue
            except SpectralMiss:
                continue
            
            nextGapFraction = len([x for x in nextPeaks if x[0] == None]) / float(len(nextPeaks))
            
            if nextGapFraction <= noneGapFraction * 0.9:
                return assignMassCharge(peakPoints, mass*multiple)
            
            
    finalMass = refineMass(peaks, mass)
    finalPeaks = list(zip(*obtainPeaksForMass(peakPoints, finalMass, finalTolerance)))
    return [x for x in finalPeaks if x[0]], finalMass
                
        
        
 
    


# This is the centerpoint of the deconvolution algorithm.
# It happens in three stages:
# - First, peak regions are detected with iterativeFitBreak and various
#   thresholding things are done to try and winnow out the unimportant ones.
# - Second, a "potential mass spectrum" is made by multiplying every
#   peak by every possible charge.  This results in essentially the same
#   thing as the Mann algorithm output.  In order to be useful, these
#   potential-mass-peaks are binned over small mass increments, and the
#   total intensity of each bin is taken.
# - Third, bins are taken in order of highest intensity and their associated
#   masses are used as guesses for ouptut masses.  Lots of things are done
#   to refine and vet each guess; in particular, due to the periodic nature of
#   the Mann algorithm chart, no two guesses are allowed to be multiples/fractions
#   of each other.  For the rest of these conditions and procedures, see
#   assignMassCharge().
def assignAllMasses(data, removalArea = None, requiredLength = 0, 
                    peakIterations = 2, massRange = None):
    peaks = iterativeFitBreak(data, line, peakIterations)
    print(("Peaks found. %s" % time.clock()))
    
    if not peaks:
        print("No peaks detected! ")
        return []
    
    peakPoints = []
    for peak in peaks:
        peakPoints.append(max(peak, key = lambda x: x[1]))
        

    # Filter out very low intensity "peaks."  Should be replaced by a formal
    # signal-to-noise threshold, most likely.
    peakIntensities = [x[1] for x in peakPoints]
    avgInt = sum(peakIntensities) / len(peakIntensities)
    peakIntensities = [x for x in peakIntensities if x <= avgInt]
    avgInt2 = sum(peakIntensities) / len(peakIntensities)
    peakIntensities = [x for x in peakIntensities if x <= avgInt2]
    threshold = (sum(peakIntensities) / len(peakIntensities)) * thresholdCoefficient
    threshold = max([threshold, minimumThreshold])

    peakPoints = [x for x in peakPoints if x[1] >= threshold]
    
    # The above may still leave an overwhelming number of peaks.  Ideally would
    # take out peaks based on overcrowding in given mz regions.
    if len(peakPoints) > 300:
        peakPoints = sorted(peakPoints, key = lambda x: x[1], reverse = True)[:300]
        
    if not peakPoints:
        print("No peaks detected!")
        return []    
        
    peakPoints.sort(key = lambda x: x[1], reverse = True)
    
    lowBound, highBound = (min([x[0] for x in peakPoints]), 
                           max([x[0] for x in peakPoints]))
    
    #print "Projecting possible mass values. %s" % time.clock()
    potentialMassSpec = []
    for peak in peakPoints:
        mz = peak[0]
        intensity = peak[1]
        for charge in range(minPossibleCharge, maxPossibleCharge):
            potentialMassSpec.append(((mz*charge) - (charge*proton), intensity, mz, charge))
    
    #if not massRange:
    highMassBound = int(highBound*maxPossibleCharge)
    lowMassBound = 0
    #else:
        #highMassBound = massRange[1]
        #lowMassBound = massRange[0]
    
    
    potentialMassSpec.sort(key = lambda x: x[0])
    print((len(potentialMassSpec)))
    potentialMassSpec = (x for x in potentialMassSpec) # Scopes and closures!
    
    print(("Binning mass values. %s" % time.clock()))
    massBins = []
    low = 0
    nextEl = next(potentialMassSpec)
    for i in range(lowMassBound, highMassBound, binSize):
        if not i: 
            continue
        high = i
        
        massBin = []
        while nextEl[0] < high:
            massBin.append(nextEl)
            try:
                nextEl = next(potentialMassSpec)
            except StopIteration:
                massBin.append(nextEl)
                break
            
        if not massBin:
            continue

        avgMass = sum(unzip(massBin)[0]) / len(massBin)
        std = np.std([x[0] for x in massBin])
        core = [x for x in massBin if avgMass - std <= x[0] <= avgMass + std]
        avg2 = sum([x[0] for x in core]) / len(core)
        
        
        #intensity = sum(unzip(massBin)[1])
        massBins.append([massBin, avg2, 0])
        low = high

    binIntensities = [sum([x[1] for x in xs]) for xs,_,_ in massBins]
    #print sum(binIntensities)
        
    massBins[0][2] = binIntensities[0]+(binIntensities[1]/2)
    for i in range(1, len(massBins) - 1):
        massBins[i][2] = (binIntensities[i-1]/2)+binIntensities[i]+(binIntensities[i+1]/2)
    massBins[-1][2] = (binIntensities[-2]/2)+binIntensities[-1]

    massBins.sort(key = lambda x: x[2])
    
    peakPoints.sort(key = lambda x: x[0])
    chargeSequences = []
    masses = []
    while massBins and peakPoints:
        highPoint = massBins.pop()
        highBin, highMass, highIntensity = highPoint
        
        mass = highMass
        invalidMass = False
        for previousMass in masses:
            multiple = mass / previousMass
            divisible = previousMass / mass
            if abs(round(multiple) - multiple) < 0.01:
                #print "Multiple of previous mass."
                invalidMass = True
                break
            if abs(round(divisible) - divisible) < 0.01:
                #print "Fraction of previous mass."
                invalidMass = True
                break
        if invalidMass:
            continue

        try:
            chargePeaks, finalMass = assignMassCharge(peakPoints, mass)
        except EmptyMass as err:
            print(err)
            continue
        except (SpectralMiss, ZeroDivisionError) as err:
            continue
            
        for peak, _ in chargePeaks:
            try:
                peakPoints.remove(peak)
            except ValueError:
                pass
            if removalArea:
                peakPoints = [x for x in peakPoints if abs(peak[0] - x[0]) > removalArea]
            
        if len(chargePeaks) >= requiredLength and ((not massRange) or massRange[0] <= finalMass <= massRange[1]):
            chargeSequences.append((chargePeaks, finalMass))
        masses.append(finalMass)

      
    #print "Done. %s" % time.clock()
   
    return chargeSequences
                

##def mannEvaluation(chargeSequence):
    ##sequence, mass = chargeSequence
    
    #points = []
    #for index in range(1, len(sequence)):
        #if not sequence[index-1][1]: continue
        #i = sequence[index-1][1]
        #if not sequence[index][1]: continue
        #j = sequence[index][1] - i
        #imz = sequence[index-1][0][0]
        #jmz = sequence[index][0][0]
        
        
        #x = 1.0 / i
        #y = ((imz / jmz) - 1) / j
        
        ##errors.append(abs(x - y))
        
    
        #points.append((x, y))
    
    ##pyt.scatter(unzip(points)[0], unzip(points)[1])
    ##pyt.show()
    ##print "Done."

def getAllMW(mzSpectrum,
             speciesCount = None,
             removalArea = None,
             minimumPeaks = 0,
             peakIterations = 2,
             smoothWidth = 0,
             massRange = None):
    transferWidth = 3
    chargeSequences = assignAllMasses(mzSpectrum, 
                                      removalArea = removalArea,
                                      requiredLength = minimumPeaks,
                                      peakIterations = peakIterations,
                                      massRange = massRange)
    
    mwLabels = [] # (mass-intensity, label) pairs used to label the middle graph.
    #mzSpectrumPeaks = # Peak sequences on the mz spectrum, to 
    mzLabelSequences = [] # [(mass-intensity, label)] sequences to label the bottom graph.
    zcRawSpectrum = set() # Aggregator for the zero-charge spectrum.
    
    alphabet = [chr(65 + x) for x in range(0, 26)]    
    
    masses = []
    for (sequence, mass), letter in zip(chargeSequences, alphabet):
        masses.append((mass, letter))
        peakLabels = []
        for peak, charge in sequence:
            #massPeak = (peak[0] * charge) - (charge * proton), peak[1]
            #zcSpectrum.append(massPeak)
            startMZ, endMZ = peak[0] - transferWidth, peak[0] + transferWidth
            for mz, intensity in [x for x in mzSpectrum if startMZ < x[0] < endMZ]:
                mass = (mz * charge) - (charge * proton)
                zcRawSpectrum.add((mass, intensity))
                
            peakLabels.append((peak, letter))
        mzLabelSequences.append(peakLabels)
    
    zcRawSpectrum = (x for x in sorted(list(zcRawSpectrum), key = lambda x: x[0]))
    
    zcSpectrum = []
    try:
        nextEl = next(zcRawSpectrum)
    except StopIteration:
        raise EmptySpectrum("No valid features were found in the data.")
    
    baseline = True
    done = False
    for mass in range(10001, 100000, 1):
        points = []
        while nextEl[0] < mass:
            points.append(nextEl)
            try:
                nextEl = next(zcRawSpectrum)
            except StopIteration:
                points.append(nextEl)
                done = True
                break
        if done: break
        
        if points:
            if baseline:
                basePoint = mass - 1, 0
                zcSpectrum.append(basePoint)
                baseline = False
            avgPoint = mass, sum([x[1] for x in points]) / len(points)
            zcSpectrum.append(avgPoint)
        elif not baseline:
            basePoint = mass, 0
            zcSpectrum.append(basePoint)
            baseline = True
    zcSpectrum.append((mass, 0))
    
    mwLabels = []
    for mass, _ in masses:
        try:
            #massPoint = max([x for x in zcSpectrum if abs(x[0] - mass) < 10], 
                            #key = lambda x: x[1])
            massPoint = mass, min(zcSpectrum, key = lambda x: abs(x[0] - mass))[1]
        except ValueError: continue # Why does this happen?
        mwLabels.append(massPoint)
    
    mwLabels = list(zip(mwLabels, alphabet))
    if speciesCount:
        speciesCount = int(speciesCount)
        mwLabels = mwLabels[:speciesCount]
        mzLabelSequences = mzLabelSequences[:speciesCount]
        
    subMassSpectra = []
    subSpec = []
    for i in range(0, len(zcSpectrum)-1):
        before, after = zcSpectrum[i], zcSpectrum[i+1]
        subSpec.append(before)
        if after[0] - before[0] > 100:
            subMassSpectra.append(subSpec)
            subSpec = []
    if subSpec:
        subMassSpectra.append(subSpec)
    
    zcSpectrum2 = []
    for subSpec in subMassSpectra:
        if smoothWidth:
            if smoothWidth % 2 == 0:
                smoothWidth -= 1
            mzs, intensities = unzip(subSpec)
            subSpec = list(zip(mzs, savitzkyGolay(intensities, smoothWidth, 3)))
            subSpec = subSpec + [(zcSpectrum[-1][0]+1, 0)]
        zcSpectrum2 += subSpec
    
    
    return mwLabels, mzLabelSequences, zcSpectrum




def mannAlgorithm(mzSpectrum, massRange, chargeRange):
    massSpectrum = defaultdict(int)
    for point in mzSpectrum:
        if point[1] < 3: continue
        for charge in range(chargeRange[0], chargeRange[1]):
            #massSpectrum.append(((point[0]*charge) - (charge*proton), point[1]))
            mass = (point[0]*charge) - (charge*proton)
            if mass >= massRange[1]: break
            elif mass <= massRange[0]: continue
            else: massSpectrum[round(mass)] += point[1]
    
    #print "Stage 2 %s" % time.clock()
    massSpec = []
    for massPt in sorted(massSpectrum.keys()):
        massSpec.append((massPt, massSpectrum[massPt]))
        
    #print "Stage 3 %s" % time.clock()
    #reducedSpec = []
    #step = 10
    #for i in range(step, len(massSpec)-step, step):
        #avgMZ = 0
        #avgInt = 0
        #for j in range(i - step, i + step):
            #avgMZ += massSpec[j][0]
            #avgInt += massSpec[j][1]
        #avgMZ = avgMZ / step * 2
        #avgInt = avgInt / step * 2
        #reducedSpec.append((avgMZ, avgInt))

    smoothedMassSpec = list(zip(unzip(massSpec)[0],
                           savitzkyGolay(unzip(massSpec)[1], 75, 5)))

    
    return smoothedMassSpec
    













## trainingPoints is a sequence of feature sequences; trainingOuput is a sequence of
## dependent features;
## dataPoints is a sequence of sequences same as trainingPoints.
## If trainingPoints == dataPoints, trainingOutput ~~ return value, if it works correctly.
#def emulateFunctionBasedOnExamples(trainingPoints, trainingOutput, dataPoints,
                                   #normFactor = 1, paramLength = 16):
    #from scipy.optimize import leastsq

    #def normalize(shift, scale, data):
        #return (data - shift) / scale    
    
    #def model(param, features):
        #output = 0
        #for i in range(0, len(param)/len(features)):
            #degree = i
            #for c, f in enumerate(features):
                #output += param[(i*len(features))+c]*(f**i)
        #return output
    
    #def costFunc(param, features, target, normFactor):
        #features = unzip([list(f) for f in features])
        
        #iteratee = []
    
        #for i in range(0, len(target)):
            #iteratee.append((target[i], features[i]))
            
        
        #cost = []
        #for idealOutput, featurePoint in iteratee:
            #cost.append(model(param, featurePoint) - idealOutput)
        
        #for p in param:
            #cost.append(p * normFactor)
        
        #return cost    
    
    
    #features = len(trainingPoints)
    #featureScales = []
    #featureShifts = []
    #TP = []
    #for feature in trainingPoints:
        #shift = np.average(feature)
        #scale = np.std(feature)
        #TP.append([normalize(shift, scale, f) for f in feature])
        #featureScales += [scale]
        #featureShifts += [shift]
    ## Scaling output isn't necessary, is it?
    ## Don't be confused by not having done this!
    
    ## Replace explicit list set with variable length list of ones?
    ##initialParameters = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    #initialParameters = [1,1,1,1,1,1]
    #initialParameters = [1] * paramLength
    #parameters, _ = leastsq(costFunc, initialParameters, args = (TP, trainingOutput, normFactor))
    
    #print parameters
    
    #DP = []
    #for feature, scale, shift in zip(dataPoints, featureScales, featureShifts):
        #DP.append([normalize(shift, scale, f) for f in feature])
    
    #DOutput = [model(parameters, dat) for dat in zip(*DP)]

    
    #return DOutput









# Takes in datafile path and produces mz/intensity pairs;
# sum (average?) over the most intense RT region in the datafile.
def findESIPeaks(datafile, timeWindow = None):
    if datafile.lower().endswith('.txt'):
        return readPeakList(datafile)
    
    data = mzFile(datafile)
    
    scanStart, scanStop = data.time_range()
    
    if not timeWindow:
        xic = data.xic(scanStart, scanStop, 100, 2200)
        
        highestPoint = max(xic, key = lambda x: x[1])
        halfMaximum = highestPoint[1] / 2
        
        elution = []
        level = 0
        for i in range(len(xic)-1, 0, -1):
            level = (level/2) + xic[i][1]
            if level >= halfMaximum:
                elution.append(xic[i])
            elif elution and highestPoint[0] > xic[i][0]:
                break
            else:
                elution = []
    else:
        startWindow, endWindow = timeWindow
        assert startWindow >= round(scanStart) and endWindow <= round(scanStop), "Invalid time window."
        
        headers = [x for x in data.headers() if x[3] == 'MS1']
        elution = []
        for header in headers:
            if startWindow > header[0]: continue
            elif endWindow < header[0]: break
            else: elution.append((header[0], None))

    sampleScan = data.scan(1)
    scanLength = len(sampleScan)
    spectrum = [list(x) for x in zip(unzip(sampleScan)[0], [0] * scanLength)]
    for time, _ in elution:
        scan = data.scan(time)
        for i in range(0, scanLength):
            #assert spectrum[i][0] == scan[i][0]
            spectrum[i][1] += scan[i][1]

  
    highestSpectralPoint = max(spectrum, key = lambda x: x[1])
    scaleFactor = highestSpectralPoint[1] / 100.0
    for i in range(0, scanLength):
        spectrum[i][1] = spectrum[i][1] / scaleFactor
    
    return spectrum



def writePeakList(spectrum, outputName):
    peakList = ["m/z\tIntensity\n\n"]
    for peak in spectrum:
        newLine = "{0:.2f}\t{1:.2f}\n".format(peak[0], peak[1])
        peakList += newLine    
    
    output = open(outputName, "w")
    for line in peakList:
        output.write(line)
    output.close()
    
def readPeakList(inputName):
    infile = open(inputName, "r")
    spectrum = []
    for line in infile:
        try:
            point = line.split()
            mz = float(point[0])
            intensity = float(point[1])
        except:
            continue
        spectrum.append((mz, intensity))
        
    infile.close()
        
    return spectrum

def compilePeakLists(inputFiles):
    subSpectra = [readPeakList(x) for x in inputFiles]
    
    spectrum = [list(x) for x in subSpectra[0]]
    for subspectrum in subSpectra[1:]:
        for i in range(0, len(subspectrum)):
            assert spectrum[i][0] == subspectrum[i][0]
            spectrum[i][1] += subspectrum[i][1]
    
    return spectrum

# From SciPy recipe.
def savitzkyGolay(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial

    y = np.array(y)

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')



def smoothSpectrum(spectrum, window = None, order = 3):
    mzs, intensities = unzip(spectrum)
    
    if not window:
        timediffs = [abs(mzs[i] - mzs[i+1]) for i in range(0, len(mzs)-1)]
        timePerPoint = np.average(timediffs)
        window = 6 / timePerPoint
        if window % 2 == 0: window += 1

    smoothIntensities = savitzkyGolay(intensities, window, order)
    return list(zip(mzs, smoothIntensities))





def processZCSpectrum(spectra):
    intensityAverages = []
    intensities = []
    for massPoints, _, _ in spectra:
        intensityAverages.append(np.average(unzip(massPoints)[1]))
        intensities += unzip(massPoints)[1]
    overallAverage = sum(intensities)/len(intensities)
        
    relativeSpectra = defaultdict(list)
    for (massPoints, massPeak, charge), averageInt in zip(spectra, intensityAverages):
        peakRelativePoints = []
        for mass, intensity in massPoints:
            relativeMass = mass #- massPeak[0]
            if intensity:
                scaledIntensity = intensity * (overallAverage/averageInt)
            else:
                scaledIntensity = 0
            assert not isnan(relativeMass)
            assert not isnan(scaledIntensity)
            peakRelativePoints.append((relativeMass, scaledIntensity))
        peakRelativePoints.sort(key = lambda x: x[0])
        relativeSpectra[charge] += peakRelativePoints
    
    zcSpectrum = sorted(sum(list(relativeSpectra.values()), []), key = lambda x: x[0])
    
    def extrapolateZeroes(spectrum):
        filledSpectrum = []
        recentMZ = None
        gap = 10
        for mz, intensity in spectrum:
            if not recentMZ:
                filledSpectrum.append(((mz - gap), 0))
            elif (recentMZ + gap) < mz:
                for newMZ in range(int(np.ceil(recentMZ)), int(np.floor(mz)), gap):
                    filledSpectrum.append((newMZ, 0))
        
            filledSpectrum.append((mz, intensity))
            recentMZ = mz
        filledSpectrum.append((mz + gap, 0))
        return filledSpectrum
    
    assert not any([isnan(x[1]) for x in zcSpectrum]) 
    zeroed = extrapolateZeroes(zcSpectrum)

    assert not any([isnan(x[1]) for x in zeroed])     
    return smoothSpectrum(zeroed, 51, 3)



#END CHARGE TRANSFORMATION ALGORITHMS

#BEGIN GRAPHING/INTERFACE




def placeLabels(dataPoints, highestPoint, labelColor, axis,
                figure = None):

    if not figure:
        ystart, yend = pyt.ylim()
        xstart, xend = pyt.xlim()
    else:
        #ax = figure.get_axes()
        ystart, yend = axis.get_ylim()
        xstart, xend = axis.get_xlim()
    letterHeight = (yend - ystart) / 20
    letterWidth = (xend - xstart) / 60
    figTop = highestPoint
    
    def deriveLimits(x, y, label):
        textHeight = letterHeight
        textWidth = len(label) * letterWidth
        
        top = y + textHeight
        bot = y
        left = x - (textWidth/2.0 + letterWidth)
        right = x + (textWidth/2.0 + letterWidth)
        
        return top, bot, left, right



    # If there is no collision, returns None.  Otherwise,
    # returns adjustment to x and y for each that result in no collision.                     
    def deCollide(first, second):
        top1, bot1, left1, right1 = deriveLimits(*first)
        top2, bot2, left2, right2 = deriveLimits(*second)
        
        if left2 < right1 < right2 or left1 < right2 < right1:
            firstIsLeftmost = left1 < left2
            if firstIsLeftmost:
                deflection = abs(right1 - left2)/1.8
                #return deflection * -1, deflection
            else:
                deflection = abs(right2 - left1)/1.8
                #return deflection, deflection * -1
            return deflection
    
    labels = []
       
    for x, y, label in dataPoints:
        labels.append((x, figTop, label))    
    
    iterations = 0
    collisions = True
    while collisions:
        if iterations > 100:
            break
        else:
            iterations += 1
        collisions = False
        
        for i in range(0, len(labels)):
            for j in range(0, len(labels)):
                if i == j: continue
                x1, y1, label1 = labels[i]
                x2, y2, label2 = labels[j]
                collision = deCollide((x1, y1, label1), (x2, y2, label2))
                
                if collision:
                    #corr1, corr2 = collision
                    if dataPoints[i][0] < dataPoints[j][0]:
                        corr1 = collision * -1
                        corr2 = collision
                    else:
                        corr1 = collision
                        corr2 = collision * -1
                        
                    if xstart < x1 + corr1 < xend:
                        labels[i] = x1 + corr1, y1, label1
                    if xstart < x2 + corr2 < xend:
                        labels[j] = x2 + corr2, y2, label2
                    collisions = True
            
    
    arrowProperties=dict(arrowstyle="-", #linestyle="dashed",
                        color='red',
                        shrinkA=2, #shrinkB=5,
                        #patchA=None,
                        #patchB=None,
                        connectionstyle="angle,angleA=-90,angleB=180,rad=5",
                        )    
    
    for (dpX, dpY, _), (x, y, label) in zip(dataPoints, labels):
        axis.annotate(label, xy = (dpX, dpY), xytext = (x, y),
                      arrowprops= arrowProperties, #{'facecolor':'red', 'width' : 1, 'connectionstyle' : },
                      horizontalalignment = 'center',
                      verticalalignment = 'bottom')
        
    return axis
        
        
        


# spectrum is a list of separator/intensity pairs.
# peakPoints is a list of (separator/intensity)/label pairs;
# each point is printed with label (if any) on top of intensity (if chosen.)
# Separation axis is X, as per usual.
# margin is the fraction of each dimension length that should be a left between furthest out features
# and the edge of the graph.
def plotSpectrum(spectrum, peakPoints, wholeFigure, place, places, xRange = None, yRange = None, 
                 printIntensity = True, yMargin = 5, xMargin = 5, xUnits = "Mass", 
                 yUnits = None, size = (20, 7), labelColor = 'black', smallTicks = False,
                 recoverPeaks = False, carefulLabelling = False, labelParams = None,
                 caption = ''):
    #for peak, label in peakPoints:
        #assert peak in spectrum
    
   
    
    if not peakPoints:
        leastSpectrum = min(spectrum, key = lambda x: x[0])
        mostSpectrum = max(spectrum, key = lambda x: x[0])         
        xRange = (leastSpectrum[0], mostSpectrum[0])
        yRange = (0, max(spectrum, key = lambda x: x[1])[1] + yMargin)
    else:
 
        if not xRange:
            if len(peakPoints) > 1:
                leastPeak = min(peakPoints, key = lambda x: x[0][0])[0]
                mostPeak = max(peakPoints, key = lambda x: x[0][0])[0]
                xMargin = (mostPeak[0] - leastPeak[0]) / xMargin
            else:
                leastPeak = (1, 1)
                mostPeak = (10, 1)
        
            xRange = leastPeak[0] - xMargin, mostPeak[0] + xMargin
            if xRange[1] - xRange[0] < 100:
                xRange = leastPeak[0] - 50, mostPeak[0] + 50
                
        
            
        if not yRange:
            mostPeak = max(spectrum, key = lambda x: x[1])
            yMargin = mostPeak[1] / yMargin
            
            yRange = 0, mostPeak[1] + yMargin
            

    
    spectrum.sort(key = lambda x: x[0], reverse = True)    


    ax = wholeFigure.add_subplot(places, 1, place)
    ax.set_xlim(xRange)
    ax.set_ylim(yRange)
    ax.plot(unzip(spectrum)[0], unzip(spectrum)[1], color = 'black')
    
    labels = []
    for knownPeak, label in peakPoints:
        if recoverPeaks:
            peak = max(min(spectrum, key = lambda x: abs(np.ceil(knownPeak[0]) - x[0])),
                       min(spectrum, key = lambda x: abs(np.floor(knownPeak[0]) - x[0])),
                       key = lambda x: x[1])
        else:
            peak = knownPeak
        
        
        if printIntensity and label:
            labelString = "{0}\n{1:.2f}".format(str(label), knownPeak[0])
        elif printIntensity and not label:
            labelString = "{0:.2f}".format(knownPeak[0])
        elif label and not printIntensity:
            labelString = str(label)
        else:
            labelString = ""
            
        labels.append((peak[0], peak[1], labelString))

    if carefulLabelling:
        #if labelParams:
            #placeLabels(labels, yRange[1] - yMargin, labelColor, ax, labelParams[0], labelParams[1])
        #else:
            #placeLabels(labels, yRange[1] - yMargin, labelColor, ax)                        
        placeLabels(labels, yRange[1] - yMargin, labelColor, ax, wholeFigure)                        
    else:
        for x, y, label in labels:
            ax.text(x, y, label, color = labelColor,
                    horizontalalignment = 'center',
                    verticalalignment = 'bottom')
            

    ax.set_xlabel(xUnits)
    if yUnits: ax.set_ylabel(yUnits)

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='y', left = False, right = False, labelleft = False)
    ax.tick_params(axis = 'y', which = 'minor', bottom = True, direction = 'out', length = 100)
    ax.tick_params(axis='x', top = False)
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    ax.set_title(caption, {'verticalalignment':'bottom'}, loc = 'right', color = 'r')
        
    return wholeFigure
    



            


consolidateMode = False
App = None
Pres = None
ppFilename = None
def slideWritingNonsense(outputfile, imagefile, inputfile):
    global App
    global Pres
    global ppFilename    
    global consolidateMode    
    
    if not MSO:
        errorStr = ("Microsoft Office Object Library not found.  To generate\n"
                    "this library, run the client/makepy.py script located \n"
                    "in the win32com directory of your Python packages, and"
                    "select the version appropriate to your copy of Office.")
        wx.MessageBox(errorStr, "Could not write to Powerpoint.")
        return

    if not consolidateMode:
        App = win32com.client.Dispatch("PowerPoint.Application")
        Pres = App.Presentations.Add()
        Slide = Pres.Slides.Add(1, 12)
        Slide.Shapes.AddPicture(FileName = imagefile, LinkToFile = False, 
                                SaveWithDocument = True, Left = 0, Top = 0,
                                Width = int(1024*0.66), Height = int(768*.66))
        text = Slide.Shapes.AddTextbox(1, # From MSO.constants.msoTextOrientationHorizontal,
                                       Left=10, Top=10, Width = 750, Height = 40)
        text.TextFrame.TextRange.Text = inputfile
        text.TextFrame.TextRange.Font.Size = 10
        Pres.SaveAs(outputfile)
        Pres.Close()
        App.Quit()
        Slide = None
        Pres = None
        App = None
    else:
        if not App:
            global App
            global Pres
            App = win32com.client.Dispatch("PowerPoint.Application")
            Pres = App.Presentations.Add()
            ppFilename = outputfile
        Slide = Pres.Slides.Add(1, 12)
        Slide.Shapes.AddPicture(FileName = imagefile, LinkToFile = False, 
                                SaveWithDocument = True, Left = 0, Top = 0,
                                Width = int(1024*0.66), Height = int(768*.66))
        text = Slide.Shapes.AddTextbox(MSO.constants.msoTextOrientationHorizontal,
                                       Left=10, Top=10, Width = 750, Height = 40)
        text.TextFrame.TextRange.Text = inputfile
        text.TextFrame.TextRange.Font.Size = 10        
        
def powerPointCleanup():
    global App
    global Pres    
    global ppFilename
    print(ppFilename)
    if App: print((5))
    if Pres: print((6))
    assert App and Pres and ppFilename
    Pres.SaveAs(ppFilename)
    Pres.Close()
    App.Quit()
    
    
    
        
        
            


#def getMaxima(spectrum, minElevation = 2):
    #spectrum.sort(key = lambda x: x[0])
    #interpolationPoints = unzip(spectrum)[0]
    
    #spectralCurve = emulateFunctionBasedOnExamples([unzip(spectrum)[0]], unzip(spectrum)[1], 
                                                      #[interpolationPoints], normFactor = 1,
                                                      #paramLength = 5)
    #curveSpectrum = zip(interpolationPoints, spectralCurve)
    
    #if curveSpectrum[0][1] > spectrum[0][1]:
        #peakRegion = False
    #else:
        #peakRegion = True
    #peaks = []
    #previousIntersection = 0
    #for curvePeak, truePeak in zip(curveSpectrum, spectrum):
        #assert curvePeak[0] == truePeak[0]
        #if peakRegion and curvePeak[1] > truePeak[1]:
            #peakRegion = not peakRegion
            #if abs(previousIntersection - curvePeak[0]) > 0.5:
                #peaks.append((previousIntersection, curvePeak[0]))
            #previousIntersection = curvePeak[0]
        #elif not peakRegion and curvePeak[1] < truePeak[1]:
            #peakRegion = not peakRegion
            #previousIntersection = curvePeak[0]
            
    #maxima = []
    #for start, end in peaks:
        #try:
            #trueMaximum = max([x for x in spectrum if start < x[0] < end], key = lambda x: x[1])
            #curveMaximum = max([x for x in curveSpectrum if start < x[0] < end], key = lambda x: x[1])
        #except ValueError:
            #continue
        
        #if trueMaximum[1] > curveMaximum[1] + minElevation:
            #maxima.append(trueMaximum)
        
            
    #return maxima
            

    
def zoomOnMostIntense(spectrum, maxLabels, label):
    assert spectrum, "Empty mass spectrum!"
    
    spectrum.sort(key = lambda x: x[0])
    labelPoint = label[0]
    highPoint = max(spectrum, key = lambda x: x[1])
    
    accumulator = []
    for point in spectrum:
        if point[1]:
            accumulator.append(point)
        else:
            try:
                if min([x[0] for x in accumulator]) < labelPoint[0] < max([x[0] for x in accumulator]):
                    break
                else:
                    accumulator = []
            except ValueError:
                accumulator = []
    #assert highPoint in accumulator
    
    peakRegions = [x for x in iterativeFitBreak(accumulator, curve, 0) 
                   if len(x) > 20 and max([y[1] for y in x]) * 5 > highPoint[1]]
        
    
    if peakRegions:
        peaks = [max(region, key = lambda x: x[1]) for region in peakRegions]
    else:
        # In case iterativeFitBreak failed to find anything.
        # (From the above assertion it is known accumulator is non-empty.)
        peaks = [max(accumulator, key = lambda x: x[1])]
        
    peaks.sort(key = lambda x: x[1], reverse = True)
    
    
    
    if maxLabels == None: 
        maxLabels = 99
    return ((accumulator[0][0], accumulator[-1][0]), 
            list(zip(peaks, [None]*maxLabels)))
    
    
    
    
    

#def zoomOnMostIntense(spectrum, peaks, maxLabels):
    #mostIntense = max(peaks, key = lambda x: x[1])
    #zoomPoint = mostIntense[0]
    
    #threshold = np.average([i for m, i in spectrum])
    #nearSpectrum = []
    #accumulator = 0
    #spectrum.sort(key = lambda x: x[0])
    #accs = []
    #for point in spectrum:
        #accumulator = (accumulator/2) + point[1]
        #if accumulator > threshold:
            #nearSpectrum.append(point)
        #elif point[0] > zoomPoint:
            #break
        #else:
            #nearSpectrum = []
        #accs.append(accumulator)
            
    #maxima = getMaxima(nearSpectrum)
    #zoomBegin = nearSpectrum[0][0]
    #zoomEnd = nearSpectrum[-1][0]
    
    #if len(maxima) > maxLabels:
        #maxima.sort(key = lambda x: x[1], reverse = True)
        #maxima = maxima[:maxLabels]
        
    
    #return (zoomBegin, zoomEnd), zip(maxima, [None] * len(maxima))


#def chargeDecWholeRunMode(dataFile, outputfile, snThreshold,
                          #mzTolerance, maxCharge, mzRange, segmentTime,
                          #segmentStart, consolidate, speciesCount):
def chargeDecWholeRunMode(dataFile, outputfile, mzRange, massRange,
                          speciesCount, removalArea, requiredLength,
                          peakIterations, segmentTime, segmentStart,
                          consolidate, zcLabelCount, zoomLabelCount,
                          recalibrationSettings, smoothWindow):
    data = mzFile(dataFile)
    outputBase = '.'.join(outputfile.split(".")[:-1])
    outputExtension = outputfile.split(".")[-1]
    
    scanStart, scanStop = data.time_range()
    scanStart = max(scanStart, segmentStart)
    
    allScanTimes = [x[0] for x in data.scan_info(scanStart, scanStop) if x[3] == "MS1"]

    for time in numberSequence(scanStart, scanStop, segmentTime):
        start = time
        end = min([time + segmentTime, scanStop])
        segmentScans = [data.scan(x) for x in allScanTimes if start <= x <= end]
        if not segmentScans: continue
        
        sumSpectrum = [[x, y] for (x, y) in segmentScans[0]]
        for spectrum in segmentScans[1:]:
            for i in range(0, len(spectrum)):
                assert spectrum[i][0] == sumSpectrum[i][0]
                sumSpectrum[i][1] += spectrum[i][1]
        
        highestSpectralPoint = max(sumSpectrum, key = lambda x: x[1])
        scaleFactor = highestSpectralPoint[1] / 100.0
        for i in range(0, len(sumSpectrum)):
            sumSpectrum[i][1] = sumSpectrum[i][1] / scaleFactor
            
        if consolidate:
            segmentOutput = outputBase + "." + outputExtension
        else:
            segmentOutput = (outputBase + ".{0:.2f}-{1:.2f}.".format(start, end) 
                             + outputExtension)
        
        try:
            chargeDeconvolution(sumSpectrum, segmentOutput, mzRange, massRange,
                                speciesCount, removalArea, requiredLength,
                                peakIterations, zcLabelCount, zoomLabelCount, 
                                recalibrationSettings, smoothWindow)
        
            print(("Completed segment {0:.2f}-{1:.0f} of {2}".format(start, end, outputfile)))
        except AssertionError as err:
            print(err)
            print(("Failed on segment {0:.2f}-{1:.2f} of {2}".format(start, end, outputfile)))
    
    
        
        
        
    
    

#def chargeDeconvolution(mzSpectrum, outputfile, snThreshold,
                        #mzTolerance, maxCharge, mzRange, speciesCount,
                        #zcLabelCount, zoomLabelCount):
def chargeDeconvolution(mzSpectrum, outputfile, 
                        mzRange = (300, 2000), 
                        massRange = (10000, 100000),
                        speciesCount = 5,
                        removalArea = 3,
                        minimumPeaks = 8,
                        peakIterations = 2,
                        zcLabelCount = None,
                        zoomLabelCount = None,
                        slopeIntercept = (0.0, 0.0),
                        smoothWindow = False,
                        consolidateSet = False,
                        returnFigure = False):                        
    
    rSlope, rIntercept = slopeIntercept
    if consolidateSet:
        global consolidateMode
        consolidateMode = True
    
    #mzSpectrum = findESIPeaks(datafile)
    if outputfile.split(".")[-1] == "txt":
        writePeakList(mzSpectrum, outputfile)
        return
    
    # Filter by mzRange.
    mzSpectrum = [x for x in mzSpectrum if mzRange[0] < x[0] < mzRange[1]]
    # Apply correction factors.
    mzSpectrum = [(mz + (mz * rSlope) + rIntercept, intensity) for (mz, intensity) in mzSpectrum]
    # Rescale intensity values.
    highIntensity = max([x[1] for x in mzSpectrum]) / 100
    mzSpectrum = [(mz, intensity/highIntensity) for (mz, intensity) in mzSpectrum]
   
   
    try:
        mwLabels, labelledMZPeaks, zcSpectrum = getAllMW(mzSpectrum, speciesCount,
                                                         removalArea, minimumPeaks, 
                                                         peakIterations, smoothWindow,
                                                         massRange)
    
        closeupRange, spectralMaxima = zoomOnMostIntense(zcSpectrum, zoomLabelCount, mwLabels[0])
    except EmptySpectrum:
        mwLabels, labelledMZPeaks, zcSpectrum = [], [], [(0,0)]
        closeupRange, spectralMaxima = None, []
    
    triFigure = pyt.figure()
    triFigure = plotSpectrum(mzSpectrum, sum(labelledMZPeaks, []), triFigure, 3, 3,
                             xRange = mzRange, printIntensity = False, xUnits = "M/Z", 
                             labelColor = 'red', smallTicks = True)
    
    triFigure = plotSpectrum(zcSpectrum, mwLabels, triFigure, 2, 3,
                             xRange = massRange, recoverPeaks = True,
                             carefulLabelling = True)

    triFigure = plotSpectrum(zcSpectrum, spectralMaxima, triFigure, 1, 3,
                             xRange = closeupRange, printIntensity = True,
                             xMargin = 0.5, carefulLabelling = True)
    
    triFigure.subplots_adjust(hspace = 0.5)

    if returnFigure:
        return triFigure
    else:
        if outputfile.split(".")[-1] == "ppt" or outputfile.split(".")[-1] == "pptx":
            inputfile = ".".join(outputfile.split(".")[:-1])
            triFigure.savefig("tempDeconvolutionImage.png")
            slideWritingNonsense(outputfile, os.path.join(os.getcwd(), "tempDeconvolutionImage.png"), inputfile)
            os.remove("tempDeconvolutionImage.png")
        else:
            triFigure.savefig(outputfile, bbox_inches='tight')
    
        triFigure.clear()    
        print("Done!")
    






        



            
        

def chargeTransform(inputFile,
                    outputFile = None, 
                    mzRangeStart = 300,
                    mzRangeEnd = 2000,
                    massRangeStart = 10000,
                    massRangeEnd = 100000,
                    speciesCount = 5,
                    removalArea = 2,
                    minimumPeaks = 8,
                    peakIterations = 2,
                    zcLabelCount = 10,
                    zoomLabelCount = 5,
                    segmentMode = False,
                    startSegment = 0,
                    segmentSize = 5,
                    consolidate = False,
                    recalSlope = 0,
                    recalIntercept = 0,
                    smoothWidth = 0
                    ):
    """
    Performs charge state transformation on an MS data file.  If outputFile is
    specified, the output is saved as an image of the given type; otherwise the
    function returns a matplotlib figure instance containing the transformation
    rendering.
    
    
    mzRangeStart, mzRangeEnd <- Beginning and ending of the m/z range that will
    be considered when looking for ion peaks.    
    
    massRangeStart, massRangeEnd <- Beginning and ending of the possible analyte
    mass range.  Molecules outside of this range will not be found, or may produce
    spurious results.

    speciesCount <- Maximum number of chemical species that will be 
    investigated in a single run.
    
    removalArea <- The M/Z region around each peak that is removed from
    consideration when the peak is accounted for.
    
    minimumPeaks <- Minimum peaks corresponding to a given mass that must be
    found for that mass to be considered a valid hit.
    
    peakIterations <- Iterations of the peak-finding algorithm performed before
    analysis.  2 is typically a good number.
    
    zcLabelCount <- Maximum number of chemical species that will be labelled
    on the zero-charge spectrum results.
    
    segmentMode <- If False, chargeTransform takes a spectrum averaged over
    the region of greatest intensity in the time dimension.  If True,
    chargeTransform performs multiple analyses, each over a different sequential
    segment of scan time.
    
    segmentSize <- Length in minutes of a analysis segment (only used if 
    segmentMode == True.)
    
    startSegment <- The segment at which chargeTransform will start performing
    analyses (only used if segmentMode == True.)
    
    consolidate <- If True, all output files will be put in the same output
    (only used if output type is .ppt and segmentMode == True.)
    
    recalSlope <- Slope of the correction line applied to scan mzs; 0 applies no
    mz-varying correction.
    
    recalIntercept <- Y-intercept point of the correction line applied to the scan
    mzs; 0 applies no constant correction.
    
    smoothWidth <- The width of the window used in the Savitzky-Golay smoothing
    process on the zero-charge spectrum.  Enter 0 or None for no smoothing.
    """
    
    outputFile = os.path.abspath(outputFile)
    
    assert not (segmentMode and (not outputFile)), "Cannot return segment-mode analysis as a Figure object!"
    
    if not segmentMode:
        spectrum = findESIPeaks(inputFile)
        if outputFile:
            return chargeDeconvolution(spectrum, outputFile, (mzRangeStart, mzRangeEnd), (massRangeStart, massRangeEnd),
                                       speciesCount, removalArea, minimumPeaks, peakIterations,
                                       zcLabelCount, zoomLabelCount, (recalSlope, recalIntercept),
                                       smoothWidth)
        else:
            return chargeDeconvolution(spectrum, None, (mzRangeStart, mzRangeEnd), (massRangeStart, massRangeEnd),
                                       speciesCount, removalArea, minimumPeaks, peakIterations,
                                       zcLabelCount, zoomLabelCount, (recalSlope, recalIntercept),
                                       smoothWidth,
                                       returnFigure = True)            
    else:
        return chargeDecWholeRunMode(inputFile, outputFile, (mzRangeStart, mzRangeEnd), (massRangeStart, massRangeEnd),
                                     speciesCount, removalArea, minimumPeaks, peakIterations,
                                     segmentTime, segmentStart, consolidate,
                                     zcLabelCount, zoomLabelCount, (recalSlope, recalIntercept))
    
        
 