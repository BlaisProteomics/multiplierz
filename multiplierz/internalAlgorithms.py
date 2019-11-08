from math import floor, ceil, sqrt
from collections import defaultdict, deque
import gzip

import os, sys

import bisect
from itertools import chain
from functools import reduce

try:
    from numpy import average
except ImportError:
    # Someday, hopefully, numpy will no longer be a multiplierz dependency.        
    def average(xs, weights = None):
        if weights:
            assert len(weights) == len(xs), "Weights must be same length as input sequence."
            return (sum([x*w for x, w in zip(xs, weights)]) 
                    / 
                    float(len(xs)) * reduce((lambda x, y: x * y), weights, 1))
        else:
            return sum(xs) / float(len(xs))



  

def iselect(i, sequence):
    return (x[i] for x in sequence)
def select(i, sequence):
    return [x[i] for x in sequence]



def counts(sequence):
    return [(x, sequence.count(x)) for x in sorted(set(sequence))]

def matching_file(filename, directory = None, tag = None):
    if directory is None:
        directory = os.path.dirname(filename)
    
    filename = filename.lower()
    files = [x.lower() for x in os.listdir(directory)]
    stem = os.path.basename(filename).split('.')[0]
    match_files = [x for x in files if os.path.basename(x).split('.')[0] == stem
                   and x != filename]
    if tag:
        match_files = [x for x in match_files if x.endswith(tag.lower())]
    if len(match_files) > 1:
        raise IOError("Ambigous file matches. %s" % match_files)
    elif not match_files:
        raise IOError("No file match for %s" % filename)
    else:
        return os.path.join(directory, match_files[0])

def average_scan(scans):
    # Probably move to spectral_processing.py when I'm more sure about it.
    
    dists = list(chain(*[[s[i+1][0] - s[i][0] for i in range(len(s)-1)]
                         for s in scans]))
    max_width = min(dists) - 0.000001
    aggregated_scan = aggregate_points(list(chain(*scans)), 
                                       MAX_WIDTH = max_width)
    
    
    avg_scan = []
    for agg_pts in aggregated_scan:
        avg_mz = average([x[0] for x in agg_pts])
        
        avg_int = sum([x[1] for x in agg_pts]) / len(scans)
        # So that points only in a few scans aren't preferentially treated.
        
        avg_scan.append((avg_mz, avg_int))
    
    return sorted(avg_scan)
    
    


def aggregate_points(pointlist, distance_function = None, MAX_WIDTH = 0.025):
    """
    Given a list of points (which may be floats or tuples with the
    position-relevant value in the first position), returns a list of
    grouped points such that all points within a group are within MAX_WIDTH
    of each other.
    """
    peaks = []
    agg = []
    
    if distance_function is not None:
        distance = distance_function
    elif isinstance(pointlist[0], tuple):
        def distance(x, y):
            return x[0] - y[0]
        pointlist.sort()
    else:
        def distance(x, y):
            return x - y
        pointlist.sort()

    def hard_case(group):
        # Should be refactored so that the distance function is only called
        # once per consecutive pair, in case this is being used on objects with
        # expensive distance calculations.
        if len(group) <= 2:
            return [group] # Within limit by construction.
        elif distance(group[-1], group[0]) < MAX_WIDTH:
            return [group]
        else:
            splitpt = max(list(range(len(group)-1)),
                          key = lambda x: distance(group[x+1], group[x])) + 1
            assert group[splitpt:] and group[:splitpt]
            return chain(hard_case(group[splitpt:]), hard_case(group[:splitpt]))
        
 
    for pt in pointlist:
        if (not agg) or distance(pt, agg[-1]) < MAX_WIDTH:
            agg.append(pt)
        else:
            if distance(agg[-1], agg[0]) < MAX_WIDTH:
                peaks.append(agg) # The simple case.
            else:
                peaks += list(hard_case(agg))
            agg = [pt]
    if agg:
        if distance(agg[-1], agg[0]) < MAX_WIDTH:
            peaks.append(agg)
        else:
            peaks += hard_case(agg)
    
    return peaks
        

def peak_in(spectrum, low_mz, high_mz):
    # Returns biggest peak within given range, or None if there are no peaks in range.
    # Uses that bisect follows standard sort-by-first-item conventions.
    low_i = bisect.bisect_left(spectrum, low_mz)
    high_i = bisect.bisect_right(spectrum, hi_mz)
    try:
        return max(spectrum[low_i:high_i], key = lambda x: x[1])
    except ValueError:
        return None

def peak_at(spectrum, mz, tol):
    # Returns biggest peak within given range, or None if there are no peaks in range.
    # Uses that bisect follows standard sort-by-first-item conventions.    
    low_i = bisect.bisect_left(spectrum, (mz-tol, None))
    high_i = bisect.bisect_right(spectrum, (mz+tol, None))
    try:
        return max(spectrum[low_i:high_i], key = lambda x: x[1])
    except ValueError:
        return None

def table_lookup(rows, key):
    return dict((row[key], row) for row in rows)

# Even dill can't pickle generators, though, so this doesn't work.
def async_generator(generator, backlog_size = 10):
    from multiprocess import Process, Queue
    pipe = Queue(maxsize = backlog_size)
    done_flag = 'FO0B@R'
    def gen_reader(pipe):
        for thing in generator:
            pipe.put(thing) # And smoke it.
        pipe.put(done_flag)
    
    reader_proc = Process(target = gen_reader, args = (pipe,))    
    reader_proc.start()
    
    while True:
        next_thing = pipe.get()
        if next_thing == done_flag:
            break
        yield next_thing
        
    reader_proc.join()

# Overhead from task management/external loops can be a slowdown if the
# same thread is handling outside stuff and tasks at the same time; new
# version that spawns a separate process to manage the process pool?
# (Or a separate thread, possibly?)
# Try to have subtasks return their data directly to the main process
# without passing through the manager process; inter-thread communication of
# large data is also a slowdown.
def assign_multiprocess_ext(function, data, pool_args = {}, **task_args):
    from multiprocess import Queue, Process, cpu_count
    from Queue import Full, Empty
    from time import sleep
    process_count = pool_args.get('processes', cpu_count()-1)
    input_pipe, output_pipe, control_pipe = (Queue(process_count), Queue(process_count),
                                             Queue(process_count))
    stop_signal = hash('OK STOP NAO.')
    def multiprocessor(inpipe, outpipe, controlpipe):
        def returner_process(inp, outp, task):
            args, kwargs = inp.get()
            outpipe.put(task(*args, **kwargs))
            return True
            
        jobs = []
        while True:
            done = [x for x in jobs if x.ready()]
            if done:
                jobs = [x for x in jobs if x not in done] # Avoids race condition!       
            else:
                sleep(0.1)
                
            for thing in done:
                thing.successful()
                assert thing.get()
            while len(jobs) < process_count:
                cmd = controlpipe.get()
                if cmd == stop_signal:
                    break
                elif cmd == True:
                    newjob = Process(target = returner_process, 
                                        args = (inpipe, outpipe))
                    newjob.start()
                    jobs.append(newjob)
                                # I *think* the pipes have to be passed explicitly,
                                # but I haven't checked.
                else:
                    raise Exception
        outpipe.put(stop_signal)
    
    multiproc_proc = Process(target = multiprocessor,
                             args = (input_pipe, output_pipe, control_pipe))
    multiproc_proc.start()

    if isinstance(data, list):
        data = (x for x in data)
    nexttask = next(data)
    while True:
        try:
            input_pipe.put_nowait(nexttask)
            control_pipe.put_nowait(True)
            nexttask = next(data)
        except Full:
            pass
        except StopIteration:
            break
        try:
            yield output_pipe.get_nowait()
        except Empty:
            sleep(0.1)

    control_pipe.put(stop_signal)
    while True:
        try:
            out = output_pipe.get()
            if out == stop_signal:
                break
            else:
                yield out
        except Empty:
            sleep(0.1)
    
    multiproc_proc.join()
        
            
        
                
            
                

def assign_multiprocess(function, data, pool_args = {}, **task_args):
    # Pool.map() is convenient, but leads to situations where one last
    # job-batch takes way longer than the others, or each iteration leaves
    # all but one process waiting for the last to finish. Distributing jobs
    # one-by-one is much more efficient in some cases.
    #
    # WARNING- RESULTS ARE NOT RETURNED IN ANY PARTICULAR ORDER.
    # 
    # Because I've been reading about Haskell, "function" can now be a list of
    # functions that will be applied as a single composed function.
    #import multiprocessing
    from time import sleep
    
    if 'processes' in task_args:
        import warnings
        warnings.warn("You probably meant to put 'processes' in pool_args.")
    
    if isinstance(function, list):
        # Allows multiprocessing over non-base-level functions.
        import pathos.multiprocessing as multiprocessing
        import traceback
        
        assert all(map(callable, function)), "[Function/list of functions], [list of data objects]!"
        assert not task_args, "Can't include kwargs to composed functions!"
        function_list = function
        def composed_function(*data_item):
            for func in function_list:
                if isinstance(data_item, str):
                    data_item = (data_item,)
                try:
                    data_item = func(*data_item)
                except Exception as err:
                    print('\n\n\n################')
                    traceback.format_exc()
                    print('################\n\n\n')
                    print(("In: " + str(func) + '\n'))
                    raise err
            return data_item
        function = composed_function
    else:
        import multiprocessing
    from multiprocessing import TimeoutError    
    
    process_count = pool_args.get('processes', multiprocessing.cpu_count()-1)
    
    workforce = multiprocessing.Pool(**pool_args)
    try:
        tasks = data[:]
    except TypeError:
        tasks = data
    jobs = []
    results = []
    try:
        returnquota = len(data)
    except TypeError:
        returnquota = None
    returncount = 0
    while tasks or jobs:
        # Should launch new jobs *before* yielding.
        done = [x for x in jobs if x.ready()]
        if done:
            jobs = [x for x in jobs if x not in done] # Avoids race condition!
        while tasks and len(jobs) < process_count:
            if isinstance(tasks, list):
                newtask = tasks.pop()
            else:
                newtask = next(tasks)
            if isinstance(newtask, str) or len(newtask) == 1:
                newtask = (newtask,)
            jobs.append(workforce.apply_async(function, args = newtask,
                                              kwds = task_args))
            
        for done_thing in done:
            done_thing.successful()
            try:
                yield done_thing.get(10) 
                returncount += 1
            except TimeoutError:
                print("Failure to return in 10 seconds!")
                jobs.append(done_thing)
        sleep(1)
    
    if returnquota:
        assert returncount == returnquota
    
    workforce.close()
    workforce.join()
    
    #return results





def floatrange(start, stop, step = 1):
    cur = start
    while cur < stop:
        yield cur
        cur += step


# For converting profile data or XIC to an "official"
# index-corresponds-to-physical-position intensity vector.
def pts_to_bins(pts, bincount):
    pts.sort()
    startpt, stoppt = pts[0][0], pts[-1][0]
    binwidth = (stoppt - startpt)/float(bincount)
    bins = defaultdict(float)
    prevpt = pts.pop(0)
    nextpt = pts.pop(0)
    for lmz in floatrange(startpt, stoppt, binwidth):
        rmz = lmz + binwidth
        #if rmz > nextpt[0] and lmz < nextpt[0]:
            ## In-between- could be more complicated.
            #bins[lmz] += nextpt[1]
        #else: #rmz < nextpt[0]:
        while pts and rmz > nextpt[0]:
            bins[lmz] += prevpt[1]
            prevpt = nextpt
            nextpt = pts.pop(0)                
        bins[lmz] += prevpt[1]
            
    
    return sorted(bins.items())

import numpy as np
def linear_bin(pts, start, stop, bincount):
    binpts = (x for x in np.linspace(start, stop, bincount+1))
    pts.sort()
    bins = np.array(shape = (bincount, 1))
    binpt = next(binpts)
    for x, y in pts:
        while binpt < x:
            binpt = next(binpts)
        bins[index] += y
    return bins
            
            
def psm_assignment(psm):
    try:
        return psm['Peptide Sequence'], psm['Variable Modifications'], psm['Charge']
    except KeyError:
        return psm['Peptide']
    

def splitOnFirst(string, char):
    """
    Splits on only the first instance of a character.  Like the .split()
    method, removes the splitting character, but other instances of the
    character are retained.  Always returns a 2-tuple, with the whole
    string in the first element if the char is not found.
    """
    index = string.find(char)
    if index >= 0:
        return string[:index], string[index+1:]
    else:
        return string, ''
        
        
def insert_tag(filename, tag):
    words = filename.split('.')
    if words[-1].lower() in ['gz']:
        return '.'.join(words[:-2] + [tag, words[-2], words[-1]])
    else:
        return '.'.join(words[:-1] + [tag, words[-1]])
        

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
    if not os.path.exists(directory):
        raise IOError('%s not found.' % directory)
    if not recursive:
        return [os.path.join(directory, x) for x in
                os.listdir(directory) if x.lower().endswith(ext.lower())
                and x[0] != '~']
    else:
        output = []
        for path, subdirs, files in os.walk(directory):
            for filename in files:
                if filename.lower().endswith(ext.lower()):
                    if 'xls' in ext.lower() and '~' in filename: # Obnoxious Excel thing.
                        continue
                    output.append(os.path.join(path, filename))
        return output



def retryTimes(function, times, exceptions):
    for _ in range(times-1):
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

def PPM_bounds(ppm, mass):
    val = (mass/1000000.0)*ppm
    return mass - val, mass + val
    
    
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

    uncorrected = numbers[plex+2:2*plex+2]
    corrected = numbers[1:plex+1]

    datafile = words[0] # Usually.
    scanNum = words[1] # Again, usually.

    mz = desc.split('-')[1]

    return {'datafile':datafile,
            'scan':scanNum,
            'mz':mz,
            'corrected':corrected,
            'uncorrected':uncorrected}

    
class NaiveProximitySequence(object):
    def __init__(self, sequence, indexer = (lambda x: x)):
        self.sequence = [(indexer(x), x) for x in sequence]
        self.indexer = indexer
    
    def __getitem__(self, thing):
        return min(self.sequence, key = lambda x: abs(thing - x[0]))[1]
    
    def remove(self, thing):
        index = self.indexer(thing)
        del self.sequence[self.sequence.index([x for x in self.sequence if x[0] == index][0])]
        
    def asList(self):
        return list(zip(*self.sequence))[1]
    
    
    
    
    

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
        if isinstance(key, str):
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
    
    threshold = average(list(zip(*scan))[1])
    peaks = []
    peak = []
    for pt in scan:
        if pt[1] > threshold:
            peak.append(pt)
        elif peak:
            centroid = average(list(zip(*peak))[0], weights = list(zip(*peak))[1]), max(list(zip(*peak))[1])
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
                  
# THIS COPY SHOULD BE REMOVED is redundant with copy in spectral_process.py.
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
            raise NotImplementedError("Called scan_features on a profile-mode spectrum.")

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
        for chg, chargeSet in list(activeSets.items()):
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
    # would correctly extend it.
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
            raise NotImplementedError("Called scan_features on a profile-mode spectrum.")

    if enforce_isotopic_ratios:
        global isotopicRatios
    else:
        isotopicRatios = defaultdict(lambda: (-1000000, 1000000))
    lowC13Rat, highC13Rat = isotopicRatios[0, 1]
    
    
    if isinstance(max_charge, list):
        chargeFractions = [(x, 1.0/x) for x in max_charge]
        max_charge = max(max_charge)
    else:
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
                            for oldPt in reversed(finishedSet['envelope']): # Reverse so as to not simply re-create the old envelope.
                                possibleNexts = [(z, oldPt[0] + cF) for z, cF in chargeFractions]# if oldPt[0] + cF > pmz]
                                if possibleNexts:
                                    #re_attached = False
                                    #for actPt in sorted(activePts.keys(), key = lambda x: x[1], reverse = True):
                                        #for z, pn_mz in possibleNexts:
                                            #if (inPPM(tolerance, actPt[0], pn_mz) and 
                                                #lowC13Rat <= (oldPt[1] / actPt[1]) <= highC13Rat):
                                                
                                                #nextMZ = pmz + chargeFracDict[z]
                                                #newActiveSet = {'envelope':[oldPt, actPt],
                                                                #'charge':z,
                                                                #'next':nextMZ,
                                                                #'prevInt':actPt[1]}
                                                #activeSets[z].append(newActiveSet)  
                                                #del activePts[actPt]
                                                #re_attached = True                                                
                                                #break
                                    #if not re_attached:    
                                        activePts[oldPt] = possibleNexts
                                else:
                                    unassigned.append(oldPt)
                        else:
                            unassigned += finishedSet['envelope']                    

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
        for actPt in sorted(list(activePts.keys()), key = lambda x: x[1], reverse = True):
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
    for charge, things in list(activeSets.items()):
        for thing in things:
            if len(thing['envelope']) >= min_peaks:
                envelopes[charge].append(thing['envelope'])
            else:
                unassigned += thing['envelope']            
                
    unassigned += list(activePts.keys())
    


    ## Cleanup phase- each envelope checks for unassigned peaks that
    ## would correctly extend it.
    if cleanup:
        raise NotImplementedError

    if correction != None:
        return dict(envelopes), unassigned, newCorrection
    else:
        return dict(envelopes), unassigned            
    




class ProximityIndexedSequence(object):
    # Turns out there's a python recipe to do just about this but better, but
    # now I have a bunch of legacy code that uses this interface!  Ah well.
    def __init__(self, sequence, indexer = lambda x: x, *args, **kwargs):
        
        self.indexer = indexer
        self.lookup = collectByCriterion(sequence, indexer)
        self.indexes = sorted(self.lookup.keys())
        
    def seal(self):
        pass
    
    # add and remove are linear time, is the only problem with this version.
    # Pretty sure in most use cases this is acceptable.
    def add(self, value):
        index = self.indexer(value)
        bisect.insort_left(self.indexes, index)
        if index not in self.lookup:
            self.lookup[index] = [value]
        else:
            self.lookup[index].append(value)
    def remove(self, value): 
        index = self.indexer(value)
        self.indexes.remove(index)
        self.lookup[index].remove(value)
        
    def __getitem__(self, index):
        try:
            neardex = bisect.bisect_left(self.indexes, index)
            founddex =  min(self.indexes[neardex-1:neardex+1],
                            key = lambda x: abs(index - x))
            return self.lookup[founddex][0]
        except IndexError:
            return self.lookup[self.indexes[-1]][0]
        except ValueError:
            return self.lookup[self.indexes[0]][0]
    def asList(self):
        return list(chain(*(self.lookup[i] for i in self.indexes)))
    def __iter__(self):
        for value in chain(*(self.lookup[i] for i in self.indexes)):
            yield value
    def returnRange(self, begin, stop):
        startdex = bisect.bisect_left(self.indexes, begin)
        stopdex = bisect.bisect_right(self.indexes, stop)
        output = []
        for index in self.indexes[startdex:stopdex]:
            output += self.lookup[index]
        return output
    
    def rebalance(self, *args, **kwargs):
        pass
        
        
ProxSeq = ProximityIndexedSequence # Less of a mouthful.



    
# Just for debugging.
def deisochart(envelopes, unassigned):
    import matplotlib.pyplot as pyt
    from random import seed, shuffle
    
    seed(1)
    
    def lines(thing, charge = None, **etc):
        pyt.vlines(list(zip(*thing))[0], [0]*len(thing), list(zip(*thing))[1], **etc)
        if charge:
            for pt in thing:
                pyt.text(pt[0], pt[1], str(charge))
    
    lines(unassigned, color = 'b')
    length = len(sum(list(envelopes.values()), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in list(envelopes.items()):
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
        pyt.vlines(list(zip(*thing))[0], [0]*len(thing), 
                   [x*upness for x in zip(*thing)[1]], **etc)
        if charge:
            for pt in thing:
                pyt.text(pt[0], upness * pt[1], str(charge))
    
    lines(unassigned1, color = 'b')
    length = len(sum(list(envelopes1.values()), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in list(envelopes1.items()):
        for env in chgEnvs:
            label = "%s(%s)" % (i, chg)
            lines(env, charge = label, color = colors[i], linewidth = 3)
            i += 1
            
            
    lines(unassigned2, upness = -1, color = 'b')
    length = len(sum(list(envelopes2.values()), []))
    colors = pyt.cm.Set1([x/float(length) for x in range(0, length)])
    shuffle(colors)
    i = 0
    for chg, chgEnvs in list(envelopes2.items()):
        for env in chgEnvs:
            label = "%s(%s)" % (i, chg)
            lines(env, charge = label, upness = -1, color = colors[i], linewidth = 3)
            i += 1
    
    
    pyt.show()  
    
    

