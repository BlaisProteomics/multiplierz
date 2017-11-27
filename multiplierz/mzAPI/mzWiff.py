from comtypes.client import CreateObject
import os
from collections import defaultdict
import warnings
from multiplierz.mzAPI import mzScan, mzFile as mzAPImzFile
from multiplierz.internalAlgorithms import centroid as centroid_func, ProximityIndexedSequence


__author__ = 'William Max Alexander'

debug = True


# DEV NOTE: Attempts to reconcile access to the ABSCIEX idea of how MS works
# with the abstract RAW-esque way that it works keep hitting difficulties. It
# looks like the best solution to that would be to have two separate WIFF
# mzFile variants, one of which abstracts over to provide a RAW-like
# interface (nicely compatible, hopefully, interleaving experiments by cycle)
# and the other which demands sample/experiment specifications rigorously
# (making no attempt to push those under the rug, as the current
# implementation does.)
# _explicit_numbering will not keep track of "current" sample and
# experiment numbers at all, and require each to be set explicitly on each
# relevant call.  _implicit_numbering will derive what sample and experiment
# to access at each call.  The underlying COM object does what it can to 
# avoid unnecessary sample/experiment switching, but that will be left out
# of mzFile-level code for the time being.
# _explicit_numbering also has decrement-to-zero-indexed duty.


# Some terminology: "cycle" should mean a given run through the set of
# experiments (though not all experiments are present in every cycle!),
# 'scan' should mean a specific measured spectrum, numbered in order of
# acquisition.


class mzFile_implicit_numbering(mzAPImzFile):
    """
    mzAPI interface class for WIFF files.
    
    When initialized, the 'sample' argument sets what sample will be accessed
    by this instance of the object; all calls will pertain only to data from
    that sample. If not specified, this defaults to sample 1. To switch
    samples, change the .sample attribute appropriately.
    
    Scans in a WIFF file are indexed according to time, and scans from
    experiment 0 are assumed to be MS1-level scans. All other experiments are
    assumed to be MS2-level scans. So, scan number 1 will generally be an MS1
    scan, scans 2-N will be MS2 scans, and so on. .scan_info() will return
    information on which particular scans are on what level throughout the
    file.
    
    In order to access scans by experiment and sample number explicitly, initialize
    mzFile with the 'experiment_numbering' argument set to True.
    """
    
    def __init__(self, datafile, sample = 1, **etc):
        self.data_file = datafile
        self.file_type = 'wiff'
        
        self.data = mzFile_explicit_numbering(datafile,
                                              sample = sample,
                                              experiment = 1)
        
        self.sample = sample
        self.exp_num = self.data.source.GetExperiments(sample-1)

        scans_present = [(c, e) for c, e, s in zip(*self.data.scan_info())[2] if s == self.sample]
        self.make_explicit = dict(enumerate(scans_present))
        self.make_implicit = dict((y, x) for x, y in self.make_explicit.items())
        
    def scan(self, scan, **kwargs):
        if isinstance(scan, float):
            if scan != int(scan):
                raise RuntimeError, ("Scan time specification is ambiguous "
                                     "with implicit experiment numbering; use "
                                     "explicit_numbering mode mzFile.")
            else:
                scan = int(scan)
                
        cycle, experiment = self.make_explicit[scan]
        return self.data.scan(cycle, experiment = experiment, sample = self.sample,
                              **kwargs)
    
    
    def scan_info(self, start_scan = 0, stop_scan = 999999, start_mz = 0, stop_mz = 999999):
        exp_info = self.data.scan_info(start_scan, stop_scan, sample = self.sample)

        return [(rt, precM, self.make_implicit[cycle, experiment], level, centroid)
                for (rt, precM, (cycle, experiment, sample), level, centroid)
                in exp_info
                if sample == self.sample]
    
    def xic(self, start_time = 0, stop_time = 999999, start_mz = 0, stop_mz = 2000,
            filter = None, experiment = 1):
        """
        Gets the eXtracted Ion Chromatogram of the given time- and mz-range.
        
        For a standard MS1 XIC, leave 'experiment' set to 1 (for usual methods.)
        """
    
        return self.data.xic(start_time, stop_time, start_mz, stop_mz, 
                             sample = self.sample, experiment=experiment, 
                             filter = filter)
    
    def time_range(self):
        """
        Gets the total retention-time range for the data.
        """
        
        return self.data.time_range(self.sample)
    
    def scan_range(self):
        # Source gives cycle count.
        start, stop = self.data.scan_range(sample = self.sample)
        return start, len(self.make_implicit)
    
    def scan_for_time(self, rt):
        """
        Gets the scan index for the specified retention time.
        """
        
        return self.make_implicit[self.data.scan_for_time(rt, sample = self.sample)[:2]]
    
    def time_for_scan(self, scan):
        cycle, exp = self.make_explicit[scan]
        return self.data.time_for_scan(cycle, exp, self.sample)
        
    def filters(self):
        return self.data.filters()
    
    def scan_modes(self):
        return self.data.scan_modes()
    
    def headers(self):
        return self.data.scan_info()
    
    def MRM_info(self):
        return self.data.MRM_info()
    
    def MRM_channels(self):
        return self.data.MRM_channels()
    
    def MRM_scan(self, cycle):
        return self.data.MRM_scan(cycle)
    
    def close(self):
        self.data.close()
        

class mzFile_explicit_numbering(mzAPImzFile):
    """
    mzAPI interface class for WIFF files.
    
    Scans in a WIFF file are indexed by sample, experiment, and cycle.  Briefly,
    each sample corresponds to a MS run contained in the file.  Each cycle contains
    an MS1 scan and the MS2 scans performed from that point, each of which has
    a different experiment number.  Note that as a result, all MS1 scans are 
    typically in experiment 1.    
    """
    
    def __init__(self, data_file, sample = 1, experiment = None, **etc):
        self.file_type = 'wiff'
        self.data_file = data_file
        
        #try:
        self.source = CreateObject("{9eabbbb3-5a2a-4f73-aa60-f87b736d3476}")
        #except WindowsError as err:
            #print "WiffReaderCOM.dll not found in registry."
            #raise err        
    
        if not os.path.exists(data_file + '.scan'):
            raise IOError, "%s.scan not found!" % data_file
    
        self.source.OpenWiffFile(os.path.abspath(data_file))
        
        self.sample = sample
            
        self.sample_count = self.source.GetSamples()
        self.experiment_count = self.source.GetExperiments(sample-1)
        
        # Will only store list for default parameters, 
        # since that's what's used by headers() and filters().
        self._scan_info = {} 
        
        # Until we get around to fixing the COM-level time/scan 
        # translation functions, can use the scan_info data to
        # provide a lookup.
        scan_rt = [(scan, rt) for rt, _, scan, _, _ in self.scan_info()]
        self._scan_to_rt = dict(scan_rt)
        self._rt_to_scan = ProximityIndexedSequence(scan_rt, indexer = lambda x: x[1])
        self._rt_to_scan.seal()
        
        
        
       
    
    def scan(self, scan_name, experiment, sample = None, centroid = False):
        """
        Retrieve a spectrum by the scan name.
        
        scan_name may either be the (int) cycle number, a (float) retention
        time value, a double (cycle number, experiment number) or a triple
        (cycle number, experiment number, sample number).
        
        If .set_sample and .set_experiment have not been called for a given mzFile
        instance, calls to .scan must specify these values.  Sample
        and experiment numbers given to .scan are not saved.
        """
        
        if sample == None:
            sample = self.sample            
        
        #if centroid:
            #print 'WARNING- Centroid argument to scan() does not currently work for WIFF files.'
                
        if isinstance(scan_name, int):
            cycle = scan_name  
        elif isinstance(scan_name, float):
            cycle = self.scan_for_time(scan_name, experiment, sample)
        else:
            raise NotImplementedError, "scan_name must be float or int."
     
        scan = zip(*self.source.GetSpectrumIndex(sample-1, experiment-1, cycle-1))
        if centroid:
            scan = centroid_func(scan)
        return scan
    
    
    def scan_info(self, start_cycle = 0, stop_cycle = 999999, experiment = None, sample = None):
        """
        Returns a list of [(time, mz, scan_name, scan_type, scan_mode)] in 
        the time and mz range provided in the sample and experiment specified
        (or all experiments and previously set sample, if not specified.)
        
        If sample is omitted, results are returned for the first sample.  If
        experiment is omitted, results are returned for all experiments
        in the sample.
        """
        
        if not sample:
            sample = self.sample
        
        if sample-1 in self._scan_info:
            return [x for x in self._scan_info[sample-1] if
                    start_cycle <= x[2][0] <= stop_cycle and
                    (experiment == None or x[2][1] == experiment)]
               
        # Zero- versus one-indexing gets dicy here.
        
        if sample:
            samples = [sample-1]
        else:
            raise Exception
        
        if experiment:
            expCounts = {sample-1 : [experiment]}
        else:
            expCounts = dict([(x, range(0, self.source.GetExperiments(x))) for x in samples])
        
        scaninfo = []
        cycleInfo = {}
        # samples is always length 1, due to the above checks, so this
        # is currently more complicated than it needs to be.
        for sample in samples:
            #expInfo = {}
            #for exp in expCounts[sample]:
                #expInfo[exp] = self.source.ExperimentInfo(sample, exp)
                
            ## Obnoxious to have to call this simply for percursor masses.
            #cycleData = self.source.GetSampleData(sample) 
            #for exp, cyc, rt, mass, cole in zip(*cycleData):
                #cycleInfo[sample, int(exp), int(cyc)] = rt, mass, cole                
                
            cycles = self.source.GetNumCycles(sample)
            start_pt = max([0, start_cycle])
            stop_pt = min([cycles, stop_cycle])
            for cycle in range(start_pt, stop_pt):
                for exp in expCounts[sample]:
                    (level, precM, centroid, tof,
                     colE, minInt, maxInt) = self.source.GetScanData(sample, exp, cycle)
                    rt = self.source.GetRTOfScan(sample, exp, cycle)
                    
                    level = 'MS%d' % int(level)
                    precM = precM if precM > 0 else 0
                    centroid = 'p' if centroid == 0 else 'c'
                    
                    if (level != 'MS1') and not precM:
                        if debug:
                            pass
                            #assert not self.scan((cycle, exp, sample))
                        # It seems like this means it wasn't a real scan?
                        continue
                    
                    # Slight cheat to resolve issue where multiplie scans
                    # wind up with the precise same RT, which causes unexpected
                    # behavior in many use cases.
                    if scaninfo and rt == scaninfo[-1][0]:
                        rt += 0.000001
                    
                    scaninfo.append((rt, precM, (cycle+1, exp+1, sample+1), level, centroid))
                    
        
        if (start_cycle == 0 and stop_cycle == 999999 and
            experiment == None):
            self._scan_info[sample] = scaninfo
        
        return scaninfo
                    
                
                    
    def xic(self, start_time = 0, stop_time = 999999, start_mz = 0, stop_mz = 2000,
            sample = None, experiment = 1, filter = None):
        """
        Get the eXtracted Ion Chromatogram of the given time- and mz-range in
        the given sample data.  If sample is not specified, the previously set
        default (1, unless changed by set_sample()) is used.
        
        Experiment is by default 1 (the experiment of MS1 scans.)
        """
        if filter and filter.strip().lower() not in ['full ms', 'full ms2']:
            raise NotImplementedError, "Filter strings are not compatible with WIFF files. %s" % filter
        
        if not sample:
            sample = self.sample
        
        xic = zip(*self.source.XicByExp(sample-1, experiment-1, float(start_mz), float(stop_mz)))
        return [x for x in xic if start_time <= x[0] <= stop_time]
    
    def tic(self, start_time = None, stop_time = None, sample = None, experiment = None):
        """
        Get the Total Ion Chromatogram of the given time range in the given sample
        data.  This is not dependent on the experiment setting.
        """
    
        if not sample:
            sample = self.sample
            
        tic = zip(*self.source.TicByExp(sample-1, 0))
        if start_time:
            tic = [x for x in tic if x[0] >= start_time]
        if stop_time:
            tic = [x for x in tic if x[0] <= stop_time]
        return tic
    
    def time_range(self, sample = None):
        """
        Gets the total retention time range of the sample.  By default uses the
        MS1 sample (1).
        """
        if sample == None:
            sample = self.sample
        
        return tuple(self.source.GetRTRange(sample-1))
    
    def time_for_scan(self, cycle, experiment = None, sample = None):
        """
        Returns the retention time of a given cycle.
        """
        if sample == None:
            sample = self.sample
        
        #return self.source.GetRTOfScan(sample-1, experiment-1, cycle-1)
        return self._scan_to_rt[(cycle, experiment, sample)]
    
    def scan_for_time(self, rt, experiment = None, sample = None):
        """
        Returns the (cycle, experiment, sample) at a given retention time, if any.
        """        
        if sample == None:
            sample = self.sample
        
        #return int(self.source.GetIndexOfRT(sample-1, experiment-1, rt))
        return self._rt_to_scan[rt][0]
    
    def scan_name_from_scan_time(self, rt, experiment = None, sample = None):
        return self.scan_for_time(rt, experiment, sample)
    
    def scan_time_from_scan_name(self, cycle, experiment = None, sample = None):
        return self.time_for_scan(cycle, experiment, sample)
    
    def scan_range(self, sample = None, experiment = None):
        """
        Returns the beginning and ending cycle numbers of the given sample.
        """
        if sample == None:
            sample = self.sample
            
        return (0, self.source.GetNumCycles(sample-1))
    
    
    def filters(self):
        """
        XCalibur-style MS filter strings, for back-cross-compatibility
        with scripts designed to work with raw.py.
        """
        
        expInfo = {}
        for sample in range(0, self.source.GetSamples()):
            for experiment in range(0, self.source.GetExperiments(sample)):
                expInfo[sample, experiment] = self.source.ExperimentInfo(sample, experiment)
        
        
        self._filters = []
        for rt, mz, (cycle, exp, sample), level, mode in self.scan_info():
            try:
                mzrange = map(int, expInfo[sample-1, exp-1][1:3]) # ...Perhaps?
            except ValueError:
                assert expInfo[sample-1, exp-1][1] == 'MRM' and expInfo[sample-1, exp-1][1] == 'MRM'
                mzrange = 0, 0
            if level == 'MS1':
                levelstr = 'ms'
                locstr = ''
            elif level == 'MS2':
                levelstr = 'ms2'
                locstr = '%.2f@00.00 ' % mz
            else:
                raise Exception, level
            
            detector = expInfo[sample-1, exp-1][3]
            
            filterstr = "%s + %s NSI Full %s %s[%d.00-%d.00]" % (detector, mode, levelstr, locstr,
                                                                 mzrange[0], mzrange[1])
            
            self._filters.append((rt, filterstr))
        
        return self._filters
    
    def scan_modes(self):
        modes = []
        for sample in range(0, self.source.GetSamples()):
            for experiment in range(0, self.source.GetExperiments(sample)):
                modes.append((experiment+1, self.source.ExperimentInfo(sample, experiment)[3]))
        return modes
    
    def headers(self, *etc):
        self._headers = self.scan_info()
        return self._headers
    
    def MRM_info(self):
        values = self.source.GetMRMInfo(self.sample-1, 0)
        keys = ['ExperimentType', 'SpectrumType', 'RawDataType',
                'Polarity', 'IDAType', 'SourceType', 'NumberOfScans']
        return dict(zip(keys, values))
    
    def MRM_channels(self):
        return list(self.source.GetMRMChannels(self.sample-1, 0))
    
    def MRM_scan(self, cycle):
        return list(self.source.MRMScan(self.sample-1, 0, cycle))
    
    def close(self):
        # No close function in the COM object yet.  Should there be one?
        pass