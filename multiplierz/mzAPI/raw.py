import clr
import sys, os

dll_path = 'thermo_dlls'
dlls = ['ThermoFisher.CommonCore.Data',
        'ThermoFisher.CommonCore.RawFileReader',
        'ThermoFisher.CommonCore.BackgroundSubtraction',
        'ThermoFisher.CommonCore.MassPrecisionEstimator',
        'System.Collections']

# DLLs cannot be loaded from network mounts!

sys.path += [os.path.join(os.path.dirname(__file__), dll_path)]
#base_path = os.getcwd()
#os.chdir(dll_path)
for dll in dlls:
    clr.AddReference(dll)
#os.chdir(base_path)

from ThermoFisher.CommonCore.Data import ToleranceUnits
from ThermoFisher.CommonCore.Data.Business import (ChromatogramSignal, ChromatogramTraceSettings,
                                                   DataUnits, Device, GenericDataTypes, SampleType,
                                                   Scan, TraceType, Range)
from ThermoFisher.CommonCore.Data.FilterEnums import IonizationModeType, MSOrderType
from ThermoFisher.CommonCore.Data.Interfaces import (IChromatogramSettings, IScanEventBase,
                                                     IScanFilter, RawFileClassification)
from ThermoFisher.CommonCore.MassPrecisionEstimator import PrecisionEstimate
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter

from ThermoFisher.CommonCore.Data import Extensions 

def besttype(x):
    try:
        return float(x)
    except ValueError:
        return x

class mzFile(object):
    def __init__(self, filename, *etc, **etcetc):
        self.source = RawFileReaderAdapter.FileFactory(filename)
        if self.source.InstrumentCount > 1:
            print("%s has data from %d instruments, but only one is supported." 
                  % (os.path.basename(filename), self.source.InstrumentCount))
        self.source.SelectInstrument(Device.MS, 1)
        
        self.data_file = filename
        self.file_type = 'RAW'
        
        self._filters = None
        self._info = None
        
        self.time_from_scan = self.source.RetentionTimeFromScanNumber
        self.scan_from_time = self.source.ScanNumberFromRetentionTime
        self.timeForScan = self.time_from_scan
        self.scanForTime = self.scan_from_time
    
    def scan_range(self):
        return (self.source.RunHeaderEx.FirstSpectrum,
                self.source.RunHeaderEx.LastSpectrum)
    
    def time_range(self):
        return (self.source.RunHeaderEx.StartTime,
                self.source.RunHeaderEx.EndTime)
    
    def filters(self):
        # In RawReader it's an additional O(N) step to get the RTs to
        # match each filter.  It's recommended to use .filter_list() in
        # future code, which doesn't include RTs.
        if not self._filters:
            scanrange = self.scan_range()
            scans = list(range(scanrange[0], scanrange[1]+1))
            self._filters = []
            for scan in scans:
                rt = self.time_from_scan(scan)
                filt = self.specific_filter(scan)
                self._filters.append((rt, filt))
        
        return self._filters
    
    def scan_info(self, *etc, **etcetc):
        if not self._info:
            self._info = []
            scanrange = self.scan_range()
            for scannum in range(scanrange[0], scanrange[1]+1):
                time = self.time_from_scan(scannum)
                info = self.source.GetFilterForScanNumber(scannum)
                level = 'MS%d' % info.MSOrder
                if level != 'MS1':
                    mz = info.GetMass(0) # I guess >0 is for MS3 etc?
                else:
                    mz = 0.0
                pol = 'p' if info.Polarity == 1 else '-'
                self._info.append((time, mz, scannum, level, pol))
                
        return self._info
    
    def headers(self):
        return self.scan_info()
    
    def filter_list(self, *etc, **etcetc):
        if not self._filters:
            self._filters = [IScanFilter(x).ToString() for x in 
                             self.source.GetFilters()]
        return self._filters
    
    def specific_filter(self, scannum):
        return IScanFilter(self.source.GetFilterForScanNumber(scannum)).ToString()
    
    def scan(self, scannum, centroid = False):
        # There's also apparently a "Scan.FromFile" method of getting scan data?
        
        scan_stats = self.source.GetScanStatsForScanNumber(scannum)
        # Does IsCentroidScan indicate that profile data is not available?
        if centroid or scan_stats.IsCentroidScan:
            if centroid == False:
                raise IOError("No profile data for scan %s" % scannum)
            
            stream = self.source.GetCentroidStream(scannum, False)
            if stream.Masses is not None and stream.Intensities is not None:
                return list(zip(stream.Masses, stream.Intensities))
            else:
                # Fall back on "profile" mode, which seems to usually turn
                # out centroid data for some reason.  The format is confused.
                scan = self.source.GetSegmentedScanFromScanNumber(scannum, scan_stats)
                return list(zip(scan.Positions, scan.Intensities))                
        
        else: # Profile-only scan.
            scan = self.source.GetSegmentedScanFromScanNumber(scannum, scan_stats)
            scan = list(zip(scan.Positions, scan.Intensities))
            return scan
        
    def lscan(self, scannum):
        stream = self.source.GetCentroidStream(scannum, False)
        return list(zip(stream.Masses, stream.Intensities, 
                        stream.Noises, stream.Charges))
    
    def rscan(self, scannum):
        stream = self.source.GetCentroidStream(scannum, False)
        return list(zip(stream.Masses, stream.Intensities, 
                        stream.Resolutions))    
    
    def average_scan(self, start_scan, stop_scan, filter = 'Full ms', centroid = False):
        average_scan = Extensions.AverageScansInScanRange(self.source,
                                                          start_scan, stop_scan,
                                                          filter)
        
        #options = self.source.DefaultMassOptions()
        #options.ToleranceUnits = ToleranceUnits.ppm
        #options.Tolerance = 5.0
        #scanfilter = self.source.GetFilterForScanNumber(start_scan)
        #average_scan = self.source.AverageScansInScanRange(start_scan, stop_scan,
                                                           #scanfilter, options)
                                                           
        if not average_scan.HasCentroidStream:
            raise IOError("Could not retrieve average scan %d - %d" % (start_scan,
                                                                        stop_scan))
            # May still be able to get a profile mode scan out of it...?
        if centroid:
            return list(zip(average_scan.CentroidScan.Masses,
                            average_scan.CentroidScan.Intensities))
        else:
            return list(zip(average_scan.SegmentedScan.Positions,
                            average_scan.SegmentedScan.Intensities))
    
    def xic(self, start_time = 0, stop_time = 99999,
            start_mz = 0, stop_mz = 99999,
            filter = 'Full ms'):
        start_scan, stop_scan = list(map(self.scan_from_time, [start_time, stop_time]))
        #settings = ChromatogramTraceSettings(TraceType.BasePeak)
        settings = ChromatogramTraceSettings(filter, [Range.Create(start_mz, stop_mz)])
        
        xic_data = self.source.GetChromatogramData([settings], start_scan, stop_scan)
        xic_trace = ChromatogramSignal.FromChromatogramData(xic_data)[0]
        return list(zip(xic_trace.Times, xic_trace.Intensities))
    
    def extra_info(self, scan):
        trailer = self.source.GetTrailerExtraInformation(scan)
        # Labels come with trailing ':'s that are annoying; those are removed.
        return {l[:-1]:besttype(v) for l, v in zip(trailer.Labels, trailer.Values)}
    
    def scanInjectionTime(self, scan):
        trailer = self.source.GetTrailerExtraInformation(scan)
        return float(trailer.Values[list(trailer.Labels)
                                    .index('Ion Injection Time (ms):')])
    
    def scanPrecursor(self, scan):
        trailer = self.source.GetTrailerExtraInformation(scan)
        return (float(trailer.Values[list(trailer.Labels).index('Monoisotopic M/Z:')]),
                int(trailer.Values[list(trailer.Labels).index('Charge State:')]))
                              
    def close(self):
        self.source.Dispose()
    
            
