from comtypes.client import CreateObject
from multiplierz.mzAPI import mzScan, mzFile as mzAPImzFile


__author__ = 'William Max Alexander'


class mzFile(mzAPImzFile):
    def __init__(self, file_name, *args, **kwargs):
        self.file_type = 'raw'
        self.data_file = file_name
        
        try:
            self.source = CreateObject("{10729396-43ee-49e5-aa07-85f02292ac70}")
        except WindowsError as err:
            print "RawReader.dll not found in registry."
            raise err
        self.source.OpenRawFile(file_name)
        
        self._filters = None
        self._scaninfo = None
        
        # A bunch of functions are pure C#.
        self.scan_range = self.source.scan_range
        self.time_range = self.source.time_range
        self.time_for_scan = self.source.time_from_scan
        self.scan_for_time = self.source.scan_from_time
        self.scanForTime = self.source.scan_from_time
        self.timeForScan = self.source.time_from_scan
        self.scan_time_from_scan_name = self.source.time_from_scan
        
        
    
    def scan(self, scan_number, centroid = False):
        if isinstance(scan_number, float): # Taken as a scan time.
            scan_number = self.source.scan_from_time(scan_number)
        if centroid:
            return [x[:2] for x in self.source.centroid_scan(scan_number)]
        else:
            return list(self.source.profile_scan(scan_number))
    
    def lscan(self, scan_number):
        scan = self.source.centroid_scan(scan_number)
        if len(scan[0]) < 4:
            raise IOError, "Full lscan data not available for scan %s" % scan_number
        return [x[:4] for x in scan]
    
    def rscan(self, scan):
        return list(self.source.centroid_scan(scan))
    
    def average_scan(self, startscan, stopscan, filter):
        return list(self.source.average_scan(startscan, stopscan, filter))
    
    def filters(self):
        if not self._filters:
            filterlist = zip(self.source.GetAllFilterInfoTimes(),
                             self.source.GetAllFilterInfo())
            self._filters = [(time, string) for time, string in filterlist
                             if time and string]
            
        return self._filters
    
    def extra_info(self, scan):
        info = {}
        for key, val in self.source.get_extra_scan_info(scan):
            key = key.strip(':')
            try:
                info[key] = float(val)
            except ValueError:
                info[key] = val
        return info
    
    def scanInjectionTime(self, scan):
        keys, vals = zip(*self.source.get_extra_scan_info(scan))
        try:
            return float(vals[keys.index('Ion Injection Time (ms):')])
        except ValueError:
            return float(vals[keys.index('Ion Injection Time (ms):')].replace(',', ''))
    
    def scanPrecursor(self, scan):
        keys, vals = zip(*self.source.get_extra_scan_info(scan))
        return (float(vals[keys.index('Monoisotopic M/Z:')]),
                float(vals[keys.index('Charge State:')]))
        
    
    def scan_info(self, start_time=0, stop_time=None, start_mz=0, stop_mz=100000):
        if not self._scaninfo:
            self._scaninfo = []
            
            for scan in range(*self.scan_range()):
                info = self.source.GetFilterInfoForScan(scan)
                time = self.source.time_from_scan(scan)
                self._scaninfo.append((time, float(info[1]), scan,
                                       info[0].upper() if info[0].upper() != 'MS' else 'MS1',
                                       'p' if info[2] == 'Profile' else 'c'))
                
        if start_time:
            start_scan = self.scanForTime(start_time)
        else:
            start_scan = self.scan_range()[0]
        if stop_time:
            stop_scan = self.scanForTime(stop_time)
        else:
            stop_scan = self.scan_range()[1]        
        return [x for x in self._scaninfo if 
                start_scan <= x[2] <= stop_scan and
                start_mz <= x[1] <= stop_mz]
    
    def headers(self):
        return self.scan_info()
               
    def xic(self, start_time = -1, stop_time = -1, start_mz = -1, stop_mz = -1, filter=None):
        if start_time == -1 and stop_time == -1 and start_mz == -1 and stop_mz == -1 and filter==None:
            return self.tic()
        if not filter:
            filter = 'Full ms '
        if start_time != -1:
            start_time = self.scanForTime(start_time)
        if stop_time != -1:
            stop_time = self.scanForTime(stop_time)
        return list(self.source.get_xic(start_time, stop_time, start_mz, stop_mz, filter))
    
    def tic(self):
        return list(self.source.get_tic())
    
    def close(self):
        pass # Check for memory leaks!
    

    
    def centroid(self, scan, *foo, **bar):
        print "mzFile.centroid is deprecated, use mzFile.scan(..., centroid = True) instead."
        return self.scan(scan, centroid = True)
    

                
    
#if __name__ == '__main__':
    #import os
    #import glob
    #import subprocess
    #os.chdir(r'C:\Users\Max\Documents\Visual Studio 2017\Projects\ConsoleApp3\ConsoleApp3\bin\Debug')
    
    #netDir = sorted(glob.glob('c:/Windows/Microsoft.NET/Framework%s/v[34]*/RegAsm.exe' % '64'))[-1]
    #dllpath = r'C:\Users\Max\Documents\Visual Studio 2017\Projects\ConsoleApp3\ConsoleApp3\bin\Debug\ConsoleApp3.dll'
    #ret = subprocess.call([netDir, dllpath, "/tlb", "/codebase"])    