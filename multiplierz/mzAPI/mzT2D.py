from multiplierz import vprint
import multiplierz.mzAPI
#from comtypes.client import CreateObject
from win32com.client import Dispatch
import os






class mzFile(multiplierz.mzAPI.mzFile):
    """
    mzFile class for T2D files.  Since T2D, unlike the other supported
    formats, only contains one scan, the API is different: scan() only requires
    a centroid= argument, and info() retrieves a dict of other information
    about the scan.
    
    Since T2D files do not contain MS-level or charge information, 
    """
    def __init__(self, data_file):
        self.datafile = data_file
        #self.source = CreateObject("{7e3450b1-75e7-49b2-9be7-64cbb2458c56}")
        self.source = Dispatch("{7e3450b1-75e7-49b2-9be7-64cbb2458c56}")
        self.source.OpenFile(data_file)
        
        cole, level, _, _, desc, explevel, exppar = self.source.GetInfo
        
        descwords = desc.split()
        filename = descwords[2].strip(',')
        if filename.split('_')[1].upper() == 'MSMS':
            level = 2
            precursor = float(filename.split('_')[2])
        elif filename.split('_')[1].upper() == 'MS':
            level = 1
            precursor = '-'
        else:
            vprint("Unparsable T2D file name; MS level and precursor mass unavailable.")
            level = -1
        
        mzrange = float(descwords[-3].strip('(')), float(descwords[-1].strip(')'))
            
        
        self._info = {'Collision Energy' : float(cole),
                      'MS Level' : int(level),
                      'Precursor' : float(precursor) if precursor != '-' else None,
                      'Range' : mzrange}
        
    
    def scan(self, scan_name = None, centroid = False):
        """
        Get the spectrum from the T2D file.  Since there's only one spectrum,
        scan_name is ignored.
        """
        # Apparently zero-argument methods in COM objects are
        # treated as data.  Odd!
        if centroid:
            return list(zip(*self.source.GetCentroidScan))
        else:
            return list(zip(*self.source.GetScan))
        
    
    def info(self):
        return self._info

