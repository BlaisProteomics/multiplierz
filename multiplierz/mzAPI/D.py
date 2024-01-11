from multiplierz import myHome, vprint
from multiplierz.mzAPI import mzFile as baseFile
from multiplierz.mzAPI import mzScan

import clr
import os, sys
import numpy as np
from math import floor, ceil

dll_path = 'agilentdlls'
sys.path += [os.path.join(os.path.dirname(__file__), dll_path)]
dlls = ['MassSpecDataReader', 'BaseCommon', 'BaseDataAccess']
for dll in dlls: clr.AddReference(dll)
import Agilent
from Agilent.MassSpectrometry.DataAnalysis import *

#pythonnet does not auto-transfer methods/vars from implemented classes.
#Only transferred those relevant to classes used, though more functionality
#could potentially be added. Note: .dll files inspectable with ILDASM via VS
for item in dir(IMsdrDataReader): 
    try: setattr(MassSpecDataReader, item, getattr(IMsdrDataReader, item))
    except: pass
for item in dir(IBDAChromFilter): 
    try: setattr(BDAChromFilter, item, getattr(IBDAChromFilter, item))
    except: pass
for item in dir(IMsdrChargeStateAssignmentFilter): 
    try: setattr(MsdrChargeStateAssignmentFilter, item, getattr(IMsdrChargeStateAssignmentFilter, item))
    except: pass

class mzFile(baseFile):
    
    def __init__(self, datafile, **kwargs):
        self.file_type= '.d'
        self.data_file = datafile
        self._filters = None
        self.ticObj = None
        self.info = None
        self.scanRange = None
        self.noFilter = MsdrPeakFilter()
        self.source = MassSpecDataReader()
        s = self.source.OpenDataFile(self.source, datafile)
        if not s: raise IOError("Error opening %s" % datafile)
        
    def close(self):
        self.source.CloseDataFile(self.source)
        
    def time_range(self):
        if not self.ticObj: self.ticObj = self.source.GetTIC(self.source)
        assert ticObj.AcquiredTimeRange.Length == 1, "Multiple time ranges are not supported"
        return ranges.GetValue(0).Min, ranges.GetValue(0).Max
    
    def scan_range(self):
        return 0, self.source.FileInformation.MSScanFileInformation.TotalScansPresent
    
    def scan_info(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None):
        if self.info == None:
            self.info = []
            for index in range(self.source.FileInformation.MSScanFileInformation.TotalScansPresent):
                infoObj = self.source.GetScanRecord(self.source, index)
                rt = infoObj.RetentionTime
                mz = infoObj.MZOfInterest
                if not rt: break
                if start_time != None and rt <= start_time: continue
                if stop_time != None and rt >= stop_time: continue
                if start_mz != None and mz <= start_mz: continue
                if stop_mz != None and mz >= stop_mz: continue
                level = 'MS%d' % int(infoObj.MSLevel)
                polarity = str(infoObj.IonPolarity)
                #scanType = str(infoObj.MSScanType)
                self.info.append((rt, mz, index, level, polarity))
        return self.info
    
    def headers(self):
        return self.scan_info()
    
    #Thermo-style filter strings for all spectra; legacy compatibility
    def filters(self):
        ionization = str(self.source.MSScanFileInformation.IonModes)

        if not ionization:
            vprint("Could not determine separation/ionization; defaulting to GCMS.")
            separator = 'GC'
        elif ionization and (ionization == "CI" or ionization == "EI"):
            separator = 'GC'
        else:
            separator = 'TOF'

        colEs = self.source.MSScanFileInformation.CollisionEnergy
        if colEs.Length == 1:
            colE = colEs.GetValue(0)
        else:
            colE = None

        if not self._filters:
            self._filters = []
            for rt, mz, index, level, polarity in self.scan_info():
                scanObj = self.source.GetSpectrum(self.source, index)

                if colE: # Singular collision energy in the file.
                    energy = colE
                else:
                    energy = scanObj.CollisionEnergy

                if level != 'MS1':
                    precstr = '%.4f@%.2f' % (mz, energy)
                else:
                    precstr = ''
                string = "%s MS %s NSI Full ms%s %s[%.2f-%.2f]" % (separator, polarity, 
                                                               int(level[2]) if level != 'MS1' else '',
                                                               precstr,
                                                               (scanObj.MeasuredMassRange.Start),
                                                               (scanObj.MeasuredMassRange.End))
                self._filters.append((rt, string))

        return self._filters
    
    def scan(self, index, mode=None):
        
        """
        Returns a spectrum from the specified scan index, where type is
        controlled by mode argument; defaults to prioiritize Peak (Centroid) 
        and then try Profile. If the requested mode is not present, an empty 
        spectrum  will be returned. Alternative modes are 'profile', 'centroid' 
        or  'profileElsePeak'. 
        """
        
        if mode != None: mode = mode.lower()
        if mode == None or mode == 'peakelseprofile':
            mode = DesiredMSStorageType.PeakElseProfile
        elif mode == 'profileelsepeak':
            mode = DesiredMSStorageType.ProfileElsePeak
        elif mode == 'profile':
            mode = DesiredMSStorageType.Profile
        elif mode == 'peak':
            mode = DesiredMSStorageType.Peak
        else: 
            return []
        scanObj = self.source.GetSpectrum(self.source, index, self.noFilter, self.noFilter, mode)
        return list(zip(scanObj.XArray, scanObj.YArray))
    
    def cscan(self, index):
        
        """
        Calculates a centroided scan from profile-mode data. If a profile
        mode copy of the specified scan is not available, this raises an
        exception; in that case, you can use mzFile.scan() with mode =
        'centroid' to return the machine-centroided scan.
        """
        
        mode = DesiredMSStorageType.Profile
        scanObj = self.source.GetSpectrum(self.source, index, self.noFilter, self.noFilter, mode)
        mzs, ints = list(scanObj.XArray), list(scanObj.YArray)
        if not mzs:
            raise IOError("Profile data for scan index %s not available." % index)
        
        threshold = average(ints)
        peaks = []
        peak = []
        for pt in zip(mzs, ints):
            if pt[1] > threshold:
                peak.append(pt)
            elif peak:
                #centroid = average(zip(*peak)[0]), average(zip(*peak)[1])
                centroid = average(list(zip(*peak))[0], weights = list(zip(*peak))[1]), max(list(zip(*peak))[1])
                peaks.append(centroid)
                peak = []
                
        return peaks
    
    def xic(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None, filter=None, UV=False):
        if filter:
            assert filter.strip().lower() == 'full ms', 'Thermo-style XIC filters are not supported for Agilent files.'

        #A full XIC can be performed with an existing TIC object
        if self.ticObj and not any([start_time, stop_time, start_mz, stop_mz]):
            return list(zip(self.ticObj.XArray, self.ticObj.YArray))

        if start_time == None: start_time = 0
        if stop_time == None: stop_time = 999999
        if start_mz == None: start_mz = 0
        if stop_mz == None: stop_mz = 999999

        chromFilter = BDAChromFilter()

        chromFilter.set_MSLevelFilter(chromFilter, MSLevel.MS) #Alt value is MSLevel.MSMS
        if not UV:
            chromFilter.set_ChromatogramType(chromFilter, ChromType.ExtractedIon)
        else:
            chromFilter.set_ChromatogramType(chromFilter, ChromType.ExtractedWavelength)
        chromFilter.set_SingleChromatogramForAllMasses(chromFilter, True)

        mzRange = MinMaxRange()
        mzRange.set_Min(start_mz)
        mzRange.set_Max(stop_mz)
        chromFilter.set_IncludeMassRanges(chromFilter, (mzRange, ))

        rtRange = MinMaxRange()
        rtRange.set_Min(start_time)
        rtRange.set_Max(stop_time)
        chromFilter.set_ScanRange(chromFilter, rtRange)

        xic = self.source.GetChromatogram(self.source, chromFilter).Get(0)
        return list(zip(xic.XArray, xic.YArray))

    def uv_trace(self):
        
        #Cannot verify functionality, so leaving a couple potential methods here
        nonmsDevs  = self.source.GetNonmsDevices()
        if not nonmsDevs.Length: raise IOError("No NonmsDevices were available")
        return self.source.GetTWC(nonmsDevs[0])

        #nonmsSource = INonmsDataReader(self.source)
        #nonmsDevs = nonmsSource.GetNonmsDevices()
        #if not nonmsDevs.Length: raise IOError("No NonmsDevices were available")
        #return nonmsSource.GetTWC(nonmsDevs[0])
        
    def deisotope_scan(self, index, tolerance_da=0.0025, tolerance_ppm=7, max_charge=None, require_peptide_profile=False):
        """
        The Agilent MassHunter DAC has the neat feature of including its own
        deisotoping algorithm!  This function uses that to return the specified
        scan in deisotoped form.
        
        tolerance_da (Daltons) and tolerance_ppm (Parts-per-million) are
        ADDED TOGETHER to obtain the total tolerance value for each peak.
        
        max_charge can be set to an integer to only consider isotopic envelopes of
        charge equal to or less; 'None' performs no charge filtering.
        
        require_peptide_profile filters based on isotopic sequences having
        the relative intensity profile caused by standard relative isotopic
        abundances.
        """
        
        scanObj = self.source.GetSpectrum(self.source, index)
        deisoFilter = MsdrChargeStateAssignmentFilter()
        
        deisoFilter.set_AbsoluteTolerance(deisoFilter, tolerance_da)
        deisoFilter.set_RelativeTolerance(deisoFilter, tolerance_ppm)
        if max_charge: 
            deisoFilter.set_LimitMaxChargeState(deisoFilter, True)
            deisoFilter.set_MaximumChargeState(deisoFilter, max_charge)
        else: 
            deisoFilter.set_LimitMaxChargeState(deisoFilter, False)
        deisoFilter.set_RequirePeptideLikeAbundanceProfile(deisoFilter, require_peptide_profile)
        
        self.source.Deisotope(self.source, scanObj, deisoFilter)
        return list(zip(scanObj.XArray, scanObj.YArray))
        
        
#Not sure if this was intended for testing as file is not available
#if __name__ == '__main__':
#    foo = mzFile(r'\\rc-data1\blaise\ms_data_share\Max\mzStudio\2016-04-13-Equil1.D')
#    foo.uv_trace()