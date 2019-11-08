try:
    from comtypes.client import CreateObject, GetModule
except ImportError as err:
    import platform
    if 'Windows' not in platform.platform():
        pass
    else:
        raise err

from multiplierz import myHome, vprint
from multiplierz.mzAPI import mzFile as baseFile
from multiplierz.mzAPI import mzScan

import os, sys
from math import floor, ceil

if os.path.basename(sys.executable) == 'mzDesktop.exe':
    agilentDir = os.path.join(os.path.dirname(sys.executable),
                              'interface_modules', 'agilentdlls')
else:
    agilentDir = os.path.join(os.path.dirname(__file__), 'agilentdlls')
assert os.path.exists(agilentDir), agilentDir
# File 'MassSpecDataReader.dll' has to be regasm'd before this works.
msdr = GetModule(os.path.join(agilentDir, r'MassSpecDataReader.tlb'))
bc = GetModule(os.path.join(agilentDir, r'BaseCommon.tlb'))
bda = GetModule(os.path.join(agilentDir, r'BaseDataAccess.tlb'))




# Scan type and ion mode are actually binary flags,
# which may cause problems for combined cases.
scanTypeDict = {7951 : "All",
                15 : "AllMS",
                7936 : "AllMSN",
                4 : "HighResolutionScan",
                256 : "MultipleReaction",
                4096 : "NeutralGain",
                2048 : "NeutralLoss",
                1024 : "PrecursorIon",
                512 : "ProductIon",
                1 : "Scan",
                2 : "SelectedIon",
                8 : "TotalIon",
                0 : "Unspecified"}
scanLevelDict = {0 : "All",
                 1 : "ms",
                 2 : "ms2"}
ionModeDict = {32 : "Apci",
               16 : "Appi",
               4 : "CI",
               2 : "EI",
               64 : "ESI",
               1024 : "ICP",
               2048 : "Jetstream",
               8 : "Maldi",
               1 : "Mixed",
               512 : "MsChip",
               128 : "NanoEsi",
               0 : "Unspecified"}
scanModeDict = {1 : "m", # Mixed
                3 : "c", # Peak
                2 : "p", # Profile
                0 : "Unspecified"}
deviceTypeDict = {20 : "ALS",
                  16 : "AnalogDigitalConverter",
                  31 : "BinaryPump",
                  31 : "CANValves",
                  42 : "CapillaryPump",
                  33 : "ChipCube",
                  41 : "CTC",
                  23 : "DiodeArrayDetector",
                  14 : "ElectronCaptureDetector",
                  17 : "EvaporativeLightScatteringDetector",
                  19 : "FlameIonizationDetector",
                  10 : "FlourescenceDetector",
                  18 : "GCDetector",
                  50 : "IonTrap",
                  3 : "IsocraticPump",
                  22 : "MicroWellPlateSampler",
                  1 : "Mixed",
                  13 : "MultiWavelengthDetector",
                  34 : "Nanopump",
                  2 : "Quadrupole",
                  6 : "QTOF",
                  32 : "QuaternaryPump",
                  12 : "RefractiveIndexDetector",
                  5 : "TandemQuadrupole",
                  11 : "ThermalConductivityDetector",
                  40 : "ThermostattedColumnCompartment",
                  4 : "TOF",
                  0 : "Unknown",
                  15 : "VariableWavelengthDetector",
                  21 : "WellPlateSampler"}

ionPolarityDict = {3 : "+-",
                   1 : "-",
                   0 : "+", 
                   2 : "No Polarity"}
desiredModeDict = {'profile':0,
                   'peak':1,
                   'centroid':1,
                   'profileelsepeak':2,
                   'peakelseprofile':3}

class mzFile(baseFile):
    def __init__(self, datafile, **kwargs):
        self.file_type = '.d'
        self.data_file = datafile
        self._filters = None
        
        self.ticObj = None
        self.info = None
        self.noFilter = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.MsdrPeakFilter')
        
        self.source = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.MassSpecDataReader')
        
        s = self.source.OpenDataFile(datafile)
        if not s:
            raise IOError("Error opening %s" % datafile)
        
    def close(self):
        self.source.CloseDataFile()
        
    def time_range(self):
        if not self.ticObj:
            self.ticObj = self.source.GetTIC()
            
        # DEEP IN OBSCURE AND FORGOTTEN COMTYPES DOCUMENTATION, I DISCOVERED
        # THE LONG-LOST METHOD OF BENDING POINTER(IUnknown) TO OBEY HUMAN
        # COMMAND!
        ranges = self.ticObj.AcquiredTimeRange
        assert len(ranges) == 1, "Multiple time ranges per file not currently supported."
        timerange = ranges[0].QueryInterface(bc.IRange)
        return timerange.Start, timerange.End
    
    def scan_range(self):
        # There's probably a more efficient way of doing this.
        info = list(zip(*self.scan_info()))[2]
        return min(info), max(info)
    
    
    def filters(self):
        """
        Thermo-style filter strings for all spectra; used for compatibility with
        various legacy functions.
        """
    
        ionization = self.source.MSScanFileInformation.IonModes
        if not ionization:
            vprint("Could not determine separation/ionization; defaulting to GCMS.")
            separator = 'GC'
        elif ionization & (4|2): # Bitwise OR and AND.
            separator = 'GC'
        else:
            separator = 'TOF'
    
        colEs = self.source.MSScanFileInformation.CollisionEnergy
        if len(colEs) == 1:
            colE = colEs[0]
        else:
            colE = None
    
        if not self._filters:
            self._filters = []
            for rt, mz, index, level, polarity in self.scan_info():
                scanObj = self.source.GetSpectrum_6(index)
                rangeobj = scanObj.MeasuredMassRange.QueryInterface(bc.IRange) # Yep, definitely spectrum-specific.
                
                if colE: # Singular collision energy in the file.
                    energy = colE
                else:
                    energy = float(scanObj.CollisionEnergy)
                
                if level != 'MS1':
                    precstr = '%.4f@%.2f' % (mz, energy)
                else:
                    precstr = ''
                string = "%s MS %s NSI Full ms%s %s[%.2f-%.2f]" % (separator, polarity, 
                                                               int(level[2]) if level != 'MS1' else '',
                                                               precstr,
                                                               (rangeobj.Start),
                                                               (rangeobj.End))
                self._filters.append((rt, string))
        
        return self._filters
        
    def headers(self):
        return self.scan_info()
        
    def scan(self, scan, centroid = None):
        """
        Returns a spectrum from the specified scan index.
        
        If both centroided and profile-mode data are present in the file,
        which is returned can be controlled by setting the mode argument to
        'profile' or 'centroid'. If the requested mode is not present, an
        empty spectrum will be returned. 'ProfileElsePeak' or
        'PeakElseProfile' will return spectrum of the preferred kind if
        present, else the other.
        """
        
        if centroid == None:
            mode = desiredModeDict['PeakElseProfile'.lower()]
        elif isinstance(centroid, str):
            mode = desiredModeDict[centroid.lower()] # Usually 'profile' or 'centroid'.
        else:
            mode = desiredModeDict['centroid' if centroid else 'profile']
        
        scanObj = self.source.GetSpectrum_8(scan, self.noFilter, self.noFilter, mode)
        return list(zip(scanObj.XArray, scanObj.YArray))
    
    
    def cscan(self, scan):
        """
        Calculates a centroided scan from profile-mode data. If a profile
        mode copy of the specified scan is not available, this raises an
        exception; in that case, you can use mzFile.scan() with mode =
        'centroid' to return the machine-centroided scan.
        """
        
        mode = desiredModeDict['profile']
        scanObj = self.source.GetSpectrum_8(scan, self.noFilter, self.noFilter, mode)        
        mzs, ints = scanObj.XArray, scanObj.YArray
        if not mzs:
            raise IOError("Profile data for scan %s not available." % scan)
        
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
        
        
    
    def scan_info(self, start_time = 0, stop_time = 999999, start_mz = 0, stop_mz = 99999):
        if self.info == None:
            self.info = []
            for index in range(1000000): # Largenum.
                infoObj = self.source.GetScanRecord(index)
                rt = infoObj.RetentionTime
                mz = infoObj.MZOfInterest # I *think* this is MZ when applicable?
                if not rt:
                    break
                if not (start_time <= rt <= stop_time and start_mz <= mz <= stop_mz):
                    continue
                
                level = 'MS%d' % infoObj.MSLevel
                #scantype = infoObj.MSScanType
                polarity = infoObj.IonPolarity
                
                self.info.append((rt, mz, index, level,
                                  ionPolarityDict[polarity]))
            if index == 1000000:
                raise IOError("File too large for constant!")
            
        return self.info
    
    def xic(self, start_time = 0, stop_time = None, start_mz = 0, stop_mz = None, filter = None, UV = False):
        if filter:
            assert filter.strip().lower() == 'full ms', 'Thermo-style XIC filters are not supported for Agilent files.'
        
        # A full XIC can be performed with the TIC object that may have been
        # retrieved regardless.
        if self.ticObj and not any([start_time, stop_time, start_mz, stop_mz]):
            return list(zip(self.ticObj.XArray, self.ticObj.YArray))
        
        if stop_time == None:
            stop_time = 999999
        if stop_mz == None:
            stop_mz = 999999
               
        chromFilter = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.BDAChromFilter')
        
        chromFilter.MSLevelFilter = 0 # "All", should perhaps instead be 1 for "ms1"?
        if not UV:
            chromFilter.ChromatogramType = 7 # Extracted-Ion
        else:
            chromFilter.ChromatogramType = 4 # ExtractedWavelength
        chromFilter.SingleChromatogramForAllMasses = True        
        
        mzRange = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.MinMaxRange')
        mzRange.Min = start_mz
        mzRange.Max = stop_mz
        mzRangeIR = mzRange.QueryInterface(bc.IRange)
        chromFilter.IncludeMassRanges = (mzRange,)
        # If THAT works...!
        
        rtRange = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.MinMaxRange')
        rtRange.Min = start_time
        rtRange.Max = stop_time
        rtRangeIR = rtRange.QueryInterface(bc.IRange)
        chromFilter.ScanRange = rtRangeIR
        

        
        xic = self.source.GetChromatogram(chromFilter)[0].QueryInterface(bda.IBDAChromData)
        return list(zip(xic.XArray, xic.YArray))
    
    
    def uv_trace(self):
        nonmsSource = self.source.QueryInterface(msdr.INonmsDataReader)
        nonmsDevs = nonmsSource.GetNonmsDevices()
        return nonmsSource.GetTWC(nonmsDevs[0])
        
    
    def deisotope_scan(self, scan,
                       tolerance_da = 0.0025, tolerance_ppm = 7,
                       max_charge = None, require_peptide_profile = False):
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
        scanObj = self.source.GetSpectrum_6(scan)
        
        deisoFilter = CreateObject(r'Agilent.MassSpectrometry.DataAnalysis.MsdrChargeStateAssignmentFilter')
        if not (tolerance_da == 0.0025 and tolerance_ppm == 7
                and not (max_charge or require_peptide_profile)):
            
            deisoFilter.AbsoluteTolerance = tolerance_da
            if max_charge:
                deisoFilter.LimitMaxChargeState = max_charge
            deisoFilter.RelativeTolerance = tolerance_ppm
            deisoFilter.RequirePeptideLikeAbundanceProfile = require_peptide_profile
            
        self.source.Deisotope(scanObj, deisoFilter) # Void type, not even a success return value.
        
        return list(zip(scanObj.XArray, scanObj.YArray))
        
        
if __name__ == '__main__':
    foo = mzFile(r'\\rc-data1\blaise\ms_data_share\Max\mzStudio\2016-04-13-Equil1.D')
    foo.uv_trace()