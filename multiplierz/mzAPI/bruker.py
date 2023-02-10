#Contains modified/merged timsdata/tsfdata.py wrappers from Bruker's TDF/TimsData SDK v2.21.0.4

import numpy as np
import sqlite3
import os, sys
from ctypes import *
from collections import defaultdict
from itertools import chain
from bisect import bisect_left, bisect_right
from pathlib import Path
from enum import Enum

if sys.platform[:5] == "win32": libname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'brukerlib', "timsdata.dll")
elif sys.platform[:5] == "linux": libname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'brukerlib', "libtimsdata.so")
else: raise Exception("Unsupported platform.")
assert os.path.exists(libname), libname    
    
try: dll = cdll.LoadLibrary(libname)
except: print("WARNING: Failed to load Bruker access modules.\n(This can be due to lack of the Visual C++ Redistributable dependency.)")

dll.tims_open_v2.argtypes = [ c_char_p, c_uint32, c_uint32 ]
dll.tims_open_v2.restype = c_uint64
dll.tims_close.argtypes = [ c_uint64 ]
dll.tims_close.restype = None
dll.tims_get_last_error_string.argtypes = [ c_char_p, c_uint32 ]
dll.tims_get_last_error_string.restype = c_uint32
dll.tims_has_recalibrated_state.argtypes = [ c_uint64 ]
dll.tims_has_recalibrated_state.restype = c_uint32
dll.tims_read_scans_v2.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, c_void_p, c_uint32 ]
dll.tims_read_scans_v2.restype = c_uint32
MSMS_SPECTRUM_FUNCTOR = CFUNCTYPE(None, c_int64, c_uint32, POINTER(c_double), POINTER(c_float))
dll.tims_read_pasef_msms.argtypes = [ c_uint64, POINTER(c_int64), c_uint32, MSMS_SPECTRUM_FUNCTOR ]
dll.tims_read_pasef_msms.restype = c_uint32
dll.tims_read_pasef_msms_for_frame.argtypes = [ c_uint64, c_int64, MSMS_SPECTRUM_FUNCTOR ]
dll.tims_read_pasef_msms_for_frame.restype = c_uint32
MSMS_PROFILE_SPECTRUM_FUNCTOR = CFUNCTYPE(None, c_int64, c_uint32, POINTER(c_int32))
dll.tims_read_pasef_profile_msms.argtypes = [ c_uint64, POINTER(c_int64), c_uint32, MSMS_PROFILE_SPECTRUM_FUNCTOR ]
dll.tims_read_pasef_profile_msms.restype = c_uint32
dll.tims_read_pasef_profile_msms_for_frame.argtypes = [ c_uint64, c_int64, MSMS_PROFILE_SPECTRUM_FUNCTOR ]
dll.tims_read_pasef_profile_msms_for_frame.restype = c_uint32

dll.tims_extract_centroided_spectrum_for_frame_v2.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, MSMS_SPECTRUM_FUNCTOR, c_void_p ]
dll.tims_extract_centroided_spectrum_for_frame_v2.restype = c_uint32
dll.tims_extract_centroided_spectrum_for_frame_ext.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, c_double, MSMS_SPECTRUM_FUNCTOR, c_void_p ]
dll.tims_extract_centroided_spectrum_for_frame_ext.restype = c_uint32
dll.tims_extract_profile_for_frame.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, MSMS_PROFILE_SPECTRUM_FUNCTOR, c_void_p ]
dll.tims_extract_profile_for_frame.restype = c_uint32

class ChromatogramJob(Structure):
    _fields_ = [
        ("id",c_int64),
        ("time_begin",c_double), ("time_end",c_double),
        ("mz_min",c_double),     ("mz_max",c_double),
        ("ook0_min",c_double),   ("ook0_max",c_double)
    ]
CHROMATOGRAM_JOB_GENERATOR = CFUNCTYPE(c_uint32, POINTER(ChromatogramJob), c_void_p)
CHROMATOGRAM_TRACE_SINK = CFUNCTYPE(c_uint32, c_int64, c_uint32, POINTER(c_int64), POINTER(c_uint64), c_void_p)  
dll.tims_extract_chromatograms.argtypes = [ c_uint64, CHROMATOGRAM_JOB_GENERATOR, CHROMATOGRAM_TRACE_SINK, c_void_p ]
dll.tims_extract_chromatograms.restype = c_uint32

convfunc_argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32 ]

dll.tims_index_to_mz.argtypes = convfunc_argtypes
dll.tims_index_to_mz.restype = c_uint32
dll.tims_mz_to_index.argtypes = convfunc_argtypes
dll.tims_mz_to_index.restype = c_uint32

dll.tims_scannum_to_oneoverk0.argtypes = convfunc_argtypes
dll.tims_scannum_to_oneoverk0.restype = c_uint32
dll.tims_oneoverk0_to_scannum.argtypes = convfunc_argtypes
dll.tims_oneoverk0_to_scannum.restype = c_uint32

dll.tims_scannum_to_voltage.argtypes = convfunc_argtypes
dll.tims_scannum_to_voltage.restype = c_uint32
dll.tims_voltage_to_scannum.argtypes = convfunc_argtypes
dll.tims_voltage_to_scannum.restype = c_uint32

dll.tims_oneoverk0_to_ccs_for_mz.argtypes = [c_double, c_int32, c_double]
dll.tims_oneoverk0_to_ccs_for_mz.restype = c_double

dll.tims_ccs_to_oneoverk0_for_mz.argtypes = [c_double, c_int32, c_double]
dll.tims_ccs_to_oneoverk0_for_mz.restype = c_double


#Strangely, despite having the same arguments/signatures, tims_* and tsf_* calls can behave differently.
#On Windows, the first time tsf data is sent to a tims_ call, it generates an access violation, 
#but ignoring this and running it through a second time returns the expected data.
#For smooth operation, different object dll calls must be used depending on the fileType
dll.tsf_open.argtypes = [ c_char_p, c_uint32 ]
dll.tsf_open.restype = c_uint64
dll.tsf_close.argtypes = [ c_uint64 ]
dll.tsf_close.restype = None
dll.tsf_get_last_error_string.argtypes = [ c_char_p, c_uint32 ]
dll.tsf_get_last_error_string.restype = c_uint32
dll.tsf_has_recalibrated_state.argtypes = [ c_uint64 ]
dll.tsf_has_recalibrated_state.restype = c_uint32
dll.tsf_read_line_spectrum_v2.argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_float), c_int32 ]
dll.tsf_read_line_spectrum_v2.restype = c_int32
dll.tsf_read_line_spectrum_with_width_v2.argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_float), POINTER(c_float), c_int32 ]
dll.tsf_read_line_spectrum_with_width_v2.restype = c_int32
dll.tsf_read_profile_spectrum_v2.argtypes = [ c_uint64, c_int64, POINTER(c_uint32), c_int32 ]
dll.tsf_read_profile_spectrum_v2.restype = c_int32
convfunc_argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32 ]
dll.tsf_index_to_mz.argtypes = convfunc_argtypes
dll.tsf_index_to_mz.restype = c_uint32
dll.tsf_mz_to_index.argtypes = convfunc_argtypes
dll.tsf_mz_to_index.restype = c_uint32

#Convert 1/K0 to CCS for a given charge and mz
def oneOverK0ToCCSforMz(ook0, charge, mz):
    return dll.tims_oneoverk0_to_ccs_for_mz(ook0, charge, mz)

#Convert CCS to 1/K0 for a given charge and mz
def ccsToOneOverK0ToCCSforMz(ccs, charge, mz):
    return dll.tims_ccs_to_oneoverk0_for_mz(ccs, charge, mz)

class PressureCompensationStrategy(Enum):
    NoPressureCompensation = 0
    AnalyisGlobalPressureCompensation = 1
    PerFramePressureCompensation = 2

class TimsData:

    def __init__ (self, analysis_directory, fileType, use_recalibrated_state=False, pressure_compensation_strategy=PressureCompensationStrategy.NoPressureCompensation):
        self.fileType = fileType
        if sys.version_info.major == 2 and not isinstance(analysis_directory, unicode): raise ValueError("analysis_directory must be a Unicode string.")
        if sys.version_info.major == 3 and not isinstance(analysis_directory, str): raise ValueError("analysis_directory must be a string.")
        self.dll = dll
        
        if self.fileType == 'tdf': 
            self.handle = self.dll.tims_open_v2(analysis_directory.encode('utf-8'), 1 if use_recalibrated_state else 0, pressure_compensation_strategy.value )
            if self.handle == 0: _throwLastTimsDataError(self.dll)
            self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))
            self.buffer_size = 128
        if self.fileType == 'tsf': 
            self.handle = self.dll.tsf_open(analysis_directory.encode('utf-8'), 1 if use_recalibrated_state else 0)
            if self.handle == 0: self.throwLastDataError(self.dll)
            self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tsf"))
            self.line_buffer_size = 1024
            self.buffer_size = 1024
            
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
    
    def __del__ (self):
        self.close()
    
    #Throw the last error string as an exception
    #If you are reading this, then your error probably has nothing to do with what was printed!
    #Look back through the stack trace to the method that called this
    def throwLastDataError(dll_handle):
        if self.fileType == 'tsf': 
            len = dll_handle.tsf_get_last_error_string(None, 0)
            buf = create_string_buffer(len)
            dll_handle.tsf_get_last_error_string(buf, len)
        elif self.fileType == 'tdf':
            len = dll_handle.tims_get_last_error_string(None, 0)
            buf = create_string_buffer(len)
            dll_handle.tims_get_last_error_string(buf, len)
        raise RuntimeError(buf.value)
    
    #If already "native" format understood by DLL then avoid extra copy, otherwise convert; validate success
    def callConversionFunc(self, frame_id, input_data, func):
        if type(input_data) is np.ndarray and input_data.dtype == np.float64: in_array = input_data
        else: in_array = np.array(input_data, dtype=np.float64)
        cnt = len(in_array)
        out = np.empty(shape=cnt, dtype=np.float64)
        success = func(self.handle, frame_id, in_array.ctypes.data_as(POINTER(c_double)), out.ctypes.data_as(POINTER(c_double)), cnt)
        if success == 0: self.throwLastDataError(self.dll)
        return out
    
    def close(self):
        if hasattr(self, 'handle') and self.handle is not None:
            if self.fileType=='tsf': self.dll.tsf_close(self.handle)
            elif self.fileType=='tdf': self.dll.tims_close(self.handle)
            self.handle = None
        if hasattr(self, 'conn') and self.conn is not None:
            self.conn.close()
            self.conn = None
    
    def indexToMz(self, frame_id, indices):
        if self.fileType == 'tsf': return self.callConversionFunc(frame_id, indices, self.dll.tsf_index_to_mz)
        elif self.fileType == 'tdf': return self.callConversionFunc(frame_id, indices, self.dll.tims_index_to_mz)
    
    def mzToIndex(self, frame_id, mzs):
        if self.fileType=='tsf': return self.callConversionFunc(frame_id, mzs, self.dll.tsf_mz_to_index)
        elif self.fileType=='tdf': return self.callConversionFunc(frame_id, mzs, self.dll.tims_mz_to_index)
        
    def scanNumToOneOverK0(self, frame_id, scan_nums):
        return self.callConversionFunc(frame_id, scan_nums, self.dll.tims_scannum_to_oneoverk0)

    def oneOverK0ToScanNum(self, frame_id, mobilities):
        return self.callConversionFunc(frame_id, mobilities, self.dll.tims_oneoverk0_to_scannum)

    def scanNumToVoltage(self, frame_id, scan_nums):
        return self.callConversionFunc(frame_id, scan_nums, self.dll.tims_scannum_to_voltage)

    def voltageToScanNum(self, frame_id, voltages):
        return self.callConversionFunc(frame_id, voltages, self.dll.tims_voltage_to_scannum)

    #Outputs tuple of lists depending on specified requestType and self.fileType
    #self.fileType, requestType: 
    #'tsf', 'line'    - (peak_indices, peak_intensities)
    #'tsf', 'width'   - (peak_indices, peak_intensities, widths)
    #'tsf', 'profile' - (intensities)
    #'tdf', 'profile' - (indices, intensities)
    def readData(self, frame_id, requestType, scan_begin=-1, scan_end=-1):
        
        while True:
            cnt = int(self.buffer_size) 
            if self.fileType == 'tsf':
                len = cnt
                if requestType != 'profile': index_buf = np.empty(shape=cnt, dtype=np.float64)
                if requestType == 'width': width_buf = np.empty(shape=cnt, dtype=np.float32)
                if requestType == 'profile': intensity_buf = np.empty(shape=cnt, dtype=np.uint32)
                else: intensity_buf = np.empty(shape=cnt, dtype=np.float32)
                if requestType == 'line': required_len = self.dll.tsf_read_line_spectrum_v2(self.handle, frame_id, index_buf.ctypes.data_as(POINTER(c_double)), intensity_buf.ctypes.data_as(POINTER(c_float)), self.buffer_size)
                elif requestType == 'profile': required_len = self.dll.tsf_read_profile_spectrum_v2(self.handle, frame_id, intensity_buf.ctypes.data_as(POINTER(c_uint32)), self.buffer_size)
                elif requestType == 'width': required_len = self.dll.tsf_read_line_spectrum_with_width_v2(self.handle, frame_id, index_buf.ctypes.data_as(POINTER(c_double)), intensity_buf.ctypes.data_as(POINTER(c_float)), width_buf.ctypes.data_as(POINTER(c_float)), self.buffer_size)
            elif self.fileType == 'tdf' and requestType == 'profile':
                buf = np.empty(shape=cnt, dtype=np.uint32)
                len = 4 * cnt
                required_len = self.dll.tims_read_scans_v2(self.handle, frame_id, scan_begin, scan_end, buf.ctypes.data_as(POINTER(c_uint32)), len)
                
            #Check and grow buffer as needed up to an arbitrary limit of 16777216
            if required_len < 0: self.throwLastDataError(self.dll)
            if required_len > len:
                if required_len > 16777216: raise RuntimeError("Maximum expected frame size exceeded.")
                if self.fileType == 'tsf': self.buffer_size = required_len
                elif self.fileType == 'tdf': self.buffer_size = required_len / 4 + 1
            else: break
        
        if self.fileType == 'tsf':
            if requestType == 'line': return (index_buf[0 : required_len], intensity_buf[0 : required_len])
            elif requestType == 'profile': return intensity_buf[0 : required_len]
            elif requestType == 'width': return (index_buf[0 : required_len], intensity_buf[0 : required_len], width_buf[0 : required_len])
        elif self.fileType == 'tdf' and requestType == 'profile':
            result = []
            d = scan_end - scan_begin
            for i in range(scan_begin, scan_end):
                npeaks = buf[i-scan_begin]
                indices = buf[d:d+npeaks]
                d += npeaks
                intensities = buf[d:d+npeaks]
                d += npeaks
                result.append((indices,intensities))
            return result

    #Read MS/MS spectra according to specified 'readType'
    #'basic'         - read some peak-picked MS/MS spectra for a given list of precursors; returns a dict mapping 'precursor_id' to a pair of arrays (mz_values, area_values).
    #'frame'         - read peak-picked MS/MS spectra for a given frame; returns a dict mapping 'precursor_id' to a pair of arrays (mz_values, area_values).
    #'profile'       - read some "quasi profile" MS/MS spectra for a given list of precursors; returns a dict mapping 'precursor_id' to the profil arrays (intensity_values).
    #'frameProfile'  - read "quasi profile" MS/MS spectra for a given frame; returns a dict mapping 'precursor_id' to the profil arrays (intensity_values).
    def readPasefMsMs (self, readType, frame_id, precursor_list=None):
        if readType == 'basic' or readType =='profile':  precursors_for_dll = np.array(precursor_list, dtype=np.int64)
        result = {}

        @MSMS_SPECTRUM_FUNCTOR
        def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
            result[precursor_id] = (mz_values[0:num_peaks], area_values[0:num_peaks])
        
        @MSMS_PROFILE_SPECTRUM_FUNCTOR
        def callback_for_dll_profile(precursor_id, num_points, intensity_values):
            result[precursor_id] = intensity_values[0:num_points]
        
        if readType == 'basic': rc = self.dll.tims_read_pasef_msms(self.handle, precursors_for_dll.ctypes.data_as(POINTER(c_int64)), len(precursor_list), callback_for_dll)
        elif readType == 'frame': rc = self.dll.tims_read_pasef_msms_for_frame(self.handle, frame_id, callback_for_dll)
        elif readType == 'profile': rc = self.dll.tims_read_pasef_profile_msms(self.handle, precursors_for_dll.ctypes.data_as(POINTER(c_int64)), len(precursor_list), callback_for_dll_profile)
        elif readType == 'frameProfile': rc =  rc = self.dll.tims_read_pasef_profile_msms_for_frame(self.handle, frame_id, callback_for_dll_profile)
        
        if rc == 0: self.source.throwLastDataError(self.dll)

        return result


    #Read peak-picked spectra for a tims frame; returns a pair of arrays (mz_values, area_values).
    def extractCentroidedSpectrumForFrame (self, frame_id, scan_begin, scan_end, peak_picker_resolution=None):
        result = None

        @MSMS_SPECTRUM_FUNCTOR
        def callback_for_dll(precursor_id, num_peaks, mz_values, area_values):
            nonlocal result
            result = (mz_values[0:num_peaks], area_values[0:num_peaks])

        if peak_picker_resolution is None: rc = self.dll.tims_extract_centroided_spectrum_for_frame_v2(self.handle, frame_id, scan_begin, scan_end, callback_for_dll, None)
        else: rc = self.dll.tims_extract_centroided_spectrum_for_frame_ext(self.handle, frame_id, scan_begin, scan_end, peak_picker_resolution, callback_for_dll, None)

        if rc == 0: self.source.throwLastDataError(self.dll)

        return result

    #Read "quasi profile" spectra for a tims frame; returns the profil array (intensity_values).
    def extractProfileForFrame (self, frame_id, scan_begin, scan_end):
        result = None

        @MSMS_PROFILE_SPECTRUM_FUNCTOR
        def callback_for_dll(precursor_id, num_points, intensity_values):
            nonlocal result
            result = intensity_values[0:num_points]

        rc = self.dll.tims_extract_profile_for_frame(self.handle, frame_id, scan_begin, scan_end, callback_for_dll, None)

        if rc == 0: self.source.throwLastDataError(self.dll)

        return result

    #Extract several MS1-only extracted-ion chromatograms, where an iterator jobs streams ChromatogramJob objects (in ascending order of 'time_begin')
    def extractChromatograms (self, jobs, trace_sink):
        
        @CHROMATOGRAM_JOB_GENERATOR
        def wrap_gen(job, user_data):
            try:
                job[0] = next(jobs)
                return 1
            except StopIteration:
                return 2
            except Exception as e:
                # TODO: instead of printing this here, let extractChromatograms throw this
                print("extractChromatograms: generator produced exception ", e)
                return 0

        @CHROMATOGRAM_TRACE_SINK
        def wrap_sink(job_id, num_points, frame_ids, values, user_data):
            try:               
                trace_sink(job_id, np.array(frame_ids[0:num_points], dtype=np.int64), np.array(values[0:num_points], dtype=np.uint64))
                return 1
            except Exception as e:
                # TODO: instead of printing this here, let extractChromatograms throw this
                print("extractChromatograms: sink produced exception ", e)
                return 0

        #The function 'trace_sink' is called for each extracted trace with three arguments: job ID, numpy array of frame IDs ("x axis"), numpy array of chromatogram values ("y axis").
        #For more information, see the documentation of the C-language API of the timsdata DLL.
        unused_user_data = 0
        rc = self.dll.tims_extract_chromatograms(self.handle, wrap_gen, wrap_sink, unused_user_data)

        if rc == 0: throwLastTimsDataError(self.dll)

class mzBruker(object):
    def __init__(self, d_directory):
        self.data_file = d_directory
        tdfFile = os.path.join(d_directory, 'analysis.tdf')
        tsfFile = os.path.join(d_directory, 'analysis.tsf')
        tdfBinFile = os.path.join(d_directory, 'analysis.tdf_bin')
        tsfBinFile = os.path.join(d_directory, 'analysis.tsf_bin')
        
        if os.path.exists(tdfFile) and os.path.exists(tdfBinFile): self.fileType = 'tdf'
        elif os.path.exists(tsfFile) and os.path.exists(tsfBinFile): self.fileType = 'tsf'
        else: raise IOError(".tsf/.tdf and/or .tsf_bin/.tdf_bin file(s) could not found.")
        
        try: self.source = TimsData(str(d_directory, 'utf-8'), self.fileType, use_recalibrated_state=False)
        except TypeError: self.source = TimsData(d_directory, self.fileType, use_recalibrated_state=False)
        
        self.db_conn = self.source.conn
        self.cur = self.db_conn.cursor()
        
        frame_types = self.dbquery("SELECT Id, MsMsType, Time FROM Frames")
        self.pasef_frames = []
        self.ms1_frames = []
        for f_id, typenum, rt in frame_types:
            if typenum == 8: self.pasef_frames.append((f_id, rt))
            elif typenum == 0: self.ms1_frames.append((f_id, rt))
            else: raise IOError("Unknown frame type: %s" % typenum)
        
    def close(self):
        self.db_conn.close()
        
    def dbquery(self, command):
        results = self.cur.execute(command).fetchall()
        if all([len(x) == 1 for x in results]): return [x[0] for x in results]
        else: return results
    
    def time_range(self):
        times = [x[1] for x in self.pasef_frames + self.ms1_frames]
        return min(ms1 + pasef), max(ms1 + pasef)

    #Returns data in list of tuples for frame; (mz, k0, intensity) for .tdf and (mz, intensity) for .tsf
    def scan(self, framenum, mzIntsReturnOnly=False, start_scan = None, stop_scan = None, force = False):
        
        scan_count = self.dbquery("SELECT COUNT(*) FROM Frames")[0]
        
        #Workaround for scan count bug
        if start_scan and stop_scan:    
            if force:
                if stop_scan >= scan_count: stop_scan = scan_count - 1
                if start_scan >= stop_scan: start_scan = stop_scan - 1
                if start_scan < 0: start_scan = 0
            else:
                if stop_scan >= scan_count: raise BadIndexNumbers("Scan over scan count.")
                if start_scan >= stop_scan: raise BadIndexNumbers("Scan stop less than scan start.")
                if start_scan < 0: raise BadIndexNumbers("Scan start less than zero.")
        
        if start_scan == None and stop_scan == None:       
            scans = self.source.readData(framenum, 'profile', 0, scan_count)
            if self.fileType == 'tdf' and not mzIntsReturnOnly: k0s = self.source.scanNumToOneOverK0(framenum, list(range(scan_count)))
        else:
            scans = self.source.readData(framenum, 'profile', start_scan, stop_scan)
            if self.fileType == 'tdf' and not mzIntsReturnOnly: k0s = self.source.scanNumToOneOverK0(framenum, list(range(scan_count)))[start_scan:stop_scan]
        
        #.tdf data should be sorted in ascending m/z order before being returned to match with other formats
        if self.fileType == 'tsf': 
            return self.source.indexToMz(framenum, list(range(len(scans)))), scans
        elif self.fileType == 'tdf': 
            index_groups, int_groups = list(zip(*[x for x in scans]))
            int_seq = list(chain(*int_groups))
            mz_seq = self.source.indexToMz(framenum, list(chain(*index_groups)))
            indexes = np.argsort(mz_seq)
            mz_seq, int_seq = mz_seq[indexes], np.array(int_seq)[indexes]
            if not mzIntsReturnOnly: 
                assert len(k0s) == len(scans), (len(k0s), len(scans))
                k0_seq = list(chain(*[[k0 for _ in ig] for k0, ig in zip(list(k0s), int_groups)]))
                k0_seq = k0_seq[indexes]
                return list(zip(mz_seq, k0_seq, int_seq))
            else: 
                return mz_seq, int_seq
 
    def frame_int(framenum):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
    
        scans = self.source.readData(framenum, 'profile', start_scan, stop_scan)
        if self.fileType == 'tsf': 
            return sum(scans)
        elif self.fileType == 'tdf':
            index_groups, int_groups = list(zip(*[x for x in scans]))
            return sum(chain(*int_groups))
    
    #Returns scans for distinct PASEF frame precursors
    def pasef_frame(self, framenum, prec_num = None, include_k0 = False):

        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
        #New .dll may have changed the needed dbquery
        
        if prec_num:
            subscans = self.dbquery(("SELECT ScanNumBegin, ScanNumEnd, Precursor "
                                     "FROM PasefFrameMsMsInfo WHERE Frame = %d AND Precursor = %d")
                                    % (framenum, prec_num))            
            assert len(subscans) == 1
            partial_frame = self.scan(framenum, subscans[0][0], subscans[0][1])
            if include_k0: return partial_frame
            else: return [(x[0], x[2]) for x in partial_frame]
        else:
            subscans = self.dbquery(("SELECT ScanNumBegin, ScanNumEnd, Precursor "
                                     "FROM PasefFrameMsMsInfo WHERE Frame = %d")
                                    % framenum)
            
            scans = []
            for start, stop, prec in subscans:
                partial_frame = self.scan(framenum, start, stop)
                if include_k0: scans.append(partial_frame)
                else: scans.append([(x[0], x[2]) for x in partial_frame])
            
            return scans
    
    def pasef_spectra(self, frames, start_scan, stop_scan, aggregate = False, force = False):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
    
        subspectra = []
        for framenum in frames:
            mz_seq, int_seq = self.scan(framenum, True, start_scan, stop_scan)
            subspectra.append(list(zip(mz_seq, int_seq)))
        if aggregate: raise NotImplementedError
        else: return list(zip(frames, subspectra))
        
    def precursor_spectrum(self, precursor):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
        #New .dll may have changed the needed dbquery
    
        frames = self.dbquery(("SELECT Frame, ScanNumBegin, ScanNumEnd "
                               "FROM PasefFrameMsMsInfo WHERE Precursor=%d" 
                               % precursor))
        pts = defaultdict(float)
        for framenum, startscan, stopscan in frames:
            if self.fileType == 'tsf':
                for mz, i in self.scan(framenum, False, startscan, stopscan): pts[mz] += i
            elif self.fileType == 'tdf':
                for mz, k0, i in self.scan(framenum, False, startscan, stopscan): pts[mz] += i
            
        return sorted(pts.items())
        
    def xic(self, start_rt, stop_rt, start_mz, stop_mz, start_k0=-1, stop_k0=-1):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
        #New .dll may have changed the needed dbquery
    
        assert start_rt <= stop_rt
        assert start_mz <= stop_mz
        assert start_k0 <= stop_k0
        
        frame_list = self.dbquery(("SELECT Id, Time, NumScans FROM Frames WHERE "
                                   "MsMsType = 0 AND Time >= %d AND Time <= %d")
                                  % (start_rt, stop_rt))  
        
        xic = []
        c = 0
        for framenum, rt, scan_count in frame_list:
            c += 1
            # MZ is a function of index, and vice versa!  Also the function
            # index-to-mz is monotonically increasing, so index bounds are
            # equivalent to matching mz bounds.
            index_bounds = self.source.mzToIndex(framenum, [start_mz, stop_mz])
            start_index, stop_index = list(index_bounds)
            
            
            #If k0 are not specified, use outermost scan bounds
            if self.fileType == 'tsf':
                min_scan, max_scan = 0, self.dbquery("SELECT COUNT(*) FROM Frames")[0]
            
            #Otherwise, since scan-to-1/k0 is monotonically decreasing, scan bounds are equivalent to matching 1/k0 bounds.
            elif self.fileType == 'tdf':
                k0s = list(reversed(self.source.scanNumToOneOverK0(framenum, range(scan_count))))
                max_scan = len(k0s) - bisect_left(k0s, start_k0)
                min_scan = len(k0s) - bisect_right(k0s, stop_k0)
            
            frameint = 0
            for indexes, ints in self.source.readData(framenum, 'profile', min_scan, max_scan):
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds): frameint += np.sum(ints[in_bounds])
        
            xic.append((rt, frameint))
            
        return xic
    
    def mobiligram(self, start_rt, stop_rt, start_mz, stop_mz, start_k0=-1, stop_k0=-1):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
        #Updated .dll may change the needed dbquery; no support for self.fileType=='tsf'
        if self.fileType=='tsf': raise NotImplementedError
    
        assert start_rt <= stop_rt
        assert start_mz <= stop_mz
        assert start_k0 <= stop_k0
        
        frame_list = self.dbquery(("SELECT Id, Time, NumScans FROM Frames WHERE "
                                   "MsMsType = 0 AND Time >= %d AND Time <= %d")
                                  % (start_rt, stop_rt))  
        
        mob = defaultdict(float)
        c = 0
        for framenum, rt, scan_count in frame_list:
            # MZ is a function of index, and vice versa!  Also the function
            # index-to-mz is monotonically increasing, so index bounds are
            # equivalent to matching mz bounds.
            index_bounds = self.source.mzToIndex(framenum, [start_mz, stop_mz])
            start_index, stop_index = list(index_bounds)
            
            # Likewise, scan-to-1/k0 (where "scan" is a subset of a frame which
            # has a variable size) is monotonically decreasing, so scan bounds
            # are equivalent to matching 1/k0 bounds.
            k0s = list(reversed(self.source.scanNumToOneOverK0(framenum, range(scan_count))))
            max_scan = len(k0s) - bisect_left(k0s, start_k0)
            min_scan = len(k0s) - bisect_right(k0s, stop_k0)
            k0_lookup = dict(enumerate(reversed(k0s)))
            
            for scan_num, (indexes, ints) in enumerate(self.source.readData(framenum, 'profile', min_scan, max_scan), start = min_scan):
                k0 = k0_lookup[scan_num]
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds):
                    mob[round(k0, 2)] += np.sum(ints[in_bounds])    
        
        return sorted(mob.items())
    
    def frames_mobiligram(self, frame_list, start_mz, stop_mz, start_k0=-1, stop_k0=-1):
    
        #WARNING/TODO: THIS CODE HAS NOT BEEN VALIDATED SINCE BEING UPDATED
        #Updated .dll may change the needed dbquery; no support for self.fileType=='tsf'
        if self.fileType=='tsf': raise NotImplementedError
    
        # Sometimes convenient to target a specific set of frames instead of an
        # RT range.
        assert start_mz <= stop_mz
        assert start_k0 <= stop_k0        
    
        frame_rt_counts = self.dbquery("SELECT Id, Time, NumScans FROM Frames WHERE " + 
                                       " OR ".join(["Id = %d" % f for f in frame_list]))
        assert len(frame_rt_counts) == len(frame_list)
    
        mob = defaultdict(float)
        c = 0
        for framenum, rt, scan_count in frame_rt_counts:
            # MZ is a function of index, and vice versa!  Also the function
            # index-to-mz is monotonically increasing, so index bounds are
            # equivalent to matching mz bounds.
            index_bounds = self.source.mzToIndex(framenum, [start_mz, stop_mz])
            start_index, stop_index = list(index_bounds)
            
            # Likewise, scan-to-1/k0 (where "scan" is a subset of a frame which
            # has a variable size) is monotonically decreasing, so scan bounds
            # are equivalent to matching 1/k0 bounds.
            k0s = list(reversed(self.source.scanNumToOneOverK0(framenum, range(scan_count))))
            max_scan = len(k0s) - bisect_left(k0s, start_k0)
            min_scan = len(k0s) - bisect_right(k0s, stop_k0)
            k0_lookup = dict(enumerate(reversed(k0s)))
            
            for scan_num, (indexes, ints) in enumerate(self.source.readData(framenum, 'profile', min_scan, max_scan), start = min_scan):
                k0 = k0_lookup[scan_num]
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds): mob[round(k0, 2)] += np.sum(ints[in_bounds])    
        
        return sorted(mob.items())    
