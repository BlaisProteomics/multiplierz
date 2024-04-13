import numpy as np
import sqlite3
import os, sys
from ctypes import *
from collections import defaultdict
from itertools import chain
from bisect import bisect_left, bisect_right
import six

if sys.platform[:5] == "win32":
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'brukerlib', "timsdata.dll")
elif sys.platform[:5] == "linux":
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'brukerlib', "libtimsdata.so")
else:
    raise Exception("Unsupported platform.")
assert os.path.exists(libname), libname    
    
try:
    dll = cdll.LoadLibrary(libname)
except:
    print("WARNING: Failed to load Bruker access modules.")
    print("(This can be due to lack of the Visual C++ Redistributable dependency.)")
dll.tims_open.argtypes = [ c_char_p, c_uint32 ]
dll.tims_open.restype = c_uint64
dll.tims_close.argtypes = [ c_uint64 ]
dll.tims_close.restype = None
dll.tims_get_last_error_string.argtypes = [ c_char_p, c_uint32 ]
dll.tims_get_last_error_string.restype = c_uint32
dll.tims_has_recalibrated_state.argtypes = [ c_uint64 ]
dll.tims_has_recalibrated_state.restype = c_uint32
dll.tims_read_scans_v2.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, c_void_p, c_uint32 ]
dll.tims_read_scans_v2.restype = c_uint32

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


class BadIndexNumbers(Exception):
    pass

def throwLastTimsDataError(dll_handle):
    """Throw last TimsData error string as an exception."""

    len = dll_handle.tims_get_last_error_string(None, 0)
    buf = create_string_buffer(len)
    dll_handle.tims_get_last_error_string(buf, len)
    raise RuntimeError(buf.value)



# This is from the example file Bruker provided; will probably want to
# rewrite/replace it if this is going to be released with multiplierz.
def decodeArrayOfStrings (blob):
    if blob is None:
        return None # property not set

    if len(blob) == 0:
        return [] # empty list

    blob = bytearray(blob)
    if blob[-1] != 0:
        raise ValueError("Illegal BLOB contents.") # trailing nonsense

    if sys.version_info.major == 2:
        return str(str(blob), 'utf-8').split('\0')[:-1]
    if sys.version_info.major == 3:
        return str(blob, 'utf-8').split('\0')[:-1]
class TimsData:

    def __init__ (self, analysis_directory, use_recalibrated_state=False):

        if sys.version_info.major == 2:
            if not isinstance(analysis_directory, six.string_types):
                raise ValueError("analysis_directory must be a Unicode string.")
        if sys.version_info.major == 3:
            if not isinstance(analysis_directory, six.string_types):
                raise ValueError("analysis_directory must be a string.")

        self.dll = dll

        self.handle = self.dll.tims_open(
            analysis_directory.encode('utf-8'),
            1 if use_recalibrated_state else 0 )
        if self.handle == 0:
            throwLastTimsDataError(self.dll)

        self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))

        self.initial_frame_buffer_size = 128*64 # may grow in readScans()

    def __del__ (self):
        if hasattr(self, 'handle'):
            self.dll.tims_close(self.handle)         
            
    def __callConversionFunc (self, frame_id, input_data, func):

        if type(input_data) is np.ndarray and input_data.dtype == np.float64:
            # already "native" format understood by DLL -> avoid extra copy
            in_array = input_data
        else:
            # convert data to format understood by DLL:
            in_array = np.array(input_data, dtype=np.float64)

        cnt = len(in_array)
        out = np.empty(shape=cnt, dtype=np.float64)
        success = func(self.handle, frame_id,
                       in_array.ctypes.data_as(POINTER(c_double)),
                       out.ctypes.data_as(POINTER(c_double)),
                       cnt)

        if success == 0:
            throwLastTimsDataError(self.dll)

        return out

    def indexToMz (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_index_to_mz)
        
    def mzToIndex (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_mz_to_index)
        
    def scanNumToOneOverK0 (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_scannum_to_oneoverk0)

    def oneOverK0ToScanNum (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_oneoverk0_to_scannum)

    def scanNumToVoltage (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_scannum_to_voltage)

    def voltageToScanNum (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tims_voltage_to_scannum)

        
    # Output: list of tuples (indices, intensities)
    def readScans (self, frame_id, scan_begin, scan_end):

        # buffer-growing loop
        while True:
            cnt = int(self.initial_frame_buffer_size) # necessary cast to run with python 3.5
            buf = np.empty(shape=cnt, dtype=np.uint32)
            len = 4 * cnt

            required_len = self.dll.tims_read_scans_v2(self.handle, frame_id, scan_begin, scan_end,
                                                    buf.ctypes.data_as(POINTER(c_uint32)),
                                                    len)
            if required_len == 0:
                throwLastTimsDataError(self.dll)

            if required_len > len:
                if required_len > 16777216:
                    # arbitrary limit for now...
                    raise RuntimeError("Maximum expected frame size exceeded.")
                self.initial_frame_buffer_size = required_len / 4 + 1 # grow buffer
            else:
                break

        result = []
        d = scan_end - scan_begin
        for i in range(scan_begin, scan_end):
            npeaks = buf[i-scan_begin]
            indices     = buf[d : d+npeaks]
            d += npeaks
            intensities = buf[d : d+npeaks]
            d += npeaks
            result.append((indices,intensities))

        return result
            

class mzBruker(object):
    def __init__(self, d_directory):
        self.data_file = d_directory
        tdf = os.path.join(d_directory, 'analysis.tdf')
        tdf_bin = os.path.join(d_directory, 'analysis.tdf_bin')
        if not os.path.exists(tdf):
            raise IOError("%s not found." % tdf)
        if not os.path.exists(tdf_bin):
            raise IOError('%s not found.' % tdf_bin)
        
        self.db_conn = sqlite3.connect(tdf)
        self.cur = self.db_conn.cursor()
        try: # Could be a check for unicode type?
            self.source = TimsData(str(d_directory, 'utf-8'), 
                                   use_recalibrated_state=False) # True?
        except TypeError:
            self.source = TimsData(d_directory, 
                                   use_recalibrated_state=False) # True?            

        frame_types = self.dbquery("SELECT Id, MsMsType, Time FROM Frames")
        self.pasef_frames = []
        self.ms1_frames = []
        for f_id, typenum, rt in frame_types:
            if typenum == 8:
                self.pasef_frames.append((f_id, rt))
            elif typenum == 0:
                self.ms1_frames.append((f_id, rt))
            else:
                raise IOError("Unknown frame type: %s" % typenum)
        
    def close(self):
        self.db_conn.close()
        
    def dbquery(self, command):
        results = self.cur.execute(command).fetchall()
        if all([len(x) == 1 for x in results]):
            return [x[0] for x in results]
        else:
            return results
    
    def time_range(self):
        times = [x[1] for x in self.pasef_frames + self.ms1_frames]
        return min(ms1 + pasef), max(ms1 + pasef)

    def frame(self, framenum, start_scan = None, stop_scan = None,
              force = False):
        """
        Returns all points from the specified frame; each point is a 3D
        coordinate (mz, k0, intensity).
        """
        scan_count = self.dbquery("SELECT NumScans FROM Frames WHERE Id=%d" % framenum)[0]
        # Workaround for scan count bug.
        if start_scan and stop_scan:    
            if force:
                if stop_scan >= scan_count:
                    stop_scan = scan_count - 1
                if start_scan >= stop_scan:
                    start_scan = stop_scan - 1
                if start_scan < 0:
                    start_scan = 0                
            else:
                if stop_scan >= scan_count:
                    raise BadIndexNumbers("Scan over scan count.")
                if start_scan >= stop_scan:
                    raise BadIndexNumbers("Scan stop less than scan start.")
                if start_scan < 0:
                    raise BadIndexNumbers("Scan start less than zero.")
            
        if start_scan == None and stop_scan == None:        
            scans = self.source.readScans(framenum, 0, scan_count)
            k0s = self.source.scanNumToOneOverK0(framenum, list(range(scan_count)))
        else:
            scans = self.source.readScans(framenum, start_scan, stop_scan)
            k0s = self.source.scanNumToOneOverK0(framenum, list(range(scan_count)))[start_scan:stop_scan]
            
        assert len(k0s) == len(scans), (len(k0s), len(scans))        
        
        index_groups, int_groups = list(zip(*[x for x in scans]))
        k0_seq = list(chain(*[[k0 for _ in ig] for k0, ig in zip(list(k0s), int_groups)]))
        index_seq = list(chain(*index_groups))
        int_seq = list(chain(*int_groups))
        mz_seq = self.source.indexToMz(framenum, index_seq)
        return list(zip(mz_seq, k0_seq, int_seq))
 
    def frame_int(framenum):
        scans = self.source.readScans(framenum, start_scan, stop_scan)
        index_groups, int_groups = list(zip(*[x for x in scans]))
        return sum(chain(*int_groups))
    
    def pasef_frame(self, framenum, prec_num = None, include_k0 = False):
        """
        Returns scans for distinct precursors from a PASEF frame.
        """
        
        if prec_num:
            subscans = self.dbquery(("SELECT ScanNumBegin, ScanNumEnd, Precursor "
                                     "FROM PasefFrameMsMsInfo WHERE Frame = %d AND Precursor = %d")
                                    % (framenum, prec_num))            
            assert len(subscans) == 1
            partial_frame = self.frame(framenum, subscans[0][0], subscans[0][1])
            if include_k0:
                return partial_frame
            else:
                return [(x[0], x[2]) for x in partial_frame]
        else:
            subscans = self.dbquery(("SELECT ScanNumBegin, ScanNumEnd, Precursor "
                                     "FROM PasefFrameMsMsInfo WHERE Frame = %d")
                                    % framenum)
            
            scans = []
            for start, stop, prec in subscans:
                partial_frame = self.frame(framenum, start, stop)
                if include_k0:
                    scans.append(partial_frame)
                else:
                    scans.append([(x[0], x[2]) for x in partial_frame])
            
            return scans
    
    def pasef_spectra(self, frames, start_scan, stop_scan,
                      aggregate = False, force = False):
        subspectra = []
        for frame in frames:
            
            # Workaround for scan count bug.
            scan_count = self.dbquery("SELECT NumScans FROM Frames WHERE Id=%d" % frame)[0]
            if force:
                if stop_scan >= scan_count:
                    stop_scan = scan_count - 1
                if start_scan >= stop_scan:
                    start_scan = stop_scan - 1
                if start_scan < 0:
                    start_scan = 0                
            else:
                if stop_scan >= scan_count:
                    raise BadIndexNumbers("Scan over scan count.")
                if start_scan >= stop_scan:
                    raise BadIndexNumbers("Scan stop less than scan start.")
                if start_scan < 0:
                    raise BadIndexNumbers("Scan start less than zero.")

            index_arrs, int_arrs = zip(*self.source.readScans(frame, start_scan, stop_scan))
            indexes, ints = map(lambda x: list(chain(*x)), [index_arrs, int_arrs])
            mzs = self.source.indexToMz(frame, indexes)
            subspectra.append(list(zip(mzs, ints)))
        
        if aggregate:
            raise NotImplementedError
        else:
            return list(zip(frames, subspectra))
        
    def precursor_spectrum(self, precursor):
        frames = self.dbquery(("SELECT Frame, ScanNumBegin, ScanNumEnd "
                               "FROM PasefFrameMsMsInfo WHERE Precursor=%d" 
                               % precursor))
        pts = defaultdict(float)
        for framenum, startscan, stopscan in frames:
            for mz, k0, i in self.frame(framenum, startscan, stopscan):
                pts[mz] += i
        return sorted(pts.items())
        
        
    
    def xic(self, start_rt, stop_rt, start_mz, stop_mz, start_k0, stop_k0):
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
            
            # Likewise, scan-to-1/k0 (where "scan" is a subset of a frame which
            # has a variable size) is monotonically decreasing, so scan bounds
            # are equivalent to matching 1/k0 bounds.
            k0s = list(reversed(self.source.scanNumToOneOverK0(framenum, range(scan_count))))
            max_scan = len(k0s) - bisect_left(k0s, start_k0)
            min_scan = len(k0s) - bisect_right(k0s, stop_k0)
            
            frameint = 0
            for indexes, ints in self.source.readScans(framenum, min_scan, max_scan):
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds):
                    frameint += np.sum(ints[in_bounds])
        
            xic.append((rt, frameint))
            
        return xic
    
    def mobiligram(self, start_rt, stop_rt, start_mz, stop_mz, start_k0, stop_k0):
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
            
            for scan_num, (indexes, ints) in enumerate(self.source.readScans(framenum, min_scan, max_scan),
                                                       start = min_scan):
                k0 = k0_lookup[scan_num]
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds):
                    mob[round(k0, 2)] += np.sum(ints[in_bounds])    
        
        return sorted(mob.items())
    
    def frames_mobiligram(self, frame_list, start_mz, stop_mz, start_k0, stop_k0):
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
            
            for scan_num, (indexes, ints) in enumerate(self.source.readScans(framenum, min_scan, max_scan),
                                                       start = min_scan):
                k0 = k0_lookup[scan_num]
                in_bounds = np.where(np.logical_and(indexes >= start_index, indexes <= stop_index))[0]
                if len(in_bounds):                    
                    mob[round(k0, 2)] += np.sum(ints[in_bounds])    
        
        return sorted(mob.items())    
  
    


