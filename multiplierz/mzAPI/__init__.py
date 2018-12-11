# Copyright 2008 Dana-Farber Cancer Institute
# multiplierz is distributed under the terms of the GNU Lesser General Public License
#
# This file is part of multiplierz.
#
# multiplierz is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiplierz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiplierz.  If not, see <http://www.gnu.org/licenses/>.


"""Common API for multiple mass spectrometry instrument file access

mzAPI unifies access to MS data files by referring to scans by time.
mzAPI uses helper files that access the underlying windows libraries for
the respective instrument.
mzAPI currently supports Thermo Finnigan .RAW files, ABI .WIFF files, Agilent .D files,
mzML, and mzURL (web-based data access)

"""

__author__ = 'Jignesh Parikh, James Webber, William Max Alexander'

__all__ = ['mzFile']

import cPickle
import os
import re
import sys

from multiplierz import logger_message

# Caused problems and, really, who needs to open link files?  If necessary
# and you want to deal with win32com import errors, uncomment-out this and the
# use of it below.
#if sys.platform == 'win32':
    #import pythoncom
    #from win32com.shell import shell

    ## windows function to follow links
    #def follow_link(data_file):
        #link = pythoncom.CoCreateInstance (
            #shell.CLSID_ShellLink,
            #None,
            #pythoncom.CLSCTX_INPROC_SERVER,
            #shell.IID_IShellLink
        #)
        #link.QueryInterface(pythoncom.IID_IPersistFile).Load(data_file)

        ## GetPath returns the name and a WIN32_FIND_DATA structure
        ## which we're ignoring. The parameter indicates whether
        ## shortname, UNC or the "raw path" are to be
        ## returned.

        #data_file,_ = link.GetPath(shell.SLGP_UNCPRIORITY)

        #return data_file

#elif 'readlink' in dir(os):
    ## *nix function to follow symlinks.
    ## This is of limited utility until we can open files
    ## on *nix, of course.
    #def follow_link(data_file):
        #res = os.readlink(data_file)
        #if os.path.isabs(res):
            #return res
        #else:
            #return os.path.join(os.path.dirname(data_file), res)


def make_info_file(data_file, **kwargs):
    """Generates a text file with mapping for scan time to proprietary scan index

    """

    if data_file.lower().startswith('http://'):
        raise NotImplementedError('How do you make an info file for a URL?')

    if data_file.lower().endswith('wiff'):
        file_type = 'wiff'
    elif data_file.lower().endswith('raw'):
        file_type = 'raw'
    elif data_file.lower().endswith('mzml'):
        file_type = 'mzml'
    elif data_file.lower().endswith('mzml.gz'):
        file_type = 'mzml'


    if os.path.exists(data_file + '.mzi'):
        os.remove(data_file + '.mzi')

    if file_type == 'mzml':
        import mzML
        mzML.make_info_file(data_file)
    else:
        my_file = mzFile(data_file, **kwargs)
        (start_time, stop_time) = my_file.time_range()
        scan_list = my_file.scan_info(start_time, stop_time)
        my_file.close()

        info_list = []
        for (time, mz, scan_name, scan_type, scan_mode) in scan_list:
            my_dict = {'time': time, 'mz': mz, 'scan_name': scan_name, 'scan_type': scan_type, 'scan_mode': scan_mode}
            info_list.append(my_dict)
        info_list = mzInfoFile(info_list)

        # pickle object
        with open(data_file + '.mzi', 'w') as f:
            cPickle.dump(info_list, f)


class mzInfoFile(tuple):
    """Subclass of tuple object for storing mzInfo.

    Each element of the tuple is a dictionary that stores information about each scan
    """

    def __init__(self, s):
        tuple.__init__(self, s)
        self.index_dict = dict((d["time"], i) for i,d in enumerate(self))

    def field_list(self):
        """Returns keys for dictionary stored in each list element.

        Assumes all dictionaries have same keys.
        """

        return self[0].keys()

    def sort_by_field(self, field=None):
        if field:
            return sorted(self, key=lambda e: e[field])
        else:
            return list(self)

    def filter(self, list_key=None, value_list=None, value_range=None, sort_field=None):
        """All purpose extract data function

        list_key provides the key used to exctract data
        value_list or value_range can be used to provide the values to exctract
        value_list = [val1, val2, val3...]
        value_range = [start_val, stop_val]
        if value_list is provided, range will be ignored

        sort_field option returns the list sorted based on a specified_field
        """

        #  Extract
        temp_list = self.sort_by_field(field = list_key)
        if value_list != None:
            value_list = set(value_list)
            extracted_list = [i for i in temp_list if i[list_key] in value_list]
        elif value_range != None:
            (start_value, stop_value) = sorted(value_range)
            extracted_list = [i for i in temp_list if start_value <= i[list_key] <= stop_value]
        else:
            extracted_list = temp_list

        if sort_field != None:
            extracted_list.sort(key=lambda i: i[sort_field])

        return extracted_list

    def closest(self, key, value):
        closest_item = self[0]
        if key == "time":
            if value in self.index_dict:
                closest_item = self[self.index_dict[value]]
                return closest_item

        return min((i for i in self), key=lambda e: abs(e[key] - value))


class mzScan(list):
    """A subclass of the list object to represent a raw data scan.
    """

    def __init__(self, s, time, mode='p', mz=0.0, z=0):
        '''Create a scan object.

        - s is an iterable of (m/z, intensity) pairs
        - time is the time of acquisition
        - mode is 'p' or 'c' for profile or centroid scans respectively
        - mz is the m/z of the targeted peak (0.0 if not known/applicable)
        - z is the charge of the targeted peak (0 if not known/applicable)
        '''
        list.__init__(self, s)
        self.time = time
        self.mode = mode
        self.mz = mz
        self.z = z

    def peak(self, mz, tolerance):
        '''Returns the max intensity within a tolerance of a target m/z'''
        return max([i for m,i in self if abs(m-mz) <= tolerance] or [0])


class mzFile(object):
    """Base class for access to MS data files"""

    def __init__(self, data_file, **kwargs):
        """Initializes mzAPI and opens a new file of a specified type

        file_type can be 'raw', 'wiff', 'mzml', or 'mzurl'

        Example:
        >>> data_file = 'C:\\Documents and Settings\\User\\Desktop\\example.RAW'
        >>> mz_file = mzAPI.mzFile(data_file)

        """

        import platform
        bitness = platform.architecture()[0]
        if bitness != '64bit':
            if '32bit' in bitness:
                raise Exception, "mzAPI does not support 32-bit Python!"
            else:
                print ("WARNING- System architecture string %s not "
                      "recognized.  Is this 64-bit Windows Python?") % bitness
        
        #if data_file.lower().endswith('.lnk') or os.path.islink(data_file):
            #data_file = follow_link(data_file)

        if not (data_file.lower().startswith('http://') or os.path.exists(data_file)):
            raise IOError, "%s not found." % data_file

        if data_file.lower().startswith('http://'):
            import mzURL
            self.__class__ = mzURL.mzFile
            self.format = 'mzserver'
            mzURL.mzFile.__init__(self, data_file, **kwargs)
        elif data_file.lower().endswith('.wiff'):
            import mzWiff
            self.format = 'wiff'
            if kwargs.get('implicit_mode', False) == True:
                self.__class__ = mzWiff.mzFile_implicit_numbering
                mzWiff.mzFile_implicit_numbering.__init__(self, data_file, **kwargs)
            elif kwargs.get('implicit_mode', True) == False:
                self.__class__ = mzWiff.mzFile_explicit_numbering
                mzWiff.mzFile_explicit_numbering.__init__(self, data_file, **kwargs)
            else:
                import warnings
                warnings.warn('WIFF mzFile mode not specified; defaulting to '
                              'implicit experiment numbering.')
                self.__class__ = mzWiff.mzFile_implicit_numbering
                mzWiff.mzFile_implicit_numbering.__init__(self, data_file, **kwargs)
            #self.__class__ = mzWiff.mzFile
                
            #mzWiff.mzFile.__init__(self, data_file, **kwargs)
        elif data_file.lower().endswith('.raw'):
            import raw
            self.__class__ = raw.mzFile
            self.format = 'raw'
            raw.mzFile.__init__(self, data_file, **kwargs)
        elif (data_file.lower().endswith('.mzml') or
              data_file.lower().endswith('.mzml.gz') or
              data_file.lower().endswith('.mzmlsql')):
            import mzML
            self.__class__ = mzML.mzFile
            self.format = 'mzml'
            mzML.mzFile.__init__(self, data_file, **kwargs)
        elif data_file.lower().endswith('.d'):
            #assert sys.maxsize <= 2**32, "Agilent files are not currently supported in 64-bit Python."
            
            import D
            self.__class__ = D.mzFile
            self.format = 'd'
            D.mzFile.__init__(self, data_file, **kwargs)
        elif data_file.lower().endswith('.t2d'):
            import mzT2D
            self.__class__ = mzT2D.mzFile
            self.format = 't2d'
            mzT2D.mzFile.__init__(self, data_file, **kwargs)
        else:
            raise NotImplementedError, "Can't open %s; extension not recognized." % data_file
        
        
        
        
    def __enter__(self):
        # this method allows mzFiles to be used in with-statements
        # e.g.: with mzFile(file_name) as m: ...
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # at the end of the with-statement, the file is closed
        self.close()

    def close(self):
        """Closes the open MS data file

        Example:
        >>> mz_file.close()

        """
        raise NotImplementedError('Subclasses must implement this method')

    def scan_list(self, start_time=None, stop_time=None, start_mz=0, stop_mz=99999):
        """Gets a list of [(time,mz)] in the time and mz range provided

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_list = mz_file.scan_list(30.0, 35.0, 435.82, 436.00)

        """

        raise NotImplementedError('Subclasses must implement this method')

    def scan_info(self, start_time, stop_time=0, start_mz=0, stop_mz=99999):
        """Gets a list of [(time, mz, scan_name, scan_type, scan_mode)] in the time and mz range provided

        scan_name = number for RAW files, (cycle, experiment) for WIFF files.

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_info = my_peakfile.scan_info(30.0, 35.0, 435.82, 436.00)


        """

        raise NotImplementedError('Subclasses must implement this method')

    def scan_time_from_scan_name(self, scan_name):
        """Gets scan time for wiff (cycle, experiment) tuple or raw scan number

        Example:
        >>> # raw file
        >>> scan_time = mz_file.scan_time_from_scan_name(2165)
        >>> # wiff file
        >>> scan_time = mz_file.scan_time_from_scan_name((1252, 3))

        """

        raise NotImplementedError('Subclasses must implement this method')

    def scan(self, time):
        """Gets scan based on the specified scan time

        The scan is a list of (mz, intensity) pairs.

        Example:
        >>> scan = mz_file.scan(20.035)

        """

        raise NotImplementedError('Subclasses must implement this method')

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        """Generates eXtracted Ion Chromatogram (XIC) for given time and mz range

        The function integrates the precursor intensities for given time and mz range.
        The xic is a list of (time, intensity) pairs.

        Example:
        >>> xic = mz_file.xic(31.4, 32.4, 435.82, 436.00)

        """

        raise NotImplementedError('Subclasses must implement this method')

    def time_range(self):
        """Returns a pair of times corresponding to the first and last scan time

        Example:
        >>> time_range = mz_file.time_range()

        """

        raise NotImplementedError('Subclasses must implement this method')

    def ric(self, *args, **kwargs):
        '''This method is deprecated, use xic() instead.

        Kept for backwards compatibility.'''

        logger_message(40, 'ric() is deprecated: use xic() instead')
        return self.xic(*args, **kwargs)
