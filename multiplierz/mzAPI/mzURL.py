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

#import lxml.etree as etree
import os
import cStringIO
import tempfile

import pycurl

from multiplierz.mzAPI import mzScan, mzFile as mzAPImzFile #, mzInfoFile

# some question of how to create an info file for these..
#def make_info_file(data_file):
    #'''Makes an info file for an mzML file, so that we don't have
    #to iterate through the whole thing every time.
    #'''

    ## building our own table is about as fast as reading through
    ## the file and extracting theirs, so to simplify code we'll
    ## just do it our way (because we need that function anyway)
    #if os.path.exists(data_file + '.mzi'):
        #os.remove(data_file + '.mzi')

    #m = mzFile(data_file)
    #m._build_info_scans()

    #fh = open(data_file + '.mzi', 'wb')
    #cPickle.dump(m._info_scans, fh)
    #fh.close()

def check_mzURL(mz_server, file_name):
    '''Checks if an mzURL actually exists.

    mz_server should be the base URL of the server
    file_name is the name of the specific file (without its extension)
    '''

    if mz_server[-1] == '/':
        mz_server = mz_server[:-1]

    # Handle to libcurl object
    crl = pycurl.Curl()

    # set some general options
    crl.setopt(pycurl.FOLLOWLOCATION, True)
    crl.setopt(pycurl.URL, str(mz_server + '/files.txt'))

    output = cStringIO.StringIO()
    crl.setopt(pycurl.WRITEFUNCTION, output.write)

    try:
        for i in range(5):
            #print 'check mzurl %d' % i
            crl.perform()
            if output.getvalue():
                break
    except pycurl.error as e:
        return False

    for f in output.getvalue().splitlines():
        if os.path.splitext(f)[0].lower() == file_name.lower():
            return True
    else:
        return False


class mzFile(mzAPImzFile):
    """Class for access to mzServer URLs--i.e. remote access to data.
    Essentially this is just a pycurl session that is requesting pages
    from the server it's been pointed to. It does some caching so that
    it's not too bad.
    """

    def __init__(self, data_file, verbose=False, **kwargs):
        self.file_type = 'mzurl'
        # strip off the final slash, if it exists
        if data_file[-1] == '/':
            data_file = data_file[:-1]
        # Likewise, html or other madness.
        if any([data_file.lower().endswith(x) for x in ['html', 'raw', 'wiff']]):
            data_file = ".".join(data_file.split(".")[:-1])
        self.data_file = data_file # actually a URL to a file
        self.verbose = verbose

        self._scans = None # cache of scan_info results for the whole file

        # A string with the name and path of an appropriate temp file
        # (varies by platform)
        fd, self.cookie_file_name = tempfile.mkstemp(text=True)
        os.close(fd)

        # Handle to libcurl object
        self.crl = pycurl.Curl()

        # set some general options
        self.crl.setopt(pycurl.COOKIEFILE, self.cookie_file_name)
        self.crl.setopt(pycurl.COOKIEJAR, self.cookie_file_name)
        self.crl.setopt(pycurl.FOLLOWLOCATION, True)
        self.crl.setopt(pycurl.VERBOSE, verbose)

        self.output = cStringIO.StringIO()
        self.crl.setopt(pycurl.WRITEFUNCTION, self.output.write)

        # how would you store an info file?
        #if os.path.exists(data_file + '.mzi'):
            #self._info_file = data_file + '.mzi'
            #info_fh = open(self._info_file)
            #self._info_scans = cPickle.load(info_fh)
            #info_fh.close()
        #else:
            #self._info_file = None

    def login(self, login, password):
        raise NotImplementedError("Secure mzServer doesn't exist (and may never)")

    def close(self):
        '''Close the 'file', i.e. pycurl session.

        If self.verbose is True, writes the output to stdout
        '''
        self.crl.close()
        if self.verbose:
            print((self.output.getvalue()))
        self.output.close()
        os.unlink(self.cookie_file_name)

    def scans(self):
        '''Gets the scan_info and caches it, to save extra trips to the server'''

        if not self._scans:
            self.crl.setopt(pycurl.HTTPGET, True)
            self.crl.setopt(pycurl.URL, str(self.data_file + '/scan_info'))

            # the scan_name is the URL of that scan, which is just the time
            scan_url = self.data_file + '/scans/%s'

            response = cStringIO.StringIO()
            self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

            #self.crl.perform()
            for i in range(5):
                #print 'scans %d' % i
                self.crl.perform()
                if response.getvalue():
                    break

            scans = response.getvalue().splitlines()

            self._scans = []
            for scan_line in scans:
                time,mz,scan_type,scan_mode = scan_line.split()
                self._scans.append((float(time),
                                    float(mz),
                                    scan_url % time,
                                    scan_type.upper(),
                                    scan_mode.lower()))

            self._scans.sort()

        return self._scans

    def scan_list(self, start_time=None, stop_time=None, start_mz=0, stop_mz=99999):
        """Gets a list of [(time,mz)] in the time and mz range provided

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_list = my_peakfile.scan_list(30.0, 35.0, 435.82, 436.00)

        """
        if not start_time or not stop_time:
            (file_start_time, file_stop_time) = self.time_range()
        if not start_time:
            start_time = file_start_time
        if not stop_time:
            stop_time = file_stop_time

        return [(t,mz) for t,mz,sn,st,sm in self.scans()
                if start_time <= t <= stop_time and (st == 'MS1' or start_mz <= mz <= stop_mz)]

    def scan_info(self, start_time, stop_time=0, start_mz=0, stop_mz=99999):
        """Gets a list of [(time, mz, scan_name, scan_type, scan_mode)] in the time and mz range provided

        scan_name = the URL of the scan (goes to the HTML version)

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_info = my_peakfile.scan_info(30.0, 35.0, 435.82, 436.00)
        """
        if stop_time == 0:
            stop_time = start_time

        return [(t,mz,sn,st,sm) for t,mz,sn,st,sm in self.scans()
                if start_time <= t <= stop_time and (st == 'MS1' or start_mz <= mz <= stop_mz)]

    def scan_time_from_scan_name(self, scan_name):
        # the scan name is just a URL, and the last part is the time.
        if scan_name.endswith('/'):
            return float(scan_name.split('/')[-2])
        else:
            return float(scan_name.split('/')[-1])

    def scan(self, time):
        scan_time,mz,scan_name,st,scan_mode = min(self.scans(), key=lambda si: abs(time - si[0]))

        self.crl.setopt(pycurl.HTTPGET, True)
        self.crl.setopt(pycurl.URL, str(scan_name + '.txt'))

        response = cStringIO.StringIO()
        self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

        for i in range(5):
            #print 'scan %d' % i
            self.crl.perform()
            if response.getvalue():
                break

        scan = response.getvalue().splitlines()

        # how to get charge for an mzURL scan?
        return mzScan([tuple(float(v) for v in s.split()) for s in scan],
                      scan_time, mode=scan_mode, mz=mz)

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        xic_url = str(self.data_file + ('/ric/%s-%s/%s-%s.txt' % (start_time, stop_time,
                                                                  start_mz, stop_mz)))

        self.crl.setopt(pycurl.HTTPGET, True)
        self.crl.setopt(pycurl.URL, xic_url)

        response = cStringIO.StringIO()
        self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

        for i in range(5):
            #print 'xic %d' % i
            self.crl.perform()
            if response.getvalue():
                break

        scan = response.getvalue().splitlines()

        return [tuple(float(v) for v in s.split()) for s in scan]

    def time_range(self):
        return (self.scans()[0][0],
                self.scans()[-1][0])

    @property
    def server_name(self):
        '''Returns the name of the server this file is hosted on.'''
        # return the base, or the 'files' part?
        return '/'.join(self.data_file.split('/')[:-2])
