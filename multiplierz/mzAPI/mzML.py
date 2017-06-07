import multiplierz.mzAPI
from multiplierz.mzml import mzmlToSqlite, demarshal
import os, sys

import tempfile
import sqlite3 as sqlite
import warnings






class mzFile(multiplierz.mzAPI.mzFile):
    """
    mzAPI access to mzFile files. Loads everything from the XML into a SQLite
    database, which is a slow startup procedure but allows fast random-access
    to the data without tremendous RAM usage. For once-through sequential
    processing, this is NOT the way to go- try the iterator functions in
    multiplierz.mzML instead.
    
    The intermediary SQLite files can be saved for later use by calling
    the .save_cache() function; these files are fully portable between
    computers and mostly portable between software versions, and also
    slightly smaller than their source mzML.  (Usually larger than
    the source raw data, however.)
    """
    
    def __init__(self, data_file, *etc, **etcetc):
        if data_file.lower().endswith('mzml'):
            flptr, flname = tempfile.mkstemp(suffix = '.mzmlsql')
            #flptr.close()
            print "Parsing mzML to SQLite... (this may take some time.)"
            mzmlToSqlite(data_file, flname)
            self.connection = sqlite.connect(flname)
            self.cursor = self.connection.cursor()
            self.cachename = flname
        elif data_file.lower().endswith('mzmlsql'):
            print "Opening previously prepared SQLite file..."
            self.connection = sqlite.connect(data_file)
            self.cursor = self.connection.cursor()
            self.cachename = data_file
        else:
            raise NotImplementedError, '%s is not a recognized file type for mzML reader!' % data_file    
    
    def save_cache(self, save_file):
        self.connection.close()
        os.move(self.cachename, save_file)
        self.connection = sqlite.connect(save_file)
        self.cursor = self.connection.cursor()
        self.cachename = save_file
        
    
    def scan_info(self, start_mz = None, stop_mz = None, start_rt = None, stop_rt = None):
        if not any([start_mz, stop_mz, start_rt, stop_rt]):
            self.cursor.execute('SELECT ind, mz, rt FROM spectra')
        else:
            self.cursor.execute('SELECT ind, mz, rt FROM spectra WHERE mz >= %s AND mz <= %s AND rt >= %s AND rt <= %s'
                                % (start_mz, stop_mz, start_rt, stop_rt))
        manifest = list(self.cursor.fetchall())
        manifest = [(r, m, i, 'MS1' if not m else 'MS2', None) for i, m, r in manifest]
        return sorted(manifest)
    
    def scan(self, scan_name, centroid = None):
        """
        Returns the scan with the specified index.  Set centroid to True for
        centroided data, or False for uncentroided data.
        """
        self.cursor.execute('SELECT data FROM spectra WHERE ind=%s' % str(scan_name))
        scandata = demarshal(self.cursor.fetchone()[0])
        
        if centroid and not scandata['centroid']:
            raise NotImplementedError, "Requires general centroiding function!"
        elif centroid == False and scandata['centroid']:
            warnings.warn("Non-centroid data requested but only centroided data available.", RuntimeWarning)
        
        return scandata['spectrum']
    
    
    def scan_for_time(self, time, tolerance = 0.0001):
        # Float inprecision spoils the exact value of RT keys when they
        # get into Python scope, requiring the addition of a tolerance factor.
        self.cursor.execute('SELECT ind, rt FROM spectra WHERE rt >= %s AND rt <= %s'
                            % (time - (tolerance/2), time + (tolerance/2)))
        return min(self.cursor.fetchall(), key = lambda x: abs(x[1] - time))[0]
    
    def time_for_scan(self, scan):
        self.cursor.execute('SELECT rt FROM spectra WHERE ind=%s' % scan)
        return self.cursor.fetchone()[0]
    
    def scan_time_from_scan_name(self, scan):
        return self.time_for_scan(scan)
    
    def scan_range(self):
        self.cursor.execute('SELECT MAX(ind) FROM spectra')
        topscan = self.cursor.fetchone()[0]
        self.cursor.execute('SELECT MIN(ind) FROM spectra')
        botscan = self.cursor.fetchone()[0]
        return botscan, topscan
    
    def time_range(self):
        self.cursor.execute('SELECT MAX(rt) FROM spectra')
        toprt = self.cursor.fetchone()[0]
        self.cursor.execute('SELECT MIN(rt) FROM spectra')
        botrt = self.cursor.fetchone()[0]
        return botrt, toprt
    
    
    def xic(self, start_time = None, stop_time = None, start_mz = None, stop_mz = None):
        if any([start_time, stop_time, start_mz, stop_mz]):
            raise NotImplementedError, "Uncertain how to deal with custom XICs from mzML data."
        else:
            self.cursor.execute('SELECT data FROM chromato WHERE ind=0 AND startmz=0 AND stopmz=0 AND startrt=0 AND stopmz=0')
            return demarshal(self.cursor.fetchone()[0])[1]
        