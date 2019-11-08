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


__author__ = 'James Webber, Jignesh Parikh'

import os
import sys

__all__ = ['reader', 'writer']

# this dictionary is for converting between shorthand abbreviations
# and the longer more useful names used as column headers in reports
multiplierzHeaders = {
    'prot_rank': 'Protein Rank',
    'prot_db': 'Protein Database',
    'acc': 'Accession Number',
    'prot_desc': 'Protein Description',
    'prot_mass': 'Protein Mass',
    'prot_matches': 'Protein Matches',
    'prot_score': 'Protein Score',
    'seq': 'Peptide Sequence',
    'var_mods': 'Variable Modifications',
    'mz': 'Experimental mz',
    'charge': 'Charge',
    'pred_mr': 'Predicted mr',
    'delta': 'Delta',
    'pep_score': 'Peptide Score',
    'pep_rank': 'Peptide Rank',
    'start': 'Start Position',
    'end': 'End Position',
    'res_before': 'Preceding Residue',
    'res_after': 'Following Residue',
    'miss': 'Missed Cleavages',
    'spec_desc': 'Spectrum Description',
    'pep_query': 'Query',
    'ion_str': 'Ion String',
    'peak_int': 'Peak Intensity',
    'peak_area': 'Peak Area',
    'peak_width': 'Peak Width (sec)',
    'peak_comment': 'Peak Comment',
    'peak_time': 'Peak Time',
    'scan_time': 'MS2 Time',
    'prot_cov': 'Protein Coverage'
    }

# These are the default columns that every mzReport should be aware of.
# They should not be changed, since these columns and their sequence is
# still hardcoded in mzSearch.mascot.interface.output_peptides() .
default_columns = ['Protein Rank', 'Protein Database', 'Accession Number',
                   'Protein Description', 'Protein Mass', 'Protein Matches', 'Protein Score',
                   'Peptide Sequence', 'Variable Modifications', 'Experimental mz', 'Charge',
                   'Predicted mr', 'Delta', 'Peptide Score', 'Peptide Rank',
                   'Start Position', 'End Position', 'Preceding Residue', 'Following Residue',
                   'Missed Cleavages', 'Spectrum Description', 'Query']

# the types of the default columns
default_types = dict((k.lower(),t)
                     for k,t in zip(default_columns,
                                    (int, str, str, float,
                                     int, float, str, str,
                                     float, int, float,
                                     float, float, int,
                                     int, int, str, str,
                                     float, str, int)))

# types for the additional columns added by extract_peaks
default_types.update((k.lower(),t)
                     for k,t in zip(('MS2 Time', 'Peak Time', 'Peak Intensity',
                                     'Peak Width (sec)', 'Peak Comment', 'Peak Area',
                                     'Rep114','Rep115','Rep116','Rep117',
                                     'Protein Coverage'),
                                    (float, float, float, float, str, float,
                                     float, float, float, float, float)))

## 'report entry' class
#class ReportEntry(dict):
    #'''Class to represent a single entry in a report. Can be indexed
    #by name (case-insensitive column header) or by position (integer).
    #The class is mutable, but mutating it won't change the report it
    #came from--just the entry itself, which can be added to a new report.

    #On initialization, this class will convert default column values to
    #their correct type
    #'''
    #def __init__(self, columns, values, typecast = True):
        #'''On initialization, the columns are converted to lower-case
        #and values are converted to their correct type if possible.
        #'''
        ##super(dict, self).__init__()
        ###self.columns = [str(c).lower() for c in columns]
        ##self.columns = None
        
        #if isinstance(values, dict):
            #values = [values[c.lower()] for c in columns]
            
        #self.columns = [c.lower() for c in columns if c]
        #self.case_columns = [c for c in columns if c]
        #for k,v in zip(self.columns, values):
            #if typecast and v and (k in default_types):
                #try:
                    #self[k] = default_types[k](v)
                #except ValueError as err:
                    #self[k] = v
                    #continue
                    ##raise ValueError, ("Mismatched values in Mascot report; \n have"
                                       ##"you selected the correct version of Mascot in Settings?")
            #else:
                #self[k] = v if v != None else ''

    #def __getitem__(self, key):
        #if isinstance(key, int):
            #return self[self.columns[key]]
        #elif isinstance(key, slice):
            #return [self[k] for k in self.columns[key]]
        #else:
            #return super(ReportEntry, self).__getitem__(key.lower())

    #def get(self, key, default):
        ## I think this would only make sense for key access, not list-like access.
        #if key.lower() in self: 
            #return super(ReportEntry, self).__getitem__(key.lower())
        #else:
            #return default
    
    ## These all assume dict-nature.
    #def keys(self):
        #return self.case_columns
    
    #def items(self):
        #return [(k, self[k.lower()]) for k in self.case_columns]
    
    #def iteritems(self):
        #return ((k, self[k.lower()]) for k in self.case_columns)
        
    #def __delitem__(self, key):
        #del self[key]
        #self.case_columns = [x for x in self.case_columns if x.lower() != key.lower()]
        #self.columns.remove(key.lower())
    
    #def iterkeys(self):
        #return (x for x in self.case_columns)
    
    
        
        
    #def __setitem__(self, key, value):
        #if isinstance(key, int):
            #self[self.columns[key]] = value
        #elif isinstance(key, slice):
            #if len(self.columns[key]) != len(value):
                #raise KeyError("Assignment to a slice must preserve length")
            #for c,v in zip(self.columns[key], value):
                #self[c] = v
        #else:
            #try:
                #if key.lower() not in self.columns:
                    #self.columns.append(key.lower())
                    #self.case_columns.append(key)
                #super(ReportEntry, self).__setitem__(key.lower(), value)
            #except AttributeError:
                ## If there's no columns, then this is some pickle/multiprocessing
                ## situation.  Best to just play along.
                #super(ReportEntry, self).__setitem__(key.lower(), value)
                
    #def __contains__(self, key):
        #return super(ReportEntry, self).__contains__(key.lower())

    ##def __getattr__(self, name):
        ##if (name in multiplierzHeaders
            ##and multiplierzHeaders[name].lower() in self):
            ##return self[multiplierzHeaders[name].lower()]
        ##elif name.lower() in self:
            ##return self[name.lower()]
        ##else:
            ##raise AttributeError("'ReportEntry' object has no attribute '%s'" % name)

    #def __setattr__(self, name, value):
        #if (name in multiplierzHeaders
            #and multiplierzHeaders[name].lower() in self):
            #self[multiplierzHeaders[name].lower()] = value
        #elif name.lower() in self:
            #self[name.lower()] = value
        #else:
            #super(ReportEntry, self).__setattr__(name,value)
            
class ReportEntry(dict):
    def __init__(self, columns, values, typecast = True):
        super(dict, self).__init__()
        if isinstance(values, dict):
            self.update(values)
        else:
            self.update(list(zip(columns, values)))


# A report is a fairly arbitrary set of data--any way of representing a list of
# rows, defined by some list of column names. The only constraint is that each
# entry conform to those headers.

# Reading and writing reports are separate tasks, and thus have separate classes
# to represent them. This constrains workflow somewhat, but that isn't necessarily
# bad--it forces code to follow a linear plan



class ReportReader(object):
    """Base class for reading multiplierz reports"""
    def __init__(self, file_name, *args, **kwargs):
        self.file_name = file_name
        # this method should define self.columns so that
        # other methods can access it

    def __iter__(self):
        raise NotImplementedError("Subclasses must override this method")

    def __enter__(self):
        # returns the open reader, no other set up necessary
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # closes the reader
        self.close()

    def close(self):
        '''For closing the report before all rows have been read'''
        raise NotImplementedError("Subclasses must override this method")


class ReportWriter(object):
    """Base class for writing multiplierz reports"""
    def __init__(self, file_name, columns=None, default_columns=False, *args, **kwargs):
        """Instantiates a Report object for access to a multiplierz report
        given a file name and file type.
        """
        self.file_name = file_name
        self.columns = columns

    def __enter__(self):
        # returns the open writer, no other set up necessary
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # closes the writer
        self.close()

    def write(self, row, metadata=None):
        """Adds a row of values to the report.

        'row' should be a sequence or dictionary of values corresponding
        to the columns of the report.

        'metadata' should be a list of tuples (col,type,data) where col is
        the column name, type is one of 'image', 'prot coverage', 'formula'.
        Not all types will be supported in all formats (csv won't support any)
        """
        raise NotImplementedError("Subclasses must override this method")

    def add_image(self, column, image):
        """Adds an image associated with the specified column to the last row
        of the report. Not supported by CSV.

        """
        raise NotImplementedError("Subclasses must override this method")

    def close(self):
        """Closes the file."""
        raise NotImplementedError("Subclasses must override this method")


# imports are down here to prevent circular import problem


#from multiplierz import myData, myTemp

# reader factory returns a ReportReader
def reader(report_file, **kwargs):
    """Returns a Report object for access to a multiplierz report
    given a file name and file type, or None if the type is invalid

    Valid file types are 'xls'/'xlsx', 'csv', 'mzid', and 'mzd' (our sqlite file extension)
    """

    file_type = os.path.splitext(report_file)[1][1:].lower()

    report_file = os.path.abspath(report_file)
        
    if not os.path.exists(report_file): raise IOError("%s not found." % report_file)

    if file_type == 'xls' or file_type == 'xlsx':
        import multiplierz.mzReport.mzSpreadsheet as mzSpreadsheet
        return mzSpreadsheet.XLSheetReader(report_file, **kwargs)
    elif file_type == 'csv' or file_type == 'gz':
        import multiplierz.mzReport.mzCSV as mzCSV
        return mzCSV.CSVReportReader(report_file, **kwargs)
    elif file_type == 'mzd':
        import multiplierz.mzReport.mzDB as mzDB
        return mzDB.SQLiteReader(report_file, **kwargs)
    elif file_type == 'mzid':
        import multiplierz.mzReport.mzIdentML as mzIdentML
        return mzIdentML.mzIdentMLReader(report_file, **kwargs)
    else:
        raise IOError('File extension %s is not supported.' % file_type)


# writer factory returns a ReportWriter
def writer(report_file, columns=None, default_columns=False, delayed_write = False, **kwargs):
    """Returns a Report writer for creating a multiplierz report
    given a file name and file type, or None if the type is invalid

    Valid file types are 'xls'/'xlsx', 'csv', and 'mzd' (the multiplierz
    sqlite file extension.)
    
    If delayed_write is True, no file-writing operations are performed until the
    .close() method is called.  (In which case, it is important to call this!)
    """
    
    report_file = os.path.abspath(report_file)

    file_type = os.path.splitext(report_file)[1][1:].lower()

    if file_type not in ('xls', 'xlsx', 'csv', 'mzd', 'gz'):
        raise NotImplementedError("Multiplierz cannot write file of type %s" % file_type)
        
        
    if delayed_write == True:
        return delayedWriter(report_file, columns, default_columns, **kwargs)

    if file_type == 'xls' or file_type == 'xlsx':
        import multiplierz.mzReport.mzSpreadsheet as mzSpreadsheet
        return mzSpreadsheet.XLSheetWriter(report_file,
                                           columns=columns,
                                           default_columns=default_columns,
                                           **kwargs)
    elif file_type == 'csv' or file_type == 'gz':
        import  multiplierz.mzReport.mzCSV as mzCSV
        return mzCSV.CSVReportWriter(report_file,
                                     columns=columns,
                                     default_columns=default_columns,
                                     **kwargs)
    elif file_type == 'mzd':
        import  multiplierz.mzReport.mzDB as mzDB
        return mzDB.SQLiteWriter(report_file,
                                 columns=columns,
                                 default_columns=default_columns,
                                 **kwargs)
    elif file_type == 'mzid':
        raise NotImplementedError('mzIdentML file writing is not supported.')
    else:
        raise IOError("File appears to be an invalid type: '.%s' writing is not supported" % file_type)


class delayedWriter(ReportWriter):
    def __init__(self, report_file, columns, default_columns, **kwargs):
        self.file_name = report_file
        if columns:
            self.columns = columns[:]
        elif default_columns:
            self.columns = default_columns + (columns or [])
        self.extraArgs = kwargs
        
        
        if len(self.columns) < len(set([x.lower() for x in self.columns])):
            import defaultdict
            counts = defaultdict(int)
            for col in self.columns:
                counts[col] += 1
                
            raise ValueError("Column titles appear more than once: %s" % [k for (k, v) in list(counts.items()) if v > 1])
    
        self.data = []
        
    def write(self, *args, **kwargs):
        self.data.append((args, kwargs))
    
    def close(self):
        trueWriter = writer(self.file_name, self.columns, delayed_write = False, **self.extraArgs)
        for args, kwargs in self.data:
            trueWriter.write(*args, **kwargs)
        trueWriter.close()
