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

import csv
#csv.field_size_limit(100000000)
import os

from multiplierz.mzReport import ReportReader, ReportWriter, ReportEntry, default_columns
from multiplierz import logger_message

class CSVReportReader(ReportReader):
    '''The CSV implementation of the mzReport class. Python's CSV module
    is actually fine for what we need, so this is mostly a wrapper of
    its functionality.'''
    def __init__(self, file_name):
        if not os.path.exists(file_name):
            raise IOError("No such file: '%s'" % os.path.basename(file_name))
        self.file_name = file_name
        self.fh = open(self.file_name, 'rb')
        self.csv = csv.reader(self.fh)
        self.columns = self.csv.next()
        
        if len(self.columns) == 1 and '\t' in self.columns[0]:
            # Retry with tab-delimiting dialect
            
            self.fh.close()
            self.fh = open(self.file_name, 'rb')
            self.csv = csv.reader(self.fh, dialect = 'excel-tab')
            self.columns = self.csv.next()[:]            
            
            

    def __iter__(self):
        self.fh.seek(0) # start at beginning of the file
        self.csv.next() # ignore headers
        for row in self.csv:
            # will convert default column values to correct type
            yield ReportEntry(columns=self.columns, values=row)

    def close(self):
        self.fh.close()


class CSVReportWriter(ReportWriter):
    '''The CSV implementation of the mzReport class. Python's CSV module
    is actually fine for what we need, so this is mostly a wrapper of
    its functionality.'''
    def __init__(self, file_name, columns=None, default_columns=False):
        if default_columns:
            self.columns = default_columns + (columns or [])
        else:
            self.columns = columns

        for i, col in enumerate(self.columns):
            if not all(map(lambda x: ord(x) < 128, col)):
                newcol = ''.join([x for x in col if ord(x) < 128])
                print "WARNING: %s has non-utf-8 characters; changing to %s" % (col, newcol)
                self.columns[i] = newcol

        if len(self.columns) > len(set(c.lower() for c in self.columns)):
            raise ValueError("Redundant columns: column headers must be unique")

        self.file_name = file_name
        self.fh = open(self.file_name, 'wb')
        self.csv = csv.writer(self.fh)

        self.csv.writerow(self.columns)

    def write(self, row, metadata=None):
        # error checking: want one value per column and nothing more
        # if row is a dict, keys should be in lower-case
        if len(row) < len(self.columns):
            raise ValueError('Must have values for each column')
        elif len(row) > len(self.columns):
            raise ValueError('Too many values')
        elif isinstance(row,dict):
            row = dict((k.lower(),v) for k,v in row.items())
            if not all(k.lower() in row for k in self.columns):
                raise ValueError('Value dictionary does not match column headers')

        if isinstance(row,dict):
            self.csv.writerow([row[col.lower()] for col in self.columns])
        else:
            self.csv.writerow(row)

        if metadata:
            logger_message(40, 'Metadata not supported in CSV files')

    def add_image(self, column, image):
        logger_message(40, 'Images are not supported in CSV files')

    def close(self):
        self.fh.close()
