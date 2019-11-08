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

import sqlite3
import os
import re
import cPickle

from multiplierz.mzReport import ReportReader, ReportWriter, ReportEntry, default_columns, default_types
from multiplierz import logger_message


# SQLite type conversion. The trailing space is for insertion into statements
sqlite_types = {int: ' INTEGER',
                int: ' INTEGER',
                float: ' REAL',
                str: ' TEXT',
                str: ' TEXT'}

# this function converts a scan into a binary blob for sqlite
# it simply uses Python's pickle functionality.
def adapt_data(data):
    return sqlite3.Binary(cPickle.dumps(data, protocol=2))

sqlite3.register_adapter(tuple, adapt_data)

# the function to convert a blob back into a scan, by un-pickling
def convert_data(data):
    return cPickle.loads(str(data))

sqlite3.register_converter('pickled', convert_data)


class SQLiteReader(ReportReader):
    '''The SQLite implementation of the mzReport reader class.'''

    def __init__(self, file_name, table_name=None, sheet_name=None):
        if sheet_name and not table_name:
            table_name = sheet_name
        elif sheet_name and table_name and (sheet_name != table_name):
            raise IOError("Doubly-specified table name: %s and %s" % (sheet_name, table_name))
        self.table_name = table_name or 'PeptideData'

        if not os.path.exists(file_name):
            raise IOError("No such file: '%s'" % os.path.basename(file_name))

        self.file_name = file_name
        self.conn = sqlite3.connect(file_name, detect_types=sqlite3.PARSE_COLNAMES)

        cursor = self.conn.execute("select name from sqlite_master where type='table' or type='view'")
        self.tables = tuple(t[0] for t in cursor.fetchall())

        if self.table_name not in self.tables:
            self.conn.close()
            raise IOError("Table %s not found in %s" % (self.table_name, self.file_name))

        cursor = self.conn.execute('select * from %s limit 1' % self.table_name)
        self.columns = [d[0] for d in cursor.description]

    @property
    def cursor(self):
        logger_message(50, "The cursor attribute is deprecated, use self.conn instead")
        return self.conn

    def __iter__(self):
        for row in self.conn.execute('select * from %s' % self.table_name):
            yield ReportEntry(self.columns, row)

    def __call__(self, query):
        cursor = self.conn.execute(query)

        columns = tuple(d[0] for d in cursor.description)
        columns = tuple((c[1:-1] if c[0] == c[-1] == '"' else c) for c in columns)

        for row in cursor:
            yield ReportEntry(columns, row)

    def close(self):
        self.conn.close()


class SQLiteWriter(ReportWriter):
    '''The SQLite implementation of the mzReport writer class.'''
    def __init__(self, file_name, table_name=None, columns=None, default_columns=False):
        self.table_name = table_name or 'PeptideData'

        if default_columns:
            self.columns = default_columns + (columns or [])
        else:
            self.columns = columns

        if len(self.columns) > len(set(c.lower() for c in self.columns)):
            raise ValueError("Redundant columns: column headers must be unique")

        self.conn = sqlite3.connect(file_name)

        col_names = ",".join("'%s'%s" % (c, sqlite_types.get(default_types.get(c.lower()), ''))
                             for c in self.columns)

        table_command = "create table %s(%s)" % (self.table_name, col_names)
        try:
            self.conn.execute(table_command)
        except sqlite3.OperationalError:
            self.conn.close()
            raise IOError("Cannot create table: %s already exists" % self.table_name)

        if self.table_name == 'PeptideData':
            self.conn.execute("create table ImageData(RowID INTEGER, "
                              "Col TEXT, Tag TEXT, PlotData BLOB, "
                              "CONSTRAINT fk_peptide_image FOREIGN KEY (RowID) "
                              "REFERENCES PeptideData(_rowid_))")
            self.conn.execute("create unique index ImageDataIdx on ImageData(rowid,col)")

        self.lastID = 0

    @property
    def cursor(self):
        logger_message(50, "The cursor attribute is deprecated, just use self.conn instead")
        return self.conn

    def write(self, row, metadata=None):
        # error checking: want one value per column and nothing more
        # if row is a dict, keys should be in lower-case
        if len(row) < len(self.columns):
            raise ValueError('Must have values for each column')
        elif len(row) > len(self.columns):
            raise ValueError('Too many values')
        elif isinstance(row,dict):
            row = dict((k.lower(),v) for k,v in list(row.items()))
            if not all(k.lower() in row for k in self.columns):
                raise ValueError('Value dictionary does not match column headers')

        # supposedly this isn't a secure way to do a SQLite command, but I don't see
        # a better way to deal with an unknown number of values
        if isinstance(row,dict):
            self.lastID = self.conn.execute("INSERT OR REPLACE into %s values (%s)" % (self.table_name, ','.join(['?']*len(row))),
                                            tuple([row[col.lower()] for col in self.columns])).lastrowid
        else:
            self.lastID = self.conn.execute("INSERT OR REPLACE into %s values (%s)" % (self.table_name, ','.join(['?']*len(row))),
                                            tuple(row)).lastrowid

        if metadata:
            if self.table_name != 'PeptideData':
                raise ValueError('ImageData can only be added to the PeptideData table')

            # insert general images
            self.conn.executemany("INSERT OR REPLACE into ImageData values (?,?,?,?)",
                                  ((self.lastID,col,t,sqlite3.Binary(open(img,'rb').read()))
                                   for (col,t,img) in metadata if t.lower() == 'image'))

            # insert scan data for our images (to generate on the fly)
            self.conn.executemany("INSERT OR REPLACE into ImageData values (?,?,?,?)",
                                  ((self.lastID,col,t,scan)
                                   for (col,t,scan) in metadata if t.lower() in ('ms1','xic','ms2')))

        # commit changes every 1000 rows to reduce memory usage
        if self.lastID % 1000 == 0:
            self.conn.commit()

    def add_image(self, column, image):
        if self.table_name != 'PeptideData':
            raise ValueError('ImageData can only be added to the Data table')

        if self.lastID:
            self.conn.execute("INSERT into ImageData values (?,?,?,?)",
                              (self.lastID,column,'image',sqlite3.Binary(open(image,'rb').read())))
        else:
            raise IndexError("No row to add an image to")

    def close(self):
        self.conn.commit()
        self.conn.close()
