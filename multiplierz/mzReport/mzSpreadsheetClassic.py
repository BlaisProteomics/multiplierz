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

import os

import win32com.client
from pythoncom import CoInitialize, CoUninitialize

#import xlrd

from multiplierz.mzReport import ReportReader, ReportWriter, ReportEntry, default_columns

win32com.client.gencache.is_readonly=False
xlApp = None

def _start_excel():
    """Launches Excel. There is no need to call this directly,
    the Excel functions call it when necessary.

    Example:
    >>> my_excel_app = _start_excel()
    >>> my_excel_app
    <COMObject Excel.Application>

    """

    global xlApp
    CoInitialize()
    if xlApp is None:
        xlApp = win32com.client.DispatchEx('Excel.Application')
        xlApp.Visible = 0
    return xlApp


def _stop_excel():
    """Quits Excel if open

    Example:
    >>> _stop_excel()

    """

    global xlApp
    if xlApp is not None:
        xlApp.Visible = 0
        xlApp.Quit()
        xlApp = None
    CoUninitialize()


def get_xl_sheet(file_name, sheet_name='Data'):
    if os.path.exists(file_name):
        xl_app = _start_excel()
        try:
            workbook = xl_app.Workbooks.Open(file_name)

            sheet_set = set(str(workbook.Worksheets(i).Name)
                             for i in range(1,workbook.Worksheets.Count + 1))

            if sheet_name not in sheet_set:
                raise IOError('File %s has no sheet %s' % (file_name, sheet_name))

            # get all data and store it in a list of lists
            data = workbook.Worksheets(sheet_name).Range("A1").CurrentRegion.Value
            if data:
                data = [list(row) for row in data]
                #metadata = []
                #for i in range(1,len(data[0])+1):
                    #if (workbook.Worksheets(sheet_name).Cells(2,i).HasFormula
                        #and workbook.Worksheets(sheet_name).Cells(2,i).FormulaR1C1[0] == '='):
                        #metadata.extend((j,i,'formula',
                                         #workbook.Worksheets(sheet_name).Cells(j,i).FormulaR1C1)
                                        #for j in range(2,len(data)+1))
            else:
                data = [[]]
                #metadata = []

            workbook.Close(SaveChanges=0)
        except Exception, e:
            raise
        finally:
            _stop_excel()

        return data #, metadata
    else:
        raise IOError('No such file: %s' % file_name)


def write_xl_sheet(file_name, sheet_name, data, metadata,
                   overwrite=False, before=None, after=None):
    '''Writes a single Excel sheet to disk.'''

    xl_app = _start_excel()

    try:
        if os.path.exists(file_name):
            workbook = xl_app.Workbooks.Open(file_name)
            new_file = False
        else:
            workbook = xl_app.Workbooks.Add()
            if os.path.splitext(file_name)[1] == '.xlsx':
                workbook.SaveAs(file_name, 51)
            else:
                workbook.SaveAs(file_name, -4143)
            new_file = True

        sheet_list = [str(workbook.Worksheets(i).Name)
                      for i in range(1,workbook.Worksheets.Count + 1)]

        if sheet_name not in sheet_list:
            workbook.Worksheets.Add().Name = sheet_name
        elif not overwrite:
            raise IOError('Sheet %s already exists' % sheet_name)
        else:
            # need to test if this will clear out comments
            workbook.Worksheets(sheet_name).Range("A1").CurrentRegion.Clear()

        # delete all other sheets if this is a new file (ie Sheet1, Sheet2, Sheet3)
        if new_file:
            for s in sheet_list:
                if s != sheet_name:
                    workbook.Worksheets(s).Delete()

        if len(data) and len(data[0]):
            address_end = workbook.Worksheets(sheet_name).Cells(len(data),
                                                                len(data[0])).Address
            address = "%s:%s" % ("A1", address_end)
            workbook.Worksheets(sheet_name).Range(address).Value = data

        if file_name.lower().endswith('.csv'):
            workbook.SaveAs(file_name,0x6)
            workbook.Close(SaveChanges=1)
            return

        for i,(x,y,t,md) in enumerate(metadata):
            if i % 100 == 0:
                workbook.Save()

            if isinstance(t, tuple) and t[0] == 'image':
                height, width = t[1:]
                workbook.Worksheets(sheet_name).Cells(x,y).ClearComments()
                workbook.Worksheets(sheet_name).Cells(x,y).AddComment(" ")
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Fill.UserPicture(md)                
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Height = height * 72/100
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Width = width * 72/100                
            if t == 'image':
                #import Image

                #wdth,hgt = Image.open(md).size



                workbook.Worksheets(sheet_name).Cells(x,y).ClearComments()
                workbook.Worksheets(sheet_name).Cells(x,y).AddComment(" ")
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Fill.UserPicture(md)
                # ugly solution but seems to work: 72% scaling to make the image sizes accurate
                #workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Height = hgt * 72/100
                #workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Width = wdth * 72/100
            
                #from PIL.Image import open
                #wdth, hgt = open(md).size                    
                #workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Height = hgt * 72/100
                #workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Width = wdth * 72/100                    
            elif t == 'prot coverage':
                workbook.Worksheets(sheet_name).Cells(x,y).ClearComments()
                workbook.Worksheets(sheet_name).Cells(x,y).AddComment(md[0])
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters().Font.Name = "courier new"
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters().Font.Size = 12
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters().Font.Bold = False
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.Fill.ForeColor.SchemeColor = 9
                for s in range(1, len(md[0]), 61):
                    workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters(s,5).Font.Bold = True
                for (s,e) in md[1]:
                    workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters(s,e-s+1).Font.Bold = True
                    workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters(s,e-s+1).Font.ColorIndex = 3
                for (s,e) in md[2]:
                    workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.Characters(s,e-s+1).Font.Underline = True
                workbook.Worksheets(sheet_name).Cells(x,y).Comment.Shape.TextFrame.AutoSize = True
            elif t == 'formula':
                workbook.Worksheets(sheet_name).Cells(x,y).FormulaR1C1 = md
            elif t == 'text':
                workbook.Worksheets(sheet_name).Cells(x,y).NumberFormat = u'@'
                try:
                    if len(md) > 32000:
                        md = "!!!Truncated!!!" + md[:32000]
                    workbook.Worksheets(sheet_name).Cells(x,y).Value = md
                except:
                    print x
                    print y
                    print t
                    print md
                    print i
                    raise ValueError("Error!")

        if before and before in sheet_list:
            before = sheet_list.index(before) + 1
        else:
            before = None

        if after and after in sheet_list:
            after = sheet_list.index(after) + 1
        else:
            after = None

        if before or after:
            workbook.Worksheets(sheet_name).Move(Before=before, After=after)

        workbook.Activate()

        workbook.Close(SaveChanges=1)
    except Exception, e:
        raise
    finally:
        _stop_excel()


def genbank_sheet_format(file_name, num_rows):
    '''Reads the Data sheet from an Excel report and retrieves GenBank information
    for the identified proteins. Assumes that 'Accession Number' is a gi number.

    file_name - the name of the report to be processed--must have a sheet called
                'GenBank_Info'
    num_rows - the number of (non-header) rows that the sheet contains. Passing
               this in is easier than trying to figure it out via COM.
    '''

    xl_app = _start_excel()

    try:
        workbook = xl_app.Workbooks.Open(file_name)

        sheet_list = [str(workbook.Worksheets(i).Name)
                      for i in range(1,workbook.Worksheets.Count + 1)]

        if 'GenBank_Info' not in sheet_list:
            raise IOError('File %s has no sheet GenBank_Info' % file_name)

        NUM_COLS = 13

        for i in range(NUM_COLS):
            workbook.Worksheets('GenBank_Info').Cells(1, i+1).Font.Bold = True

        colorIndex = 15

        for i in range(1,14,2):
            workbook.Worksheets('GenBank_Info').Cells(1,i).Interior.ColorIndex = colorIndex

        #Set Column Widths
        workbook.Worksheets('GenBank_Info').Columns.AutoFit()

        #Set Border
        workbook.Worksheets('GenBank_Info').Range("A1:M1").Borders.Weight = 3
        workbook.Worksheets('GenBank_Info').Range("C2").Select()
        workbook.Application.ActiveWindow.FreezePanes = True

        for i in range(2, num_rows+2):
            for j in range(1, NUM_COLS+1, 2):
                workbook.Worksheets('GenBank_Info').Cells(i,j).Interior.ColorIndex = colorIndex

        # set column widths
        workbook.Worksheets('GenBank_Info').Range("A1:M1").Columns.AutoFit()
        workbook.Worksheets('GenBank_Info').Range("A1").Select()

        # move sheet
        workbook.Worksheets('GenBank_Info').Move(Before=None, After=workbook.Worksheets(workbook.Worksheets.Count))
        workbook.Close(SaveChanges=1)
    except Exception, e:
        raise
    finally:
        _stop_excel()


class XLSheetReader(ReportReader):
    def __init__(self, file_name, sheet_name=None):
        if not os.path.exists(file_name):
            raise IOError("No such file: '%s'" % os.path.basename(file_name))
        self.file_name = file_name
        self.sheet_name = sheet_name or 'Data'

        self._data = get_xl_sheet(self.file_name, self.sheet_name)
        self.columns = self._data[0]

        # the reader can't read metadata--that's why we lose images when
        # we modify excel sheets

    def __iter__(self):
        '''Returns the rows of data (not the headers) for reading'''
        for i,row in enumerate(self._data):
            if i == 0:
                continue
            yield ReportEntry(columns=self.columns, values=row)

    def close(self):
        '''Doesn't do anything, because the sheet is loaded into memory
        at initialization'''
        self._data = [[]]


class XLSheetWriter(ReportWriter):
    def __init__(self, file_name, sheet_name=None, columns=None, default_columns=False):
        '''Initializes the writer. For Excel sheets, this doesn't include opening
        a file handle--all IO is done when the file is closed. This means that you
        *must* call XLSheetWriter.close() or your file won't be saved!
        '''
        
        self.file_name = file_name.replace('/', '\\')
        self.sheet_name = sheet_name or 'Data'

        self.offset = 0

        if default_columns:
            self.columns = default_columns + (columns or [])
        else:
            self.columns = columns[:]

        self._data = [self.columns]
        self._metadata = []

        if len(self.columns) > len(set(c.lower() for c in self.columns)):
            raise ValueError("Redundant columns: column headers must be unique")

    def write(self, row, metadata=None):
        # error checking: want one value per column and nothing more
        if len(row) < len(self.columns):
            raise ValueError('Must have values for each column')
        elif len(row) > len(self.columns):
            raise ValueError('Too many values')
        elif isinstance(row,dict):
            row = dict((k.lower(),v) for k,v in row.items())
            if not all(k.lower() in row for k in self.columns):
                raise ValueError('Value dictionary does not match column headers')

        if isinstance(row,dict):
            row = [row[key.lower()] for key in self.columns]

        index = len(self._data) + 1
        self._metadata.extend((index,i+1,'text',d) for i,d in enumerate(row)
                              if isinstance(d, (str,unicode)) and len(d) > 768)

        self._data.append([(d[:768] if isinstance(d, (str,unicode)) else d) for d in row])

        if metadata:
            self._metadata.extend((index,self.columns.index(c)+1,t,md) for (c,t,md) in metadata)

    def add_image(self, column, image):
        index = len(self._data) + 1
        if index:
            self._metadata.append((index, self.columns.index(column) + 1, 'image', image))
        else:
            raise IndexError, "No row to add an image to"

    def close(self, *args, **kwargs):
        '''Close the file. In this case, we really open, write, and close all in
        this step. Additional arguments can be passed to the xlfile's write method.'''

        write_xl_sheet(self.file_name, self.sheet_name,
                       self._data, self._metadata,
                       overwrite=True, **kwargs)

        # zero out the data to clean up memory
        self._data = [[]]
        self._metadata = []
