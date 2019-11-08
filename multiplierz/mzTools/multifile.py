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
import sqlite3

from collections import defaultdict

from multiplierz import mzReport
from multiplierz.mzReport import reader, writer
from functools import reduce


__all__ = ['filterJoin']

class CommaJoin(object):
    '''SQLite accumulator class to aggregate a comma-separated
    list of values. Removes duplicates and sorts the result.'''

    def __init__(self):
        self.vset = set()

    def step(self, value):
        self.vset.add(value)

    def finalize(self):
        return ','.join(sorted(self.vset))


def detect_matches(file_names, fields=[], tol_field=None, tolerance=0.0, save_file = ''):
    #mzTools.logger_message(30,'Detecting Matches...')

    #fields = fields or [mzReport.multiplierzHeaders[k]
                        #for k in ['acc','seq','var_mods']]

    detect_gen = _detect_matches(file_names, fields, tol_field, tolerance)

    cols = next(detect_gen)
    match_out = [cols]

    if save_file:
        writer = mzReport.writer(save_file, columns=cols)

    for line in detect_gen:
        match_out.append(line)
        if save_file:
            writer.write(line)

    detect_gen.close()

    if save_file:
        writer.close()

    #mzTools.logger_message(20,'Matches Detected')

    return match_out


def _detect_matches(file_names, fields, tol_field, tolerance=0.0):

    all_mzd = all((os.path.splitext(f)[1].lower() == '.mzd')
                  for f in file_names)

    # a separate set code for when all of the files are SQLite databases.
    # might as well take advantage of the ability to do a real query
    if all_mzd:
        # create a SQLite database in memory
        conn = sqlite3.connect(':memory:') # connection object

        # short-hand names for each file
        table_names = [('mzd%d' % i) for i,f in enumerate(file_names)]

        for f,t in zip(file_names, table_names):
            # attach each mzResults file to the database
            conn.execute('attach database (?) as (?)', (f, t))

        # add 'comma-join' function to aggregate file names
        conn.create_aggregate("cjoin", 1, CommaJoin)

        #if tol_field:
            #fields.append(tol_field)

        field_list = ','.join('"%s"' % f for f in fields)

        sub_query = ('select \'%s\' as "File Name",'
                     + field_list
                     + ' from %s.PeptideData')

        union_query = ' union '.join((sub_query % (t,t)) for t in table_names)

        query = ('select distinct cjoin("File Name"),%s'
                 + ' from (%s) group by %s') % (field_list, union_query, field_list)

        try:
            yield (fields
                   + [os.path.basename(f) for f in file_names]
                   + ["Detections"])

            # if tolerance is specified, need to do a second grouping step
            if tol_field:
                tol_dict = defaultdict(lambda : defaultdict(set))

                for row in conn.execute(query):
                    file_set = set(row[0].split(','))
                    t = tuple(row[1:-1])
                    tol_val = float(row[-1])

                    k = min(tol_dict[t] or [tol_val], key=lambda k: abs(k - tol_val))
                    if abs(k - tol_val) <= tolerance:
                        tol_dict[t][k].update(file_set)
                        tol_dict[t][tol_val].update(tol_dict[t][k])
                    else:
                        tol_dict[t][tol_val].update(file_set)

                for t in tol_dict:
                    for v in tol_dict[t]:
                        file_row = [int(f in tol_dict[t][v]) for f in table_names]
                        file_row.append(sum(file_row))

                        yield (t + (v,) + tuple(file_row))
            else:
                for row in conn.execute(query):
                    file_set = set(row[0].split(','))
                    file_row = [int(f in file_set) for f in table_names]
                    file_row.append(sum(file_row))

                    yield (row[1:] + tuple(file_row))
        finally:
            conn.close()
    else:
        if tol_field:
            rows = defaultdict(lambda : defaultdict(set))
        else:
            rows = defaultdict(set)

        for name in file_names:
            report = mzReport.reader(name)

            #mzTools.logger_message(10, name)

            for row in report:
                t = tuple(row.get(field.lower()) for field in fields)

                if tol_field:
                    tol_val = row.get(tol_field.lower())

                    rows[t][tol_val].add(name)

                    for k in list(rows[t].keys()):
                        if abs(k - tol_val) <= tolerance:
                            rows[t][tol_val].update(rows[t][k])
                            rows[t][k].add(name)
                else:
                    rows[t].add(name)

            report.close()

        if tol_field:
            yield (fields + [tol_field]
                   + [os.path.basename(f) for f in file_names]
                   + ["Detections"])

            for t in rows:
                for v in rows[t]:
                    file_row = [int(f in rows[t][v]) for f in file_names]
                    file_row.append(sum(file_row))

                    yield (t + (v,) + tuple(file_row))

        else:
            yield (fields
                   + [os.path.basename(f) for f in file_names]
                   + ["Detections"])

            for t in rows:
                file_row = [int(f in rows[t]) for f in file_names]
                file_row.append(sum(file_row))

                yield (t + tuple(file_row))


def filter_join(file_names, key_source_file, exclude=False,
                append=False, save_file_suffix='_filtered', ext='.xls'):
    if (not append) and (not save_file_suffix):
        raise ValueError('Save_file_suffix cannot be empty string when not combining output.')

    filtered_files = _filter_join(file_names, key_source_file, exclude, save_file_suffix)
    #if not filtered_files:
        #return

    cols = next(filtered_files)

    #Convert csv files to xls and append to manipulations log
    if append:
        combo_out_file = os.path.join(os.path.dirname(key_source_file),
                                      'Combined%s%s' % (save_file_suffix, ext))
        if os.path.exists(combo_out_file):
            os.remove(combo_out_file)

        # union of file columns--checking to see if files have the same columns
        #union_cols = reduce(set.union, [set(f[1]) for f in filtered_files])
        union_cols = reduce(set.union, [set(c) for c in cols])
        # if the union is the size of the smallest set, they are all equal
        #if len(union_cols) == min(len(f[1]) for f in filtered_files):
        if len(union_cols) == min(len(c) for c in cols):
            same_cols = True
            if 'File' not in union_cols:
                #newcols = ['File'] + filtered_files[0][1] # they're all the same so use the first one
                newcols = ['File'] + list(cols[0]) # they're all the same so use the first one
            else:
                #newcols = filtered_files[0][1]
                newcols = list(cols[0])
        else:
            # files have different columns
            same_cols = False
            if 'File' not in union_cols:
                # use the largest set of columns as initial template
                #newcols = ['File'] + max((f[1] for f in filtered_files), key=len)
                newcols = ['File'] + list(max(cols, key=len))
                # stick the remaining columns on the end, sorted as strings
                newcols += sorted(union_cols.difference(newcols))
            else:
                # use the largest set of columns as initial template
                #newcols = max((f[1] for f in filtered_files), key=len)
                newcols = list(max(cols, key=len))
                # stick the remaining columns on the end, sorted as strings
                newcols += sorted(union_cols.difference(newcols))

        if 'Filter Key' not in newcols:
            newcols.insert(1, 'Filter Key')

        combo_report = mzReport.writer(combo_out_file, columns=newcols)

        # if the columns are all the same, this is straightforward, just write them all out
        if same_cols:
            for (filename,columns,filedata) in filtered_files:
                for row in filedata:
                    row['File'] = filename # will either add File or overwrite it
                    combo_report.write(row)
        # if not, it's slightly more complicated--create dictionary for each row with None
        # if a value is missing
        else:
            blank_row = dict((c.lower(),None) for c in newcols)
            for (filename,columns,filedata) in filtered_files:
                for row in filedata:
                    new_row = blank_row.copy()
                    new_row.update(row)
                    new_row['file'] = filename # will either add File or overwrite it
                    combo_report.write(new_row)

        combo_report.close()
    else:
        for (filename,columns,filedata) in filtered_files:
            # write each sheet to a separate file
            if 'Filter Key' not in columns:
                columns = ['Filter Key'] + columns

            rep = mzReport.writer(save_file_suffix.join(os.path.splitext(filename)),
                                  columns=columns)
            for row in filedata:
                rep.write(row)

            rep.close()


def _filter_join(file_names, key_source_file, exclude, save_file_suffix = '_filtered'):

    all_mzd = (os.path.splitext(key_source_file)[1].lower() == '.mzd'
               and all((os.path.splitext(f)[1].lower() == '.mzd')
                       for f in file_names))

    if all_mzd:
        # create a SQLite database in memory
        conn = sqlite3.connect(':memory:') # connection object
        conn.execute('attach database (?) as key_table', (key_source_file,))

        cols = [d[0] for d in conn.execute('select * from key_table.PeptideData limit 1').description]
        col_set = set(cols)

        # short-hand names for each file
        table_names = [('mzd%d' % i) for i,f in enumerate(file_names)]

        col_list = []
        for f,t in zip(file_names, table_names):
            # attach each mzResults file to the database
            conn.execute('attach database (?) as (?)', (f, t))
            col_list.append(tuple(d[0] for d in conn.execute('select * from %s.PeptideData limit 1' % t).description))

        yield col_list

        try:
            for f,t,t_cols in zip(file_names, table_names, col_list):
                #mzTools.logger_message(20, f)

                t_set = set(t_cols)
                u_str = ' and '.join('(A."%s" = B."%s" or ifnull(A."%s", B."%s") is NULL)' % (c,c,c,c)
                                     for c in t_cols if c in col_set)

                query = 'select A.* from %s.PeptideData as A, key_table.PeptideData as B where %s' % (t, u_str)

                if exclude:
                    query = 'select * from %s.PeptideData except %s' % (t, query)

                cur = conn.execute(query)

                res_cols = [d[0] for d in cur.description]
                res_dict = dict((c,i) for i,c in enumerate(res_cols))

                yield (f, t_cols,
                       (dict([('Filter Key', '|'.join(str(row[res_dict[c]])
                                                      for c in cols if c in t_set))]
                             + list(zip(t_cols,(v for i,v in enumerate(row) if res_cols[i] in t_cols))))
                        for row in cur))

                conn.execute('detach database %s' % t)
        finally:
            conn.close()
    else:
        source_report = mzReport.reader(key_source_file)

        col_list = []
        for name in file_names:
            rdr = mzReport.reader(name)
            col_list.append(tuple(rdr.columns))
            rdr.close()

        yield col_list

        try:
            for name in file_names:
                #mzTools.logger_message(20, name)
                this_rep = mzReport.reader(name)

                key_cols = [col for col in source_report.columns if col in this_rep.columns]
                filter_keys = set()

                for row in source_report:
                    if int(row['Detections']) > 1:
                        filter_keys.add("|".join(str(row[col]) for col in key_cols))

                if exclude:
                    filtered_data = []
                    for row in this_rep:
                        filter_key = "|".join(str(row[col]) for col in key_cols)
                        if filter_key not in filter_keys:
                            row['Filter Key'] = filter_key
                            filtered_data.append(row)
                else:
                    filtered_data = []
                    for row in this_rep:
                        filter_key = "|".join(str(row[col]) for col in key_cols)
                        if filter_key in filter_keys:
                            row['Filter Key'] = filter_key
                            filtered_data.append(row)

                yield (name, this_rep.columns, filtered_data)

                this_rep.close()
        finally:
            this_rep.close()






# Rewrite of all of the above, so I'll know what it actually does.

def filterJoin(filenames, matchColumns, returnMode, outputKeyFile,
               combinedOutputFile = None, outputFileType = '.xlsx',
               outputTag = None,
               tolerance = None, toleranceColumn = None):
    """
    Produces a joined file, filtering out either repeat or unique
    PSMS.
    
    If returnMode is 'matched', output file contains one instance
    of each PSM group;
    if returnMode is 'unmatched', the output file contains every
    PSM that wasn't part of a larger PSM group.
    (Where 'PSM group' is a set of PSMs that are identical based on
    the given matchColumns + toleranceColumn.)
    If 'both', both kinds of output are produced.  
    
    Returns output file name(s); a tuple in the case of 'both.'
    """
    
    assert returnMode in ['matched', 'unmatched', 'both']
    #if not outputFileBase:
        #outputFileBase = filenames[0]
    
    data = []
    columnLists = []
    for filename in filenames:
        subdata = []
        inputfile = reader(filename)
        columnLists.append(inputfile.columns)
        for psm in inputfile:
            psm['Source'] = filename
            subdata.append(psm)
        data.append(subdata)
        inputfile.close()
    
    assert all([columnLists[0] == x for x in columnLists]), "Heterogeneous data columns!"
    
    datadict = defaultdict(list)
    for subdata in data:
        for psm in subdata:
            signature = tuple([psm[x] for x in matchColumns])
            datadict[signature].append(psm)
    
    if toleranceColumn:
        toldatadict = {}
        for signature, sigGroup in list(datadict.items()):
            subGroups = []
            for psm in sigGroup:
                match = False
                for subGroup in subGroups:
                    if all([abs(psm[toleranceColumn] - subpsm[toleranceColumn]) < tolerance
                            for subpsm in subGroup]):
                        match = True
                        subGroup.append(psm)
                        break
                    
                if not match:
                    subGroups.append([psm])
            
            for index, subGroup in enumerate(subGroups):
                subSig = tuple(list(signature) + index)
                toldatadict[subSig] = subGroup
        
        datadict = toldatadict
        
    if outputKeyFile:
        keyfile = writer(outputKeyFile,
                         columns = ['PSM Key'] + filenames)
        for signature, psmGroup in list(datadict.items()):
            line = {}
            line['PSM Key'] = '|'.join([str(x) for x in signature])
            line.update([(x, len([y for y in psmGroup if y['Source'] == x])) for x in filenames])
            keyfile.write(line)
        
        keyfile.close()
        
        
    outputpsms = []
    if returnMode == 'matched' or returnMode == 'both':
        #outputFileName = outputFileBase + '_matchedPSMs' + outputFileType
        #outputfile = writer(outputFileName, columns = ['Source'] + columnLists[0])
        
        for psmGroup in list(datadict.values()):
            if len(psmGroup) > 1:
                exemplar = psmGroup[0]
                sourceFiles = '; '.join(set([x['Source'] for x in psmGroup]))
                exemplar['source'] = sourceFiles
                #outputfile.write(exemplar)
                outputpsms.append(exemplar)
        #outputfile.close()
    
    if returnMode == 'unmatched' or returnMode == 'both':
        #outputFileName = outputFileBase + '_uniquePSMs' + outputFileType
        #outputfile = writer(outputFileName, columns = ['Source'] + columnLists[0])
        
        for psmGroup in list(datadict.values()):
            if len(psmGroup) == 1:
                #outputfile.write(psmGroup[0])
                outputpsms.append(psmGroup[0])
            
        #outputfile.close()
    
    outputs = []
    if outputTag:
        outputfiles = [(x, '.'.join(x.split('.')[:-1] + [outputTag, outputFileType])) for x in filenames]
        for filename, outputfile in outputfiles:
            output = writer(outputfile, columns = ['Source'] + columnLists[0])
            for psm in [x for x in outputpsms if x['Source'] == filename]:
                output.write(psm)
            output.close()
    
        outputs = [x[1] for x in outputfiles]
        
    if combinedOutputFile:
        output = writer(combinedOutputFile, columns = ['Source'] + columnLists[0])
        for psm in outputpsms:
            output.write(psm)
        output.close()
        outputs.append(combinedOutputFile)
        
    return outputs
    #if returnMode == 'matched':
        #return outputFileBase + '_matchedPSMs' + outputFileType
    #elif returnMode == 'unmatched':
        #return outputFileBase + '_uniquePSMs' + outputFileType
    #elif returnMode == 'both':
        #return (outputFileBase + '_matchedPSMs' + outputFileType, outputFileBase + '_uniquePSMs' + outputFileType)
        