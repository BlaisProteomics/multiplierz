
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

# This module is used to interact with the mascot web interface and parse mascot files

#from __future__ import division
#from __future__ import print_function
#from past.builtins import cmp
#from future import standard_library
#standard_library.install_aliases()
#from builtins import zip
#from builtins import chr
##from builtins import str
#from builtins import hex
#from builtins import range
#from past.utils import old_div
#from builtins import object
import io
import csv
#import html.entities
import os
import re
import tempfile
import time
import webbrowser
import sys

#import urllib.request, urllib.error, urllib.parse
from io import StringIO

httpPackage = None


import requests

from multiplierz.internalAlgorithms import collectByCriterion
    
#import multiplierz.msparser as ms

from collections import defaultdict
#from urllib.parse import urlencode, unquote
import multiplierz.mass_biochem as mzFunctions
import multiplierz.mzReport as mzReport
#from multiplierz.mascot.report import MascotReport

from multiplierz import myData, logger_message
from multiplierz.settings import settings

#if sys.version_info[0] == 3:
    #import multiplierz.mzSearch.mascot.msparser_py3.msparser as ms
#elif sys.version_info[0] == 2:
    #import msparser_py2.msparser as ms
#else:
    #raise Exception(sys.version_info)


default_max_hits = 999999
default_ion_cutoff = 8
default_show_same_set = True
default_show_sub_set = True
default_bold_red = False


multiplierzHeaderConversions = dict((k,mzReport.multiplierzHeaders.get(v,v)) for k,v in
                                    zip(('prot_acc', 'prot_desc', 'prot_score', 'prot_mass',
                                         'prot_matches', 'prot_hit_num', 'pep_query', 'pep_exp_mz',
                                         'pep_exp_mr', 'pep_exp_z', 'pep_calc_mr', 'pep_delta',
                                         'pep_start', 'pep_end', 'pep_miss', 'pep_score',
                                         'pep_expect', 'pep_rank', 'pep_res_before', 'pep_seq',
                                         'pep_res_after', 'pep_var_mod', 'pep_scan_title', 'qualifiers',
                                         'TotalIonsIntensity', 'NumVals', 'StringIons1', 'StringIons2',
                                         'StringIons3', 'pep_scan_title'),
                                        ('acc', 'prot_desc', 'prot_score', 'prot_mass', 'prot_matches',
                                         'prot_rank', 'pep_query', 'mz', 'mr', 'charge', 'pred_mr',
                                         'delta', 'start', 'end', 'miss', 'pep_score', 'Peptide EValue',
                                         'pep_rank', 'res_before', 'seq', 'res_after', 'var_mods',
                                         'spec_desc', 'Qualifiers', 'Total Ion Intensity', 'Ion Count',
                                         'ion_str', 'Ion String 2', 'Ion String 3', 'spec_desc')))


filetypes = {"MascotDAT" : 'dat',
             'CSV' : 'csv',
             'MZIdentML' : 'mzid'}


float_pattern = re.compile('-?[0-9]*\.[0-9]*$') # Caution- also matches .
int_pattern = re.compile('-?[0-9]+$')
exp_pattern = re.compile('-?[0-9]+\.?[0-9]*[Ee]-?[0-9]+$')

def make_type(value):
    if float_pattern.match(value):
        return float(value)
    elif exp_pattern.match(value):
        return float(value)
    elif int_pattern.match(value):
        return int(value)
    else:
        return value


# Split by commas, except within quotes.
def escapable_split(line, split, expected = None):
    splits = re.findall(r'[^"%s]*%s|"[^"]*"%s|"[^"]*"\Z' %
                        (split, split, split), line.strip())
    splits = list(map(lambda x: x.strip(',').strip('"'), splits))
    if expected is not None:
        assert len(splits) == expected
    return splits

def MascotCSVToMultiplierzReport(csv_data, outputfile = None):        
    try:
        csv_data = open(csv_data, 'r')
    except TypeError: # Probably, was passed a file pointer already.
        pass
    
    
    # There's sensitivity to the line break style here- I don't know enough
    # about string processing in Python 3 to know what the proper way to deal
    # with it is.
    header = []
    for line in csv_data:
        line = line.strip('\n')
        if outputfile is None and 'Peak list data path' in line:
            outputfile = line.split(',')[1].strip().strip('"') + '.xlsx'
        if 'Database,' in line:
            db_str = line.strip().split(',')[1].strip('"')
            db_lookup = dict([(x.split('::')[0], x) for x in db_str.split()])
        if 'Protein hits' not in line:
            header.append(line)
        else:
            header.append(line)
            break
    
    varmod_numbering = {}
    varmod_start = [i for i, h in enumerate(header) if 'Variable modifications' in h][0]
    varmod_stop = [i for i, h in enumerate(header) if 'Search Parameters' in h][0]
    for line in header[varmod_start+1:varmod_stop]:
        if not line.strip() or 'Identifier' in line: continue
        words = line.split(',')
        varmod_numbering[int(words[0])] = words[1].strip('"').split(' (')[0] # This is a bit hacky.
    
    for line in csv_data:
        if line.strip():
            break
    assert 'pep_seq' in line and 'pep_exp_mz' in line
    columns = line.strip().split(',')
    indexes = {col:i for i, col in enumerate(columns)}
    
    out_columns = set(mzReport.default_columns)
    
    rows = []
    for line in csv_data:
        values = escapable_split(line, ',', len(columns))
        
        # Most fields have direct equivalents in the Mascot report.
        row = {mz_name:make_type(values[indexes[csv_name]].strip('"')) for
               csv_name, mz_name in multiplierzHeaderConversions.items()
               if mz_name in mzReport.default_columns}
        if len(db_lookup) == 1:
            row['Protein Database'] = list(db_lookup.values())[0]
        else:
            row['Protein Database'] = db_lookup[row['Accession Number'].split('::')[0]]
            row['Accession Number'] = row['Accession Number'].split('::')[1]
            
        mods = values[indexes['pep_var_mod']]
        mod_loc_str = values[indexes['pep_var_mod_pos']]
        if mods or mod_loc_str:
            # A bunch of fiddling around has to be done to convert the obscure
            # native-Mascot modification representation to the pretty MultipilerZ
            # representation.
    
            mod_locs = [(i, x) for i, x in 
                        enumerate(mod_loc_str.split('.')[1], start = 1)
                        if x != '0']
            mod_strs = []
            for pos, modnum in mod_locs:
                modname = varmod_numbering[int(modnum)]
                aa = row['Peptide Sequence'][pos-1]
                mod_strs.append("%s%d: %s" % (aa, pos, modname))
    
            if mod_loc_str.split('.')[0] != '0':
                modname = varmod_numbering[int(lod_locs.split('.')[0])]
                mod_strs.append("Lterm: %s" % modname)
            if mod_loc_str.split('.')[-1] != '0':
                modname = varmod_numbering[int(lod_locs.split('.')[0])]
                mod_strs.append("Cterm: %s" % modname)        
            
            row['Variable Modifications'] = '; '.join(mod_strs)
        else:
            row['Variable Modifications'] = ''
        
        rows.append(row)
        assert len(row) == len(mzReport.default_columns)
    
    output = mzReport.writer(outputfile, columns = mzReport.default_columns)
    for row in rows:
        output.write(row)
    output.close()
    
    if outputfile.lower().endswith('xlsx'):
        output = mzReport.writer(outputfile, sheet_name = 'Mascot_Header',
                        columns = ["", '', '', '', ''])
        for line in header:
            vals = line.strip().split(',')
            vals = [x.strip('"') for x in vals]
            if len(vals) < 5:
                vals = vals + ['']*(5 - len(vals))
            vals = list(map(make_type, vals))
            output.write(vals)
        output.close()
    
    return outputfile









class MascotReport:
    def __init__(self, server=None, version=2.6, login=None, password=None, cleanup=True):
        if server:
            self.mascot = mascot(server, version)
            self.login = login
            self.password = password
            self.cleanup = cleanup # Set to False to leave dat file in place for further queries
            self.mascot.login(login, password)
        else:
            self.mascot = mascot()
            self.cleanup = cleanup
            print("Note: No server specified, Mascot-independent mode.")

    def mascot_header(self, report_file, mascot_header):
        '''Adds a Mascot Header page to a report (XLS or MZD)'''

        logger_message(30, 'Adding Mascot Header...')

        if report_file.lower().endswith('.mzd'):
            mascot_rep = mzReport.writer(report_file, columns=mascot_header[0], table_name='MascotHeader')
        else:
            mascot_rep = mzReport.writer(report_file, columns=mascot_header[0], sheet_name='Mascot_Header')

        for line in mascot_header[1:]:
            mascot_rep.write(line)

        mascot_rep.close()

        logger_message(20, 'Mascot Header Complete')

    def mascot_headers(self, report_file, mascot_headers):
        '''Combines a list of Mascot header pages into a single page (XLS or MZD)'''

        logger_message(30, 'Adding Mascot Headers...')

        first_columns = set()
        columns = defaultdict(dict)

        report_files = [f for f,mh in mascot_headers]

        for f,mascot_header in mascot_headers:
            first_col = []

            for line in mascot_header[1:]:
                if line[0] == ' ' or (isinstance(line[1], (str,str))
                                      and line[1].startswith('-----')):
                    continue
                columns[f][line[0]] = line[1]
                first_col.append(line[0])

            first_columns.add(tuple(first_col))

        if len(first_columns) > 1:
            logger_message(20, 'Headers seems different, will try to merge them')

            main_h = list(max(first_columns, key=len))
            main_h += sorted(reduce(set.union, first_columns, set()).difference(main_h))
        else:
            main_h = list(first_columns.pop())

        cols = ['Header'] + [(os.path.basename(f) if f else ('-' * 50))
                             for f in report_files]

        if report_file.lower().endswith('.mzd'):
            mascot_rep = mzReport.writer(report_file, columns=cols,
                                         table_name='MascotHeader')
        else:
            mascot_rep = mzReport.writer(report_file, columns=cols,
                                         sheet_name='Mascot_Header')

        for col in main_h:
            row = [col]

            for f in report_files:
                if col in columns[f]:
                    row.append(columns[f][col])
                else:
                    row.append(None)

            mascot_rep.write(row)

        mascot_rep.close()

    def retrieve_report_data(self, mascot_id, report_columns,
                             date = None,
                             rank_one = False,
                             ext = '.xlsx',
                             local_dat_file = None,
                             retain_dat_file = False,
                             **mascot_options):
        """
        Retrieves a .DAT file (or opens one local) and reads
        into tabular format, handling several issues such as
        decoy mode data, header data, and version-dependent
        rows.
        """
        
        if local_dat_file:
            datfile = os.path.abspath(local_dat_file)
        else:
            datfile = self.mascot.download_dat(myData, mascot_id, date)
            
        dat_interface = MascotDatFile(datfile, **mascot_options)
        
        header_sheet = dat_interface.mascot_header()
        data_sheet = []
        for values in dat_interface.peptide_report():
            row = mzReport.ReportEntry(mzReport.default_columns,
                                       values)
            
            if rank_one and row['Peptide Rank'] != 1:
                continue
            
            data_sheet.append(row)
        dat_interface.close()        
        
        # With 2.5 Mascot, the decoy mode test crashes Python outright.
        #if dat_interface.hasDecoyHits():
        try:
            decoy_dat_interface = MascotDatFile(datfile, 
                                                          decoyMode = True,
                                                          **mascot_options)
            for values in decoy_dat_interface.peptide_report():
                row = mzReport.ReportEntry(report_columns, values)
                data_sheet.append(row)
            decoy_dat_interface.close()       
        except Exception as err:
            print(("Reading decoy data failed with error: %s" % err))
        
        if not (retain_dat_file or local_dat_file):
            os.remove(datfile)
        
        return header_sheet, data_sheet



    def get_reports(self, mascot_ids, 
                    dates = None,
                    outputfiles = [],
                    ext = 'xlsx',
                    chosen_folder = '',
                    local_dat_files = None,
                    **report_kwargs):
        # Since msparser isn't working in Python 3.7, we can instead just convert
        # the Mascot-native CSV output to just about the familiar multiplierz report.
        
        assert not local_dat_files, "Local dat file loading won't work until msparser-based reporting is restored."
        if outputfiles:
            assert len(outputfiles) == len(mascot_ids), "Output name list length does not match mascot id list."
        else:
            outputfiles = [None] * len(mascot_ids)
        
        if dates:
            assert len(dates) == len(mascot_ids), "%d dates given for %d searches: %s" % (len(dates), 
                                                                                          len(mascot_ids), 
                                                                                          dates)
        else:
            dates = [None] * len(mascot_ids)
        
        outputs = []
        for f0_num, date, outputfile in zip(mascot_ids, dates, outputfiles):
            csv_data = self.mascot.download_file(f0_num, date = date,
                                                 filetype = 'CSV',
                                                 buffer_results = False,
                                                 **report_kwargs)
            output = MascotCSVToMultiplierzReport(csv_data,
                                                  outputfile)
            if outputfile is not None:
                assert output == outputfile
            outputs.append(output)
        
        return outputs
            
            
    
    
    def get_reports_NORMAL(self, mascot_ids,
                    dates = None,
                    outputfile = None,
                    ext = None,
                    chosen_folder = '',
                    local_dat_files = [],
                    **report_kwargs):
        

        if ext:
            ext = ext.lstrip('.')
            
        report_columns = mzReport.default_columns
        if float(self.mascot.version) >= 2.3 and 'Protein Database' not in report_columns:
            report_columns.insert(1, 'Protein Database')        
            
        if dates:
            assert len(dates) == len(mascot_ids), "Mismatched date list provided."
            mascot_searches = list(zip(mascot_ids, dates))
        else:
            mascot_searches = [(x, None) for x in mascot_ids]
            
        reports = []
        for mascot_id, date in mascot_searches:
            header, psms = self.retrieve_report_data(mascot_id, report_columns,
                                                     date,
                                                     **report_kwargs)
            datafilename = header[7][1] or mascot_id
            reports.append((mascot_id, datafilename, header, psms))
        for dat_file in local_dat_files:
            header, psms = self.retrieve_report_data(None, report_columns,
                                                     local_dat_file = dat_file,
                                                     **report_kwargs)
            datafilename = header[7][1] or mascot_id
            reports.append((None, datafilename, header, psms))     
            
        imputed_output_file_name = False
        if outputfile and ext:
            if not outputfile.lower().endswith(ext):
                outputfile += '.' + ext        
        elif not outputfile:
            if not ext:
                ext = 'xlsx'
            if reports:
                outputfile = '_'.join(zip(*reports)[1]) + '.' + ext
            elif mascot_ids:
                outputfile = '_'.join(mascot_ids) + '.' + ext.strip('.')            
            imputed_output_file_name = True
            #outputfile = '.'.join([reports[0][1], ext.strip('.')])
        elif outputfile and not ext:
            ext = outputfile.split('.')[-1]
            assert ext in ['csv', 'xlsx', 'xls', 'mzd', 'mzid']

          
        assert outputfile 
        if chosen_folder and not os.path.isabs(outputfile):
            outputfile = os.path.join(chosen_folder, outputfile)
        
        if ext and 'mzid' in ext:
            assert len(mascot_ids) == 1, ("Combined result file not supported "
                                          "for mzIdentML files.")
            self.mascot.download_mzid(mascot_id,
                                      save_file = outputfile,
                                      date = date)
        elif len(mascot_ids) <= 1:
            mascot_id, datafilename, header, psms = reports[0]
            if imputed_output_file_name or not outputfile:
                outputfile = datafilename + '.' + ext
                if chosen_folder:
                    outputfile = os.path.join(chosen_folder, outputfile)
            output = mzReport.writer(outputfile,
                                     columns = header[0],
                                     sheet_name = 'Mascot_Header')
            for line in header[1:]:
                output.write(line)
            output.close()
            
            output = mzReport.writer(outputfile, columns = report_columns,
                                     sheet_name = 'Data')
            for psm in psms:
                output.write(psm)
            output.close()
        
        else:
            #report_columns.insert(0, 'File')
            if not outputfile:
                raise IOError("Combined report file name must be specified.")
            
            if (outputfile.lower().endswith('xls') or
                outputfile.lower().endswith('xlsx') or
                outputfile.lower().endswith('mzd')):            
                for m_id, datafilepath, header, _ in reports:
                    datafilename = os.path.basename(datafilepath)
                    output = mzReport.writer(outputfile, columns = header[0],
                                    sheet_name = '%s Mascot Header' % datafilename[:17])
                                    # sheet_name has to be <= 31 characters long.
                    for line in header[1:]:
                        output.write(line)
                    output.close()
            else:
                extension = outputfile.split('.')[-1]
                print(("Omitting header tables due to %s format." % extension))
            
            output = mzReport.writer(outputfile, columns = ['File'] + report_columns,
                                     sheet_name = 'Data')
            for _, datafilepath, _, psms in reports:
                datafilename = os.path.basename(datafilepath)
                for psm in psms:
                    psm['File'] = datafilename
                    output.write(psm)
            
            output.close()


        return outputfile










class mascot(object):
    def __init__(self, server=None, version=2.1, verbose=False):
        '''Initialize a Mascot search session. This starts up Curl, but does not log in
        to the server.
        '''
        self.server = str(server) if server else server
        self.version = float(version)
        self.verbose = verbose

        self.logged_in = False
        
        self.loginToken = ''

    def login(self, login, password):
        """
        Logs in to mascot. Returns "success" or "error" based on login success
        """

        login_url = self.server.strip(r'\/') + r'//cgi/login.pl'

        login_form_seq = {'username':login,
                          'password':password,
                          'submit':'Login',
                          'display':'logout_prompt',
                          'savecookie':'1',
                          'action':'login',
                          'userid': '',
                          'onerrdisplay':'login_prompt'}
                                   
        #req = urllib.request.Request(login_url, 
                              #data = urlencode(login_form_seq))
        #f = urllib.request.urlopen(req)
        f = requests.post(login_url, data = login_form_seq)
        try:
            self.loginToken = f.cookies # f.info()['Set-Cookie']
            # I suppose in this case precise disposition of the output is less important.
        except KeyError:
            raise RuntimeError('Login was not successful; check your username and password.')
        
        
        err = None
        login_re = re.compile(r'<B><FONT COLOR=#FF0000>Error:(.*)</FONT></B>')
        for line in f.text.splitlines():
            m = login_re.search(line)
            if m:
                err = (m.group(1)[:-7] if m.group(1).endswith(' [ec=1]') else m.group(1))
        f.close()
        

        if err:
            logger_message(40, "Error logging in to server: %s" % err)
            self.logged_in = False
            return 'error'
        else:
            self.logged_in = True
            return 'success'

    def get_date(self, mascot_id):
        """
        Returns search date string corresponding to mascot ID
        
        BE CAREFUL- since the Mascot search log takes a few minutes to
        update, this is likely to fail when called against newly searched
        data.
        """

        if not self.logged_in:
            raise IOError("Not logged in to Mascot server.  (Call .login() first!)")

        mascot_id = re.sub(r'^0+', '', str(mascot_id), count=1)

        date_url = self.server + r'/x-cgi/ms-review.exe?'
        
        #req = urllib.request.Request(date_url + urlencode([('CalledFromForm', '1'),
                                                    #('f0', mascot_id)]))
        #req.add_header("Cookie", self.loginToken)
        #f = urllib.request.urlopen(req)
        ##response.write(f.read())
        res = requests.post(date_url, 
                            data = {'CalledFromForm':'1', 'f0':mascot_id},
                            cookies = self.loginToken)

        date_re = re.compile(r'.+/data/(\d+)/F%s.+' % mascot_id.zfill(6))

        matches = []
        #for line in response.getvalue().splitlines():
        for line in res.text.splitlines():
            m = date_re.match(line)
            if m:
                matches.append(m)
                
                
        if matches:
            return matches[-1].group(1)
        else:
            raise RuntimeError("Search not found in date lookup.")
        
    
    def download_dat(self, chosen_folder, mascot_id, date=None):
        mascot_id = str(mascot_id).zfill(6)
        dat_file = 'F%s.dat' % mascot_id

        dat_path = os.path.abspath(os.path.join(chosen_folder, dat_file))
        
        dat_url = self.server + r'/cgi/export_dat_2.pl?'
        #dat_url += ';'.join(['='.join([k, v]) for k, v in [('file', r'../data/%s/F%s.dat' % (date, mascot_id)),
                                                           #('export_format', 'MascotDAT'),
                                                           #('generate_file', '1'),
                                                           #('do_export', '1')]])
                                                           
        res = requests.post(dat_url,
                            data = {'do_export': '1',
                                    'export_format': 'MascotDAT',
                                    'file': r'../data/%s/F%s.dat' % (date, mascot_id),
                                    'generate_file': '1',
                                    'generate_mime': '1'},
                            cookies = self.loginToken)
        print(len(res.text))
        datOut = open(dat_path, 'w')
        datOut.write(res.text)
        datOut.close()   
        
        return dat_path 

    def download_csv(self, mascot_id, **etcetc):
        return self.download_file(mascot_id, filetype = 'CSV', **etcetc)



    def download_file(self, mascot_id, filetype = 'CSV',
                     max_hits=999999, ion_cutoff=20, bold_red=True,
                     unassigned_queries=False, show_query_data=True,
                     show_same_set=False, show_sub_set=False, date=None,
                     ion_list=False, quant=False,
                     save_file=None,
                     buffer_results= False,
                     **kwargs):
        """Downloads mascot csv file for a given id and parameters.
        Returns csv as string.

        """

        extension = filetypes[filetype]

        mascot_id = str(mascot_id).zfill(6)

        if not date:
            date = self.get_date(mascot_id)            
        assert date, "Could not retrieve date of search ID %s" % mascot_id
            
        download_url = self.server.strip(r'/') + r'/cgi/export_dat_2.pl'

        download_form_seq = [('file', r'../data/%s/F%s.dat' % (date, mascot_id)),
                             ('export_format', filetype),
                             ('REPORT', str(max_hits)),
                             ('_ignoreionsscorebelow', str(ion_cutoff)),
                             ('_requireboldred', int(bold_red)),
                             ('show_unassigned', int(unassigned_queries)),
                             ('show_same_sets', int(show_same_set)),
                             ('_showsubsets', int(show_sub_set)),
                             ('_sigthreshold', '0.5'),
                             ('pep_calc_mr', '1'), ('pep_delta', '1'),
                             ('pep_end', '1'), ('pep_exp_mr', '1'),
                             ('pep_exp_mz', '1'), ('pep_exp_z', '1'),
                             ('pep_expect', '1'), ('pep_isbold', '1'),
                             ('pep_miss', '1'), ('pep_query', '1'),
                             ('pep_rank', '1'), ('pep_score', '1'),
                             ('pep_seq', '1'), ('pep_start', '1'),
                             ('pep_var_mod', '1'), ('peptide_master', '1'),
                             ('prot_acc', '1'), ('prot_desc', '1'),
                             ('prot_hit_num', '1'), ('prot_mass', '1'),
                             ('prot_matches', '1'), ('prot_score', '1'),
                             ('protein_master', '1'), ('search_master', '1'),
                             ('show_format', '1'), ('show_header', '1'),
                             ('show_mods', '1'), ('show_params', '1'),
                             
                             ('show_pep_dupes', '1'),
                             
                             ('generate_file', 1),
                             ('sessionid', self.loginToken['MASCOT_SESSION']),
                             ('do_export', 1),
                             ('SEARCH', 'MIS'),
                             
                             ('_server_mudpit_switch', '99999999'),
                             ('pep_scan_title', int(show_query_data)),
                             ('query_master', int(ion_list)),
                             ('pep_quant', int(quant)),
                             ('percolate', '0')]
  


        
        #assert False
        # For some reason, export_dat_2.pl doesn't play properly with the
        # Python Requests package, so a hack using urllib has to be put
        # in.
        
        try:
            # Python 3
            from urllib.request import Request, urlopen
            from urllib.parse import urlencode
            from http.client import parse_headers
            encoded_form = bytes(urlencode(download_form_seq), 'utf-8')
            request_call = {'url':download_url,
                            'data':encoded_form,
                            'method':'POST'}
        except ImportError: 
            # Python 2.
            from urllib2 import Request, urlopen
            from urllib import urlencode
            parse_headers = lambda x: None
            encoded_form = urlencode(download_form_seq)
            request_call = {'url':download_url,
                            'data':encoded_form}            
            
        for i in range(10):
            #res = requests.post(download_url, params = dict(download_form_seq), cookies = self.loginToken, allow_redirects = True, stream = True, verify = False)
            
            #text = res.text
            req = Request(**request_call)
            f = urlopen(req, timeout = 999999)   
            headers = parse_headers(f)
            chunks = []
            for chunk in f:
                chunks.append(chunk.decode('utf-8'))
            text = ''.join(chunks)
            
            if ('doctype html public' not in text.lower() and '<html>' not in text and 'javascript' not in text):
                break
            
        assert ('DOCTYPE HTML PUBLIC' not in text and '<html>' not in text)
        assert text
        #assert False
        
        if buffer_results:
            result_buffer = io.StringIO()
            result_buffer.write(text)
            result_buffer.seek(0)
            return result_buffer
        else:
            if not save_file:
                save_file = os.path.join(os.getcwd(), "F%s.%s" % (mascot_id, extension))
            with open(save_file, 'w') as output:
                output.write(text)
            return save_file
    
    
    
    def download_mzid(self, mascot_id,
                     max_hits=1000, ion_cutoff=20, bold_red=True,
                     unassigned_queries=False, show_query_data=True,
                     show_same_set=False, show_sub_set=False, date=None,
                     ion_list=False, quant=False,
                     save_file=None, **kwargs):
        """Downloads mascot csv file for a given id and parameters.
        Returns csv as string.

        """

        if not self.logged_in:
            raise IOError("Not logged into to Mascot server.")

        mascot_id = str(mascot_id).zfill(6)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            raise IOError("Couldn't locate file on Mascot server.")

        if self.version == 2.1:
            download_url = self.server + r'/cgi/export_dat.pl?'
        elif self.version >= 2.2: # not sure about 2.3 here
            download_url = self.server + r'/cgi/export_dat_2.pl?'

        download_form_seq = [('file', r'../data/%s/F%s.dat' % (date, mascot_id)),
                             ('export_format', 'mzIdentML'),
                             ('REPORT', str(max_hits)),
                             ('_ignoreionsscorebelow', str(ion_cutoff)),
                             ('_requireboldred', int(bold_red)),
                             ('show_unassigned', int(unassigned_queries)),
                             ('show_same_sets', int(show_same_set)),
                             ('_showsubsets', int(show_sub_set)),
                             ('_sigthreshold', '0.05'), ('do_export', '1'),
                             ('pep_calc_mr', '1'), ('pep_delta', '1'),
                             ('pep_end', '1'), ('pep_exp_mr', '1'),
                             ('pep_exp_mz', '1'), ('pep_exp_z', '1'),
                             ('pep_expect', '1'), ('pep_isbold', '1'),
                             ('pep_miss', '1'), ('pep_query', '1'),
                             ('pep_rank', '1'), ('pep_score', '1'),
                             ('pep_seq', '1'), ('pep_start', '1'),
                             ('pep_var_mod', '1'), ('peptide_master', '1'),
                             ('prot_acc', '1'), ('prot_desc', '1'),
                             ('prot_hit_num', '1'), ('prot_mass', '1'),
                             ('prot_matches', '1'), ('prot_score', '1'),
                             ('protein_master', '1'), ('search_master', '1'),
                             ('show_format', '1'), ('show_header', '1'),
                             ('show_mods', '1'), ('show_params', '1'),
                             ('generate_file', '1')]

        if self.version == 2.1:
            download_form_seq.extend([('_mudpit', '99999999'),
                                      ('show_queries', int(show_query_data))])
        elif self.version >= 2.2: # not sure about 2.3 here...should be the same?
            download_form_seq.extend([('_server_mudpit_switch', '99999999'),
                                      ('pep_scan_title', int(show_query_data)),
                                      ('query_master', int(ion_list)),
                                      ('pep_quant', int(quant))])

        res = requests.post(download_url, 
                            data = dict(download_form_seq),
                            cookies = self.loginToken)

        if not save_file:
            save_file = os.path.join(os.getcwd(), "F%s.mzid" % mascot_id)

        with open(save_file, "w") as f:
            f.write(res.text)
        res.close()

        return save_file    

    def clean_csv(self, input_file, export_file=None, ion_list=False, mascot_var_mods=True):
        """
        Converts a mascot csv file into a tabular csv file
        """

        if not export_file:
            export_file_split = os.path.split(input_file)
            export_file = os.path.join(export_file_split[0], export_file_split[1][:-4] + "_clean.csv")

        non_tabular = ['prot_acc', 'prot_desc', 'prot_score', 'prot_mass',
                       'prot_matches', 'prot_cover', 'prot_len', 'prot_pi',
                       'prot_tax_str', 'prot_tax_id']

        if self.version == 2.1:
            hit_columns = [('prot_hit_num', int), ('prot_acc', str),
                           ('prot_desc', str), ('prot_score', int),
                           ('prot_mass', int), ('prot_matches', int),
                           ('prot_cover', float), ('prot_len', int),
                           ('prot_pi', float), ('prot_tax_str', str),
                           ('prot_tax_id', str), ('pep_query', int),
                           ('pep_exp_mz', float), ('pep_exp_mr', float),
                           ('pep_exp_z', int), ('pep_calc_mr', float),
                           ('pep_delta', float), ('pep_start', int),
                           ('pep_end', int), ('pep_miss', int),
                           ('pep_score', float), ('pep_homol', float),
                           ('pep_ident', float), ('pep_expect', float),
                           ('pep_rank', int), ('pep_res_before', str),
                           ('pep_seq', str), ('pep_res_after', str),
                           ('pep_frame', int), ('pep_var_mod', str)]

            scan_columns = [('query_number', int), ('StringTitle', str),
                            ('qualifiers', str), ('TotalIonsIntensity', float),
                            ('NumVals', int), ('StringIons1', str),
                            ('StringIons2', str), ('StringIons3', str)]
        elif self.version >= 2.2: # 2.3 should be the same here, I think
            non_tabular.append('prot_seq')
            hit_columns = [('prot_hit_num', int), ('prot_acc', str),
                           ('prot_desc', str), ('prot_score', int),
                           ('prot_mass', int), ('prot_matches', int),
                           ('prot_cover', float), ('prot_len', int),
                           ('prot_pi', float), ('prot_tax_str', str),
                           ('prot_tax_id', str), ('prot_seq', str),
                           ('pep_query', int), ('pep_rank', int),
                           ('pep_isbold', int), ('pep_exp_mz', float),
                           ('pep_exp_mr', float), ('pep_exp_z', int),
                           ('pep_calc_mr', float), ('pep_delta', float),
                           ('pep_start', int), ('pep_end', int),
                           ('pep_miss', int), ('pep_score', float),
                           ('pep_homol', float), ('pep_ident', float),
                           ('pep_expect', float), ('pep_res_before', str),
                           ('pep_seq', str), ('pep_res_after', str),
                           ('pep_frame', int), ('pep_var_mod', str),
                           ('pep_var_mod_pos', str), ('pep_num_match', int),
                           ('pep_scan_title', str)]

            scan_columns = [('query_number', int), ('moverz', float),
                            ('charge', str), ('intensity', float),
                            ('StringTitle', str), ('Scan number range', str),
                            ('Retention time range', str), ('qualifiers', str),
                            ('Peptide Mass Tolerance',str), ('Peptide Mass Tolerance Units',str),
                            ('Variable modifications', str), ('Instrument type', str),
                            ('TotalIonsIntensity', float), ('NumVals', int),
                            ('StringIons1', str), ('StringIons2', str),
                            ('StringIons3', str), ('pep_rank', int),
                            ('pep_isbold', int), ('pep_exp_mz', float),
                            ('pep_exp_mr', float), ('pep_exp_z', int),
                            ('pep_calc_mr', float), ('pep_delta', float),
                            ('pep_start', int), ('pep_end', int),
                            ('pep_miss', int), ('pep_score', float),
                            ('pep_homol', float), ('pep_ident', float),
                            ('pep_expect', float), ('pep_res_before', str),
                            ('pep_seq', str), ('pep_res_after', str),
                            ('pep_frame', int), ('pep_var_mod', str),
                            ('pep_var_mod_pos', str), ('pep_num_match', int),
                            ('pep_scan_title', str)]
        else:
            return

        #only for 2.2 (+?)
        itraq_columns = []
        itraq_row1 = 0
        itraq_ratios_per_prot_query = {}

        #Start parsing
        reader = csv.reader(open(input_file, "r"))
        headers = next(reader)

        if mascot_var_mods and self.version >= "2.2":
            while (len(headers) == 0 or (headers[0] != "Identifier" and headers[0] != 'prot_hit_num')):
                headers = next(reader)
            if headers[0] != 'prot_hit_num':
                var_mod = next(reader)
                var_mod_dict = {}
                while len(var_mod) > 0:
                    var_mod_dict[var_mod[0]] = var_mod[1].split('(')[0][:-1]
                    var_mod = next(reader)
                headers = var_mod
            else:
                var_mod_dict = None

        while len(headers) == 0 or headers[0] != "prot_hit_num":
            headers = next(reader)

        #potentially missing protein related entries which may need to be duplicated across rows...
        missing = 10
        if self.version >= 2.2: # presumably this should work for 2.3
            missing = 11
        counter = 0
        #index to scan info, if not present then no point in trying to join with query data...
        query_header = 0
        while counter < len(hit_columns):
            if counter >= len(headers):
                break
            if headers[counter] == 'pep_query':
                query_header = counter
            if headers[counter] == hit_columns[counter][0]:
                counter += 1
            else:
                if hit_columns[counter][0] in non_tabular:
                    missing -= 1
                del(hit_columns[counter])

        #Delete remaining columns from hit_columns:
        remColNum = len(hit_columns) - len(headers)
        for i in range(remColNum):
            del(hit_columns[-1])

        #No unknown columns
        assert len(headers) == len(hit_columns)

        #Because the layout is non-tabular... (human-oriented).
        missingVals = [None] * missing
        protHit = 0
        entries = []
        for row in reader:
            if len(row) == 0 or not row[0]:
                break
            if row[0] != protHit or row[1] != '':
                protHit = row[0]
                missingVals[0:missing] = row[1:(missing+1)]
            else:
                row[1:(1+missing)] = missingVals[:]

            pep_query = row[query_header]
            entry = [None] * len(hit_columns)
            # Error in Mascot: missing the last comma...
            if len(row) < len(hit_columns):
                for addCom in range(len(hit_columns)-len(row)):
                    row.append("")
            elif len(row) > len(hit_columns):
                # maybe it has some itraq info
                # Big Assumption that after all the hit columns in the row,
                # we have quant info in the format Header1,Val1,Header2,Val2
                # (ex: "115/114",11.099,"116/114",5.243,"117/114",13.046)
                itraq_info = row[len(hit_columns):]
                itraq_data = []

                for itraq_index in range(0, len(itraq_info), 2):
                    if not itraq_row1:
                        itraq_columns.append(("Quant " + itraq_info[itraq_index], float))
                    if itraq_info[itraq_index+1][-1] == "-":
                        row[len(hit_columns) + len(itraq_data)] = "0"
                        itraq_ratios_per_prot_query[(protHit, pep_query, len(itraq_data))] = 0
                    else:
                        row[len(hit_columns) + len(itraq_data)] = itraq_info[itraq_index+1]
                        itraq_ratios_per_prot_query[(protHit, pep_query, len(itraq_data))] = itraq_info[itraq_index+1]
                    itraq_data.append(itraq_info[itraq_index+1])

                entry+= [None] * len(itraq_columns)
                itraq_row1+=1

            # converting to correct format
            for counter in range(len(headers)):
                if row[counter]:
                    entry[counter] = hit_columns[counter][1](row[counter])
                else:
                    entry[counter] = ''

            # converting to correct format - itraq
            for counter in range(len(itraq_columns)):
                if len(headers)+counter >= len(entry):
                    entry.append(itraq_columns[counter][1](itraq_ratios_per_prot_query[(protHit, pep_query, counter)]))
                else:
                    entry[len(headers)+counter] = itraq_columns[counter][1](row[len(headers)+counter])

            entries.append(entry)

        headers = []
        scans = []
        try:
            while len(headers) == 0 or headers[0] != "query_number":
                headers = next(reader)
        except StopIteration:
            pass

        if headers:
            counter = 0
            while counter < len(scan_columns):
                if counter >= len(headers):
                    break

                if headers[counter] == scan_columns[counter][0]:
                    counter += 1
                else:
                    del(scan_columns[counter])

            #Delete remaining columns from hit_columns:
            remColNum = len(scan_columns) - len(headers)
            for i in range(remColNum):
                del(scan_columns[-1])

            #No unknown columns
            assert len(headers) == len(scan_columns)

            for row in reader:
                if len(row) == 0 or not row[0]:
                    break
                entry = [None] * len(scan_columns)
                # Error in Mascot: missing the last comma(s)...
                if len(row) < len(scan_columns):
                    for i in range(len(scan_columns) - len(row)):
                        row.append("")
                elif len(row) > len(scan_columns):
                    row = row[:len(scan_columns)]

                for counter in range(len(scan_columns)):
                    if row[counter] == '':
                        entry[counter] = row[counter]
                    else:
                        entry[counter] = scan_columns[counter][1](row[counter])


                scans.append(entry)
          #Apparently there is no way to close the reader...

        if len(scans) == 0 or query_header == 0:
            data = entries
            columns = hit_columns
            columns.extend(itraq_columns)  #itraq
        else:
            data = []
            columns = hit_columns
            columns.extend(itraq_columns)  #itraq
            columns.extend(scan_columns[1:])
            for entry in entries:
                data.append( entry + scans[ entry[query_header] -1][1:] )

        #write to export_file
        keep_columns = mzReport.default_columns[:]

        if ion_list:
            keep_columns.append(mzReport.multiplierzHeaders['ion_str'])

        #for iTRAQ
        keep_columns.extend(name for name,t in itraq_columns)

        #Rename columns
        temp_columns = set()
        counter = -1
        headerPos = {}
        for (name,t) in columns:
            counter+=1
            newName = name
            if name in multiplierzHeaderConversions:
                newName = multiplierzHeaderConversions[name]
            headerPos[newName] = counter
            if newName not in keep_columns:
                continue
            temp_columns.add(newName)

        #ReOrder Column Headers
        newColumns = [name for name in keep_columns if name in temp_columns]

        # process var_mods if we have position info
        if mascot_var_mods and self.version >= 2.2 and var_mod_dict:
            for entry in data:
                seq = entry[headerPos['Peptide Sequence']]
                var_mod_pos = entry[headerPos['pep_var_mod_pos']]
                if var_mod_pos:
                    entry[headerPos['Variable Modifications']] = self.parse_mod_pos(seq, var_mod_pos, var_mod_dict)

        #ReOrder Data
        newData = [[entry[headerPos[col]] for col in newColumns] for entry in data]

        with open(export_file, 'wb') as out_file:
            writer = csv.writer(out_file, dialect='excel')
            writer.writerow(newColumns)
            for entry in newData:
                writer.writerow(entry)



    def parse_mod_pos(self, sequence, var_mod_pos, var_mod_dict):
        if var_mod_pos[0] != '0':
            var_mods = ['N-term: %s' % var_mod_dict[var_mod_pos[0]]]
        else:
            var_mods = []

        seq_mods = var_mod_pos[2:-2]
        var_mods.extend('%s%d: %s' % (aa,i+1,var_mod_dict[seq_mods[i]]) for i,aa in enumerate(sequence)
                        if seq_mods[i] != '0')

        if var_mod_pos[-1] != '0':
            var_mods.append('C-term: %s' % var_mod_dict[var_mod_pos[-1]])

        return '; '.join(var_mods)

    def close(self, ):
        '''Close this Mascot search session.
        '''
        pass # There could be something to do with closing the session here.


def parse_fields(text):
    from bs4 import BeautifulSoup as BS      # Ahem.
    
    parse = BS(text, 'html.parser')
    options = parse.find_all("option")
    
    fields = {}
    for opt_type, options in collectByCriterion(options,
                                                lambda x: x.parent.attrs['name']).items():
        fields[opt_type] = [x.text for x in options]
    
    # Mod lists are stored in Javascript junk, so BS doesn't help us there,
    # but it's pretty simple.
    shortmods = []
    longmods = []
    get_mod = re.compile("'(.*)'")
    for line in text.split('\n'):
        if line.startswith('mods['):
            try:
                modstr = get_mod.search(line).group(1)
                shortmods.append(modstr)
            except AttributeError:
                print("Parse error: %s" % line)
        elif line.startswith('allMods['):
            try:
                modstr = get_mod.search(line).group(1)
                longmods.append(modstr)    
            except AttributeError:
                print("Parse error: %s" % line)
    
    fields['MODS'] = shortmods
    fields['IT_MODS'] = shortmods
    fields['ALLMODS'] = longmods
    
    return fields
    
        
        
class MascotSearcher(object):
    def __init__(self, server = None, version = None, verbose = True,
                 open_tabs = False):
        
        if server:
            self.server = str(server)
            self.version = float(version)
        else:
            self.server = settings.mascot_server
            self.version = float(settings.mascot_version)
            
        self.verbose = verbose
        self.open_tabs = open_tabs
        
        self.cookies = {}
    
    def login(self, login, password):
        login_url = self.server + r'/cgi/login.pl'            
        
        login_form_seq = [('username', login),
                          ('password', password),
                          ('submit', 'Login'),
                          ('display', 'logout_prompt'),
                          ('savecookie', '1'),
                          ('action', 'login'),
                          ('userid', ''),
                          ('onerrdisplay', 'login_prompt')]
        login_form = dict(login_form_seq)
        
        r = requests.post(login_url, data = login_form)
        self.cookies = r.cookies
        
        
        self.login_data = r.text
        
        err = None
        login_re = re.compile(r'<B><FONT COLOR=#FF0000>Error:(.*)</FONT></B>')
        for line in r.text.splitlines():
            m = login_re.search(line)
            if m:
                err = (m.group(1)[:-7] if m.group(1).endswith(' [ec=1]') else m.group(1))
    
        if err:
            logger_message(40, "Error logging in to server: %s" % err)
            self.logged_in = False
            return 1
        else:
            self.logged_in = True
            return 0 # Old-style return codes to deal with RC_MascotSearch.    
        
        
    def get_fields(self, version):
        search_form_url = self.server + r'/cgi/search_form.pl'        
        r = requests.get(search_form_url + "?" + 'FORMVER=1.01&SEARCH=MIS',
                         cookies = self.cookies)
        
        return parse_fields(r.text)
          
    
    def get_fields_better(self):
        # Doesn't seem to return a bunch of fields that RC_MascotSearch wants.
        
        class fieldParser(object):
            def __init__(self):
                self.content = io.StringIO()
            def __call__(self, content):
                self.content.write(content)
            def parse_fields(self):
                fields = {}
                curKey = None
                curValues = []
                for line in self.content.getvalue().split('\n'):
                    line = line.strip()
                    if line and line[0] == "[" and line[-1] == "]":
                        if curKey:
                            fields[curKey] = curValues
                    
                        curKey = line[1:-1]
                        curValues = []
                    elif line:
                        curValues.append(line)
                fields[curKey] = curValues
                return fields
                
        fieldFormURL = self.server + r"/cgi/get_params.pl"

        r = requests.get(fieldFormURL, cookies = self.cookies)
        fieldParserInstance = fieldParser()
        fieldParserInstance(r.text)
        
        fieldDict = fieldParserInstance.parse_fields()
        
        if 'PFA' not in fieldDict:
            fieldDict['PFA'] = list(range(0, 9))
        if 'IT_MODS' not in fieldDict:
            fieldDict['IT_MODS'] = fieldDict['MODS']
        
        return fieldDict
        
    
    def search(self, search_form, target_file = None):
        search_url = self.server + r'/cgi/nph-mascot.exe?1'
        
        search_form = dict(search_form)
        search_form.update([('FORMVER','1.01'), ('SEARCH','MIS'), ('MULTI_SITE_MODS', '1')])
        #search_form.update([('LIBRARY_SEARCH', '1')])
        if target_file:
            filename = target_file
        else:
            filename = search_form['FILE']
        assert os.path.exists(filename)
        
        searchFile = {'FILE' : open(filename, 'r')}
        
        r = requests.post(search_url,
                          data = search_form,
                          files = searchFile,
                          cookies = dict(self.cookies),
                          stream = True)
        #print r.json()['url']
        
        dat_file_id, err = (None, None)
        rec = []
        for line in r.text.splitlines():
            rec.append(line)
            res = self.parse_response(line)
            if isinstance(res,tuple):
                err = res
            elif res:
                dat_file_id = res
                
        if not dat_file_id:
            print('\n'.join(rec))
                
        return (dat_file_id,err)     
    
    def parse_response(self, s):
        m = re.search(r'A HREF="\.\.(/cgi/master_results.+/data/(\d+)/F(.+)\.dat)">Click here to see Search Report',s)
        m2 = re.search(r'.*\[M\d+\]',s)
        if m:
            if self.open_tabs:
                print("Opening web browser to result page.")
                webbrowser.open_new_tab(self.server + m.group(1))
            return ('%s:%s' % (m.group(3), m.group(2)))
        elif m2:
            return ('Search could not be performed', m2.group(0))
        return None        
    
    def close(self):
        pass
    
    
        
        
        


class MascotDatFile(object):
    def __init__(self, dat_file_path, decoyMode = False, **kwargs):
        #max_hits, ion_cutoff,
        #bold_red, unassigned_queries, show_query_data,
        #show_same_set, show_sub_set, quant

        self.dat_file_path = str(dat_file_path)
        self.res_file = ms.ms_mascotresfile(str(dat_file_path))

        if not self.res_file.isValid():
            print(dat_file_path)
            errors = []
            for i in range(1, self.res_file.getNumberOfErrors() + 1):
                errors.append(self.res_file.getErrorString(i))
                
            if any(['missing or corrupt headers' in x for x in errors]):
                errors.append("\n(Typically this occurs because Mascot couldn't find \n"
                              "the requested search and returned an invalid .DAT file.)")
            
            raise IOError("MSParser errors: \n %s" % '\n'.join(errors))
        
        self.params = self.res_file.params()
        self.args = kwargs
        if 'mascot_var_mods' not in self.args:
            self.args['mascot_var_mods'] = True

        #if self.res_file.isErrorTolerant():
            #raise ValueError("Error Tolerant searches are not supported yet")
        elif 'quant' in self.args and self.args['quant']:
            raise ValueError("Peptide Quantitation isn't supported yet")

        msres_flags = (ms.ms_mascotresults.MSRES_MAXHITS_OVERRIDES_MINPROB
                       | ms.ms_mascotresults.MSRES_GROUP_PROTEINS)

        if self.args.get('show_sub_set', default_show_sub_set):
            msres_flags = msres_flags | ms.ms_mascotresults.MSRES_SHOW_SUBSETS

        if self.args.get('bold_red', default_bold_red):
            msres_flags = msres_flags | ms.ms_mascotresults.MSRES_REQUIRE_BOLD_RED
            
        if decoyMode:
            msres_flags = msres_flags | ms.ms_mascotresults.MSRES_DECOY
            
        if self.res_file.isErrorTolerant():
            print("Experimental error-tolerant search results feature!")
            msres_flags = msres_flags | ms.ms_mascotresults.MSRES_INTEGRATED_ERR_TOL

        #mspepsumFlags = ms.ms_peptidesummary.MSPEPSUM_PERCOLATOR
        mspepsumFlags = ms.ms_peptidesummary.MSPEPSUM_NONE

        self.pep_summary = ms.ms_peptidesummary(self.res_file,
                                                msres_flags,
                                                0.05,
                                                self.args.get('max_hits', default_max_hits),
                                                '',
                                                self.args.get('ion_cutoff', default_ion_cutoff),
                                                0,
                                                '',
                                                mspepsumFlags)

        self.var_mod_dict = {}
        for i in range(1, len(self.params.getIT_MODS().split(',')) + 1):
            if i > 15:
                raise ValueError("Too many modifications!")
            self.var_mod_dict[str(hex(i))[2:].upper()] = self.params.getVarModsName(i).rsplit('(', 1)[0][:-1]

        if self.params.getQUANTITATION():
            i += 1
            while self.params.getVarModsName(i):
                if i > 15:
                    raise ValueError("Too many modifications!")                
                self.var_mod_dict[str(hex(i))[2:].upper()] = self.params.getVarModsName(i).rsplit('(', 1)[0][:-1]
                i += 1

    def close(self):
        self.res_file = None

    def getFILENAME(self):
        filenamestr = self.params.getFILENAME()
        if 'File Name' in filenamestr and 'File Path' in filenamestr:
            # Proteome Discoverer file list.
            filename = filenamestr.split(';')[0].split()[-1]
            assert filename.lower().endswith('.raw'), filename
            filename = os.path.basename(filename)
            datnum = os.path.basename(self.dat_file_path)[1:-4]
            filename = filename[:-3] + datnum + '.raw'
            
            return filename
        else:
            return filenamestr
            
    def mascot_header(self):
        mascot_header = [['Header', '-' * 50]]
        mascot_header.append([' ', ' '])

        mascot_header.append(['Search title', self.params.getCOM()])
        mascot_header.append(['Timestamp',
                              time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(self.res_file.getDate()))])
        mascot_header.append(['User', self.params.getUSERNAME()])
        mascot_header.append(['Email', self.params.getUSEREMAIL()])
        mascot_header.append(['Report Data', os.path.basename(self.dat_file_path)])

        mascot_header.append(['Peak list data path', self.getFILENAME()])
        mascot_header.append(['Peak list format', self.params.getFORMAT()])
        mascot_header.append(['Search type', self.params.getSEARCH()])
        mascot_header.append(['Mascot version', self.res_file.getMascotVer()])
        if self.res_file.getMascotVer() >= 2.3:
            mascot_header.append(['Database', ' '.join('%d::%s' % (i+1, self.params.getDB(i+1))
                                                       for i in range(self.params.getNumberOfDatabases()))])
            mascot_header.append(['Fasta file', ' '.join('%d::%s' % (i+1, self.res_file.getFastaVer(i+1))
                                                         for i in range(self.params.getNumberOfDatabases()))])
        else:
            mascot_header.append(['Database', self.params.getDB()])
            mascot_header.append(['Fasta file', self.res_file.getFastaVer()])

        mascot_header.append(['Total sequences', self.res_file.getNumSeqs()])
        mascot_header.append(['Total residues', self.res_file.getNumResidues()])
        mascot_header.append(['Sequences after taxonomy filter', self.res_file.getNumSeqsAfterTax()])
        mascot_header.append(['Number of queries', self.res_file.getNumQueries()])

        mascot_header.append([' ', ' '])
        mascot_header.append(['Variable modifications', '-' * 50])
        mascot_header.append([' ', ' '])

        mascot_header.append(['Identifier', 'Name,Delta,Neutral loss(es)'])

        #self.var_mod_dict = {}
        for i in range(1, len(self.params.getIT_MODS().split(',')) + 1):
            mascot_header.append([i, ','.join([self.params.getVarModsName(i), '%s' % self.params.getVarModsDelta(i)]
                                              + ['%s' % f for f in self.params.getVarModsNeutralLosses(i)])])
            #self.var_mod_dict[str(i)] = self.params.getVarModsName(i).split('(')[0][:-1]

        mascot_header.append([' ', ' '])
        mascot_header.append(['Search Parameters', '-' * 50])
        mascot_header.append([' ', ' '])

        mascot_header.append(['Taxonomy filter', self.params.getTAXONOMY()])
        mascot_header.append(['Enzyme', self.params.getCLE()])
        mascot_header.append(['Maximum Missed Cleavages', self.params.getPFA()])
        mascot_header.append(['Fixed modifications', self.params.getMODS()])
        if self.params.getQUANTITATION():
            mascot_header.append(['Quantitation method', self.params.getQUANTITATION()])
        else:
            mascot_header.append(['ICAT experiment', int(self.params.getICAT())])
        mascot_header.append(['Variable modifications', self.params.getIT_MODS()])
        mascot_header.append(['Peptide Mass Tolerance', self.params.getTOL()])
        mascot_header.append(['Peptide Mass Tolerance Units', self.params.getTOLU()])
        mascot_header.append(['Fragment Mass Tolerance', self.params.getITOL()])
        mascot_header.append(['Fragment Mass Tolerance Units', self.params.getITOLU()])
        mascot_header.append(['Mass values', self.params.getMASS()])
        if self.params.getSEG() > 0:
            mascot_header.append(['Protein Mass', self.params.getSEG()])
        if self.params.getINSTRUMENT():
            mascot_header.append(['Instrument type', self.params.getINSTRUMENT()])
        if self.params.getPEP_ISOTOPE_ERROR() > -1:
            mascot_header.append(['Isotope error mode', self.params.getPEP_ISOTOPE_ERROR()])

        mascot_header.append([' ', ' '])
        mascot_header.append(['Format parameters', '-' * 50])
        mascot_header.append([' ', ' '])

        #mascot_header.append(['Significance threshold', 0.05])
        mascot_header.append(['Max. number of hits', self.args.get('max_hits', default_max_hits)])
        mascot_header.append(['Ions score cut-off', self.args.get('ion_cutoff', default_ion_cutoff)])
        mascot_header.append(['Include same-set proteins', int(self.args.get('show_same_set', default_show_same_set))])
        mascot_header.append(['Include sub-set proteins', int(self.args.get('show_sub_set', default_show_sub_set))])
        mascot_header.append(['Include unassigned', int(self.args.get('unassigned_queries', False))])
        mascot_header.append(['Require bold red', int(self.args.get('bold_red', default_bold_red))])

        return mascot_header

    def parse_mod_pos(self, sequence, var_mod_pos, queryPeptide):
        if var_mod_pos[0] != '0':
            try:
                var_mods = ['N-term: %s' % self.var_mod_dict[var_mod_pos[0]]]
            except KeyError:
                assert var_mod_pos[0] == 'X'
                var_mods = ['N-term: %s' % str(self.pep_summary.getErrTolModName(*queryPeptide))]
        else:
            var_mods = []

        seq_mods = var_mod_pos[1:-1]
        #var_mods.extend('%s%d: %s' % (aa,i+1, self.var_mod_dict[seq_mods[i]])
                        #for i,aa in enumerate(sequence)
                        #if seq_mods[i] != '0')
        for i, aa in enumerate(sequence):
            try:
                if seq_mods[i] != '0':
                    var_mods.append('%s%d: %s' % (aa,i+1, self.var_mod_dict[seq_mods[i]]))
            except KeyError:
                assert seq_mods[i] == 'X' # Error-tolerant imputed modification.
                #assert seq_mods.count('X') == 1 # Pretty sure there's only one error allowed per peptide.
                var_mods.append('%s%d: %s' % (aa,i+1,
                                              str(self.pep_summary.getErrTolModName(*queryPeptide))))
                
        if var_mod_pos[-1] != '0':
            try:
                var_mods.append('C-term: %s' % self.var_mod_dict[var_mod_pos[-1]])
            except KeyError:
                assert var_mod_pos[-1] == 'X'
                var_mods.append('C-term: %s' % str(self.pep_summary.getErrTolModName(*queryPeptide)))

        return '; '.join(var_mods)

    def output_protein(self, prot):
        prot_hit_num = prot.getHitNumber()
        prot_acc = prot.getAccession()
        if self.res_file.getMascotVer() >= 2.3:
            db_idx = prot.getDB()
            prot_db = '%s::%s' % (db_idx, self.params.getDB(db_idx))
            prot_desc = self.pep_summary.getProteinDescription(prot_acc, db_idx)
            prot_mass = round(self.pep_summary.getProteinMass(prot_acc, db_idx))
        else:
            prot_desc = self.pep_summary.getProteinDescription(prot_acc)
            prot_mass = round(self.pep_summary.getProteinMass(prot_acc))
            prot_db = ''
            
        prot_matches = prot.getNumDisplayPeptides()
        prot_score = round(prot.getScore())


        return [prot_hit_num, prot_db, prot_acc, prot_desc,
                prot_mass, prot_matches, prot_score]



    def output_peptides(self, prot):
        def sort_mods(x,y):
            pos_re = re.compile('[A-Z]{1,2}(\d+).', re.IGNORECASE)
            a = int(pos_re.match(x).groups()[0])
            b = int(pos_re.match(y).groups()[0])
            return cmp(a,b)        
        
        prot_info = self.output_protein(prot)

        for j in range(1, prot.getNumPeptides()+1):
            if prot.getPeptideDuplicate(j) == prot.DUPE_DuplicateSameQuery:
                continue

            pep_query = prot.getPeptideQuery(j)

            pep = self.pep_summary.getPeptide(pep_query, prot.getPeptideP(j))

            pep_seq = pep.getPeptideStr()

            if self.args['mascot_var_mods']:
                pep_var_mod = self.parse_mod_pos(pep_seq, pep.getVarModsStr(),
                                                 (pep_query, prot.getPeptideP(j)))
                if self.res_file.getMascotVer() >= 2.4:
                    if pep.getSummedModsStr():
                        pep_sum_mod = self.parse_mod_pos(pep_seq, pep.getSummedModsStr(),
                                                         (pep_query, prot.getPeptideP(j)))
                        if pep_var_mod and pep_sum_mod:
                            pep_var_mod = '; '.join(sorted(pep_var_mod.split('; ') + pep_sum_mod.split('; '), cmp = sort_mods ))
                        else:
                            pep_var_mod = pep_var_mod or pep_sum_mod
                            
            else:
                pep_var_mod = self.pep_summary.getReadableVarMods(pep_query, pep.getRank())
            pep_exp_mz = pep.getObserved()
            pep_exp_z = pep.getCharge()
            pep_calc_mr = pep.getMrCalc()
            pep_delta = pep.getDelta()
            pep_score = pep.getIonsScore()
            pep_rank = pep.getRank() # can't use pretty-rank
            pep_start = prot.getPeptideStart(j)
            pep_end = prot.getPeptideEnd(j)
            pep_res_before = prot.getPeptideResidueBefore(j)
            pep_res_after = prot.getPeptideResidueAfter(j)
            pep_miss = pep.getMissedCleavages()
            pep_scan_title = unquote(self.res_file.getQuerySectionValueStr(pep_query, "Title"))
            
            #pep_quant = pep.getComponentStr() # Not currently used.

            # Sequence corresponds to mzReport.default_columns.
            yield (prot_info
                   + [pep_seq,
                      pep_var_mod,
                      pep_exp_mz,
                      pep_exp_z,
                      pep_calc_mr,
                      pep_delta,
                      pep_score,
                      pep_rank,
                      pep_start,
                      pep_end,
                      pep_res_before,
                      pep_res_after,
                      pep_miss,
                      pep_scan_title,
                      pep_query])

    def protein_report(self):
        for prot_hit_num in range(1, self.pep_summary.getNumberOfHits() + 1):
            prot = self.pep_summary.getHit(prot_hit_num)

            yield self.output_protein(prot)

            if self.args.get('show_same_set', default_show_same_set):
                i = 1
                while self.pep_summary.getNextSimilarProtein(prot_hit_num, i):
                    prot = self.pep_summary.getNextSimilarProtein(prot_hit_num, i)
                    yield self.output_protein(prot)
                    i += 1

            if self.args.get('show_sub_set', default_show_sub_set):
                i = 1
                while self.pep_summary.getNextSubsetProtein(prot_hit_num, i):
                    prot = self.pep_summary.getNextSubsetProtein(prot_hit_num, i)
                    yield self.output_protein(prot)
                    i += 1

    def peptide_report(self):
        for prot_hit_num in range(1, self.pep_summary.getNumberOfHits() + 1):
            prot = self.pep_summary.getHit(prot_hit_num)

            for pep in self.output_peptides(prot):
                yield pep

            if self.args.get('show_same_set', default_show_same_set):
                i = 1
                while self.pep_summary.getNextSimilarProtein(prot_hit_num, i):
                    prot = self.pep_summary.getNextSimilarProtein(prot_hit_num, i)
                    for pep in self.output_peptides(prot):
                        yield pep
                    i += 1

            if self.args.get('show_sub_set', default_show_sub_set):
                i = 1
                while self.pep_summary.getNextSubsetProtein(prot_hit_num, i):
                    prot = self.pep_summary.getNextSubsetProtein(prot_hit_num, i)
                    for pep in self.output_peptides(prot):
                        yield pep
                    i += 1
                    
                    
    def hasDecoyHits(self):
        # Depends on the documentation-contradicting behavior that a valid
        # section always returns 0 from getNumHits, and the undocumented behavior
        # that the decoy section is invalid in non-decoy searches.  Beware!!
        #return self.res_file.getNumHits(ms.ms_mascotresfile.SEC_DECOYSUMMARY) == 0
        #return self.pep_summary.getNumDecoyHitsAboveHomology(1)
        
        return self.params.getDECOY()
        



#def unescape(text):
    #def fixup(m):
        #text = m.group(0)
        #if text[:2] == "&#":
            ## character reference
            #try:
                #if text[:3] == "&#x":
                    #return chr(int(text[3:-1], 16))
                #else:
                    #return chr(int(text[2:-1]))
            #except ValueError:
                #pass
        #else:
            ## named entity
            #try:
                #text = chr(html.entities.name2codepoint[text[1:-1]])
            #except KeyError:
                #pass
        #return text # leave as is
    #return re.sub("&#?\w+;", fixup, text)


