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

import cStringIO
import csv
import htmlentitydefs
import os
import re
import tempfile
import time
import webbrowser

import urllib2

httpPackage = None
#try:
    #import pycurl
    #httpPackage = 'pycurl'
#except ImportError:
    #print 'PyCURL not detected.'
    
#try:
    #import requests
    #httpPackage = 'requests'
#except ImportError:
    #print 'Requests package not detected.'

import requests.packages.urllib3 # For the executable, due to py2exe silliness.


import requests

    
#import multiplierz.msparser as ms

from collections import defaultdict
from urllib import urlencode, unquote

import multiplierz.mass_biochem as mzFunctions
import multiplierz.mzReport as mzReport
#from multiplierz.mascot.report import MascotReport

from multiplierz import myData, logger_message
from multiplierz.settings import settings



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
                                         'pep_res_after', 'pep_var_mod', 'StringTitle', 'qualifiers',
                                         'TotalIonsIntensity', 'NumVals', 'StringIons1', 'StringIons2',
                                         'StringIons3', 'pep_scan_title'),
                                        ('acc', 'prot_desc', 'prot_score', 'prot_mass', 'prot_matches',
                                         'prot_rank', 'pep_query', 'mz', 'mr', 'charge', 'pred_mr',
                                         'delta', 'start', 'end', 'miss', 'pep_score', 'Peptide EValue',
                                         'pep_rank', 'res_before', 'seq', 'res_after', 'var_mods',
                                         'spec_desc', 'Qualifiers', 'Total Ion Intensity', 'Ion Count',
                                         'ion_str', 'Ion String 2', 'Ion String 3', 'spec_desc')))


class CurlForm(list):
    def add_formfield(self, name, field):
        self += (name, (pycurl.FORM_CONTENTS, field)),

    def add_file(self, name, filename):
        self += (name, (pycurl.FORM_FILE, filename)),

filetypes = {"MascotDAT" : 'dat',
             'CSV' : 'csv',
             'MZIdentML' : 'mzid'}


# From http://pymotw.com/2/urllib2/
# Allows uploading of files via urllib2.
# Remove if switch is made to Requests or some other more civilized library.
import itertools
import mimetools
import mimetypes
from cStringIO import StringIO
import urllib
class MultiPartForm(object):
    """Accumulate the data to be used when posting a form."""

    def __init__(self):
        self.form_fields = []
        self.files = []
        self.boundary = mimetools.choose_boundary()
        return
    
    def get_content_type(self):
        return 'multipart/form-data; boundary=%s' % self.boundary

    def add_field(self, name, value):
        """Add a simple field to the form data."""
        self.form_fields.append((name, value))
        return

    def add_file(self, fieldname, filename, fileHandle, mimetype=None):
        """Add a file to be uploaded."""
        body = fileHandle.read()
        if mimetype is None:
            mimetype = mimetypes.guess_type(filename)[0] or 'application/octet-stream'
        self.files.append((fieldname, filename, mimetype, body))
        return
    
    def __str__(self):
        """Return a string representing the form data, including attached files."""
        # Build a list of lists, each containing "lines" of the
        # request.  Each part is separated by a boundary string.
        # Once the list is built, return a string where each
        # line is separated by '\r\n'.  
        parts = []
        part_boundary = '--' + self.boundary
        
        # Add the form fields
        parts.extend(
            [ part_boundary,
              'Content-Disposition: form-data; name="%s"' % name,
              '',
              value,
            ]
            for name, value in self.form_fields
            )
        
        # Add the files to upload
        parts.extend(
            [ part_boundary,
              'Content-Disposition: file; name="%s"; filename="%s"' % \
                 (field_name, filename),
              'Content-Type: %s' % content_type,
              '',
              body,
            ]
            for field_name, filename, content_type, body in self.files
            )
        
        # Flatten the list and add closing boundary marker,
        # then return CR+LF separated data
        flattened = list(itertools.chain(*parts))
        flattened.append('--' + self.boundary + '--')
        flattened.append('')
        return '\r\n'.join(flattened)





class mascot(object):
    def __init__(self, server=None, version='2.1', verbose=False):
        '''Initialize a Mascot search session. This starts up Curl, but does not log in
        to the server.
        '''
        self.server = str(server) if server else server
        self.version = str(version) # So that James' comparison-as-str works correctly.
        self.verbose = verbose

        self.logged_in = False

        # A string with the name and path of an appropriate temp file
        # (varies by platform)
        fd, self.cookie_file_name = tempfile.mkstemp(text=True)
        os.close(fd)

        # Buffer for output
        self.output = cStringIO.StringIO()
        
        self.loginToken = ''

    def login(self, login, password):
        """
        Logs in to mascot. Returns "success" or "error" based on login success
        """

        login_url = self.server.strip(r'\/') + r'//cgi/login.pl'
        login_form = CurlForm()

        # URL encodings of form data
        login_form_seq = [('username', login),
                          ('password', password),
                          ('submit', 'Login'),
                          ('display', 'logout_prompt'),
                          ('savecookie', '1'),
                          ('action', 'login'),
                          ('userid', ''),
                          ('onerrdisplay', 'login_prompt')]
                                   
        req = urllib2.Request(login_url, 
                              data = urlencode(login_form_seq))
        f = urllib2.urlopen(req)
        try:
            self.loginToken = f.info()['Set-Cookie']
            # I suppose in this case precise disposition of the output is less important.
        except KeyError:
            raise RuntimeError, 'Login was not successful; check your username and password.'
        
        
        err = None
        login_re = re.compile(r'<B><FONT COLOR=#FF0000>Error:(.*)</FONT></B>')
        for line in f.read().splitlines():
            self.output.write(line)
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
            raise IOError, "Not logged in to Mascot server.  (Call .login() first!)"

        mascot_id = re.sub(r'^0+', '', str(mascot_id), count=1)

        date_url = self.server + r'/x-cgi/ms-review.exe?'
        
        req = urllib2.Request(date_url + urlencode([('CalledFromForm', '1'),
                                                    ('f0', mascot_id)]))
        req.add_header("Cookie", self.loginToken)
        f = urllib2.urlopen(req)
        #response.write(f.read())

        

        date_re = re.compile(r'.+/data/(\d+)/F%s.+' % mascot_id.zfill(6))

        matches = []
        #for line in response.getvalue().splitlines():
        for line in f.read().splitlines():
            m = date_re.match(line)
            if m:
                matches.append(m)
                
                
        if matches:
            return matches[-1].group(1)
        else:
            raise RuntimeError, "Search not found in date lookup."
        
        
    def download_dat(self, chosen_folder, mascot_id, date=None):
        mascot_id = str(mascot_id).zfill(6)
        dat_file = 'F%s.dat' % mascot_id

        dat_path = os.path.abspath(os.path.join(chosen_folder, dat_file))
        #if os.path.exists(dat_path):
            #return dat_path

        if not date:
            date = self.get_date(mascot_id)

        dat_url = self.server + r'/x-cgi/ms-status.exe?'

        dat_dict = {'Autorefresh': 'false',
                    'Show': 'RESULTFILE',
                    'DateDir': date,
                    'ResJob': dat_file}
        # Adding a Percolator parameter here breaks things!
        
        req = urllib2.Request(dat_url + urlencode(dat_dict))
        req.add_header("Cookie", self.loginToken)
        f = urllib2.urlopen(req)
        
        datOut = open(dat_path, 'w')
        datOut.write(f.read())
        datOut.close()

        return dat_path       

    def download_csv(self, mascot_id, **etcetc):
        return self.download_file(mascot_id, filetype = 'CSV', **etcetc)



    def download_file(self, mascot_id, filetype = 'CSV',
                     max_hits=1000, ion_cutoff=20, bold_red=True,
                     unassigned_queries=False, show_query_data=True,
                     show_same_set=False, show_sub_set=False, date=None,
                     ion_list=False, quant=False,
                     save_file=None, **kwargs):
        """Downloads mascot csv file for a given id and parameters.
        Returns csv as string.

        """

        extension = filetypes[filetype]

        if not self.logged_in:
            return

        mascot_id = str(mascot_id).zfill(6)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            return

        if self.version == '2.1':
            download_url = self.server + r'/cgi/export_dat.pl?'
        elif self.version >= '2.2': # not sure about 2.3 here
            download_url = self.server + r'/cgi/export_dat_2.pl?'

        download_form_seq = [('file', r'../data/%s/F%s.dat' % (date, mascot_id)),
                             ('export_format', filetype),
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

        if self.version == '2.1':
            download_form_seq.extend([('_mudpit', '99999999'),
                                      ('show_queries', int(show_query_data))])
        elif self.version >= '2.2': # not sure about 2.3 here...should be the same?
            download_form_seq.extend([('_server_mudpit_switch', '99999999'),
                                      ('pep_scan_title', int(show_query_data)),
                                      ('query_master', int(ion_list)),
                                      ('pep_quant', int(quant)),
                                      ('percolate', '0')])

        #download_url += urlencode(download_form_seq)
        
        if not save_file:
            save_file = os.path.join(os.getcwd(), "F%s.%s" % (mascot_id, extension))


        #self.crl.setopt(pycurl.HTTPGET, True)
        #self.crl.setopt(pycurl.URL, download_url)

        #with open(save_file, "w") as f:
            #self.crl.setopt(pycurl.WRITEFUNCTION, f.write)
            #self.crl.perform()

        #self.crl.setopt(pycurl.WRITEFUNCTION, self.output.write)
        
        req = urllib2.Request(download_url, data = urlencode(download_form_seq))
        req.add_header("Cookie", self.loginToken)
        f = urllib2.urlopen(req)
        
        with open(save_file, 'w') as output:
            output.write(f.read())

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
            raise IOError, "Not logged into to Mascot server."

        mascot_id = str(mascot_id).zfill(6)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            raise IOError, "Couldn't locate file on Mascot server."

        if self.version == '2.1':
            download_url = self.server + r'/cgi/export_dat.pl?'
        elif self.version >= '2.2': # not sure about 2.3 here
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

        if self.version == '2.1':
            download_form_seq.extend([('_mudpit', '99999999'),
                                      ('show_queries', int(show_query_data))])
        elif self.version >= '2.2': # not sure about 2.3 here...should be the same?
            download_form_seq.extend([('_server_mudpit_switch', '99999999'),
                                      ('pep_scan_title', int(show_query_data)),
                                      ('query_master', int(ion_list)),
                                      ('pep_quant', int(quant))])

        #download_url += urlencode(download_form_seq)

        #self.crl.setopt(pycurl.HTTPGET, True)
        #self.crl.setopt(pycurl.URL, download_url)
        req = urllib2.Request(download_url, data = urlencode(download_form_seq))
        req.add_header("Cookie", self.loginToken)
        mzidData = urllib2.urlopen(req)

        if not save_file:
            save_file = os.path.join(os.getcwd(), "F%s.mzid" % mascot_id)

        with open(save_file, "w") as f:
            #self.crl.setopt(pycurl.WRITEFUNCTION, f.write)
            #self.crl.perform()
            f.write(mzidData.read())
        mzidData.close()

        #self.crl.setopt(pycurl.WRITEFUNCTION, self.output.write)

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

        if self.version == '2.1':
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
        elif self.version >= '2.2': # 2.3 should be the same here, I think
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
        headers = reader.next()

        if mascot_var_mods and self.version >= "2.2":
            while (len(headers) == 0 or (headers[0] != "Identifier" and headers[0] != 'prot_hit_num')):
                headers = reader.next()
            if headers[0] != 'prot_hit_num':
                var_mod = reader.next()
                var_mod_dict = {}
                while len(var_mod) > 0:
                    var_mod_dict[var_mod[0]] = var_mod[1].split('(')[0][:-1]
                    var_mod = reader.next()
                headers = var_mod
            else:
                var_mod_dict = None

        while len(headers) == 0 or headers[0] != "prot_hit_num":
            headers = reader.next()

        #potentially missing protein related entries which may need to be duplicated across rows...
        missing = 10
        if self.version >= '2.2': # presumably this should work for 2.3
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
                headers = reader.next()
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
        if mascot_var_mods and self.version >= '2.2' and var_mod_dict:
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

    def get_range(self, a):
        '''Takes a list of integers and outputs a list of (start,end) pairs.

        Example:
        > get_range([1, 2, 3, 5, 6, 7])
        [(1, 3), (5, 7)]
        '''

        if not a:
            return []

        a.sort()
        c = [[a[0],a[0]]]
        for i in a[1:]:
            if i > c[-1][-1] + 1:
                c.append([i,i])
            else:
                c[-1][-1] = i

        return c

    def get_coverage(self, acc, mascot_id=None, date=None, cutoff=None, db_idx=1):
        if not self.logged_in:
            return

        mascot_id = str(mascot_id).zfill(6)
        acc = str(acc)
        cutoff = str(cutoff) if cutoff else "20"
        db_idx = str(db_idx or 1)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            return

        target_url = self.server + "/cgi/protein_view.pl?"

        form_seq = [('file', '../data/%s/F%s.dat' % (date, mascot_id)),
                    ('hit', acc),
                    ('px', '1'),
                    ('_ignoreionsscorebelow', cutoff)]

        if self.version >= '2.3':
            form_seq.extend((('db_idx', str(db_idx)),
                             ('_msresflags2', '4')))

        #target_url += urlencode(form_seq)

        #self.crl.setopt(pycurl.HTTPGET, True)
        #self.crl.setopt(pycurl.URL, target_url)

        #response = cStringIO.StringIO()
        #self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

        #self.crl.perform()

        #(coverage, start_end, cov_per) = self.parse_protein_view(response.getvalue())
        
        req = urllib2.Request(target_url, data = urlencode(form_seq))
        req.add_header("Cookie", self.loginToken)
        formResponse = urllib2.urlopen(req)
        
        (coverage, start_end, cov_per) = self.parse_protein_view(formResponse.read())

        #self.crl.setopt(pycurl.WRITEFUNCTION, self.output.write)

        # complex list comprehension to reformat the protein sequence in a
        # nice 10-AA-block way.
        cov_f = '\n'.join('%5d %s' % (i+1, ' '.join(coverage[i+j:i+j+10]
                                                    for j in range(0, 50, 10)))
                          for i in range(0, len(coverage), 50))

        # convert protein-index into text-index
        conv_x = lambda i: i + 6 + ((i - 1) / 50) * 11 + ((i - 1) % 50) / 10

        start_end = self.get_range([conv_x(i) for s,e in start_end for i in range(s,e+1)])

        return (cov_f, start_end, cov_per)

    def get_description(self, acc, db_idx=1, mascot_id=None, date=None):
        if not self.logged_in:
            # This fails when security is not enabled!
            return None, None

        mascot_id = str(mascot_id).zfill(6)
        acc = str(acc)
        db_idx = str(db_idx or 1)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            # There should be a warning set here!
            return None, None
        
        target_url = self.server + "/cgi/protein_view.pl?"

        form_seq = [('file', '../data/%s/F%s.dat' % (date, mascot_id)),
                    ('hit', acc),
                    ('px', '1')]

        if self.version >= '2.3':
            form_seq.append(('db_idx', db_idx))

        #target_url += urlencode(form_seq)

        #self.crl.setopt(pycurl.HTTPGET, True)
        #self.crl.setopt(pycurl.URL, target_url)

        #response = cStringIO.StringIO()
        #self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

        #self.crl.perform()

        #self.crl.setopt(pycurl.WRITEFUNCTION, self.output.write)
        
        req = urllib2.Request(target_url, data = urlencode(form_seq))
        req.add_header("Cookie", self.loginToken)
        f = urllib2.urlopen(req)

        # escape the pipes in the accession for regex processing
        acc_2 = '\|'.join(acc.split('|'))

        # multiline regex search is hopefully general enough
        desc_m = re.search((r'Match to: <B>%s</B> Score: <B>\d+</B>.+?<B>(.+?)</B>.+?'
                             'Nominal mass \(M<SUB>r</SUB>\): <B>(\d+)</B>') % acc_2,
                            #response.getvalue(),
                            f.read(),
                            flags=re.I|re.S)

        if not desc_m:
            return None, None
        else:
            return (str(desc_m.group(1)), int(desc_m.group(2)))

    def parse_protein_view(self, mascot_web_file):
        # if we're going to use regex to parse HTML we might as well go all out...

        cov_per_re = re.compile("Sequence Coverage: <B>(.+)%</B>")
        #try:
        cov_per = float(cov_per_re.search(mascot_web_file).group(1))
        #except:
            #cov_per = 0.0

        seq_re = re.compile((r'Matched peptides shown in <B><FONT COLOR=#FF0000>'
                            'Bold Red</FONT></B>(.*?)</PRE></FONT>'), flags=re.DOTALL)
        prot_re = re.compile(r'<.*?>|\s|\d')

        #try:
        coverage = prot_re.sub('', seq_re.search(mascot_web_file).group(1))
        #except:
            #coverage = ""

        #try:
        start_end = set()
        for m in re.finditer(r'<B><FONT COLOR=#FF0000>\s+(\d+) \- (\d+)', mascot_web_file):
            start_end.add((int(m.group(1)), int(m.group(2))))
        #except:
            #start_end = set()

        return (coverage, sorted(start_end), cov_per)

    def get_msms(self, query, img_file, mascot_id=None, score='',
                 pep_rank=1, date='', sequence='', mod_only=False, dat_file=True):
        if not self.logged_in:
            return

        mascot_id = str(mascot_id).zfill(6)

        if not date:
            date = self.get_date(mascot_id)

        if not date:
            return

        query = str(query)
        pep_rank = str(pep_rank)

        #Force PepRank = 1 (Hit) and then go down till sequence and score match
        pep_rank = '1'

        target_url = self.server + "/cgi/peptide_view.pl?"

        form_dict = {'file': '../data/%s/F%s.dat' % (date, mascot_id),
                    'query': query,
                    'hit': pep_rank,
                    'px': '1'}

        #self.crl.setopt(pycurl.HTTPGET, True)
        #self.crl.setopt(pycurl.URL, target_url + urlencode(form_dict.items()))

        #response = cStringIO.StringIO()
        #self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

        #self.crl.perform()
        
        req = urllib2.Request(target_url, data = urlencode(form_dict))
        req.add_header("Cookie", self.loginToken)
        f = urllib2.urlopen(req)

        #mascot_web_file = response.getvalue().splitlines()
        mascot_web_file = f.read().splitlines()

        (peptide,ionInfo,phosphoPos,
         oxidPos,img,text,info) = self.parse_msms_view(mascot_web_file, sequence, score, mod_only)

        while int(pep_rank) < 10 and info == 'Fail':
            pep_rank = str(int(pep_rank)+1)

            target_url = self.server+ '/cgi/peptide_view.pl?'

            form_dict['hit'] = pep_rank

            #self.crl.setopt(pycurl.URL, target_url + urlencode(form_dict.items()))

            #response = cStringIO.StringIO()
            #self.crl.setopt(pycurl.WRITEFUNCTION, response.write)

            #self.crl.perform()
            
            req = urllib2.Request(target_url, data = urlencode(form_dict))
            req.add_header("Cookie", self.loginToken)
            f = urllib2.urlopen(req)

            #mascot_web_file = response.getvalue().splitlines()
            mascot_web_file = f.read().splitlines()

            (peptide,ionInfo,phosphoPos,
             oxidPos,img,text,info) = self.parse_msms_view(mascot_web_file, sequence, score, mod_only)

        if mod_only:
            return (text,ionInfo,phosphoPos,oxidPos)

        target_url = self.server + '/cgi/msms_gif.pl?msms_gif_' + img

        with open(img_file, mode='wb') as fh:
            #self.crl.setopt(pycurl.URL, target_url)
            #self.crl.setopt(pycurl.WRITEFUNCTION, fh.write)
            #self.crl.perform()
            req = urllib2.Request(target_url)
            req.add_header("Cookie", self.loginToken)
            fh.write(urllib2.urlopen(req).read())
            # Honestly, I'd be surprised if this works, or if someone notices it doesn't.

        return (text,ionInfo,phosphoPos,oxidPos)

    def parse_msms_view(self, mascot_web_file, sequence, score, mod_only):
        img = ''
        text = ''
        modMatch = 0
        peptide = ''
        tableMatch = 0
        bindex = 0
        b2index = 0
        cindex = 0
        c2index = 0
        yindex = 0
        y2index = 0
        z_plus1_index = 0
        z_plus1_2index = 0
        ionindex = 0
        massesMatch = 0
        ionInfo = []
        ionArray = []
        phosphoPos = []
        oxidPos = []
        ionScore = 0

        imgRE = re.compile('SRC=\"\./msms_gif.pl\?msms_gif_(.+)\"\>')
        modRE = re.compile('\<B\>Variable modifications\:\s\<\/B\>')
        ionRE = re.compile('\<B\>Ions Score\:(.+)Expect.+')
        pepRE = re.compile('MS\/MS Fragmentation of \<B\>\<FONT COLOR\=\#FF0000\>(.+)\<\/FONT\>.+')
        masstableRE = re.compile('\s+\<TR BGCOLOR\=\#cccccc\>')
        massesMatchRE = re.compile('\s+\<TR ALIGN\=\"RIGHT\"\>')

        html_re = re.compile('<.+?>')

        for line in mascot_web_file:
            igRE = imgRE.match(line)
            mRE = modRE.match(line)
            inRE = ionRE.match(line)
            pRE = pepRE.match(line)
            mtRE = masstableRE.match(line)
            mmRE = massesMatchRE.match(line)
            if pRE:
                peptide = pRE.group(1)
                if peptide != sequence:
                    return (peptide,ionInfo,phosphoPos,oxidPos,img,text, 'Fail')

            if igRE:
                img = igRE.group(1)

            if mRE:
                modMatch = 1
                continue

            if modMatch == 1 and not inRE:
                mod = html_re.sub('', mod)
                mod = re.sub('\(shown in table\)', '', mod)
                mod = ' '.join(mod.split())

                phosphoRE = re.match('(.)(\d+).*Phospho.*', mod)
                oxidRE = re.match('(.)(\d+).*Oxidation.*', mod)
                genModRE = re.match('(.+)\s+\:\s*(.+)\s*', mod)
                if phosphoRE:
                    phosphoGroups = phosphoRE.groups()
                    phosphoPos.append(int(phosphoGroups[1])-1)
                    aa = str(phosphoGroups[0])
                    text+= str(aa)+str(phosphoGroups[1]) + ': Phospho'
                    if re.match('.*eutral.*',mod):
                        text += ' -NL'
                    text += '; '
                elif oxidRE:
                    oxidGroups = oxidRE.groups()
                    oxidPos.append(int(oxidGroups[1])-1)
                    aa = str(oxidGroups[0])
                    text += str(aa)+str(oxidGroups[1]) + ': Oxidation; '
                elif genModRE:
                    text += genModRE.group(1) + ': ' + genModRE.group(2) + '; '

            if inRE:
                modMatch = 0
                ionScore = inRE.group(1)
                ionScore = re.sub('\s+','',ionScore)
                ionScore = re.sub('\<B\>','',ionScore)
                ionScore = re.sub('\</B\>','',ionScore)
                ionScore = int(round(float(ionScore)))
                decRem = float(score) - int(float(score))
                score = int(round(float(score)))
                if decRem >= 0.5 and decRem < 0.51:
                    if abs(ionScore - score) > 1:
                        return (peptide,ionInfo,phosphoPos,oxidPos,img,text, 'Fail')
                else:
                    if ionScore != score:
                        return (peptide,ionInfo,phosphoPos,oxidPos,img,text, 'Fail')

            if mtRE:
                tableMatch = 1
                continue

            if tableMatch == 1 and not mmRE:
                if mod_only:
                    return (peptide,ionInfo,phosphoPos,oxidPos,img,text, 'Success')

                thRE = re.match('\s+<TH>(.+)</TH>',line)
                if thRE:
                    ionindex+=1
                    thGroups = thRE.groups()
                    ion = thGroups[0]
                    ion = re.sub('\<SUP\>','',ion)
                    ion = re.sub('\</SUP\>','',ion)

                    if ion == 'b':
                        bindex = ionindex
                    elif ion == 'b++':
                        b2index = ionindex
                    if ion == 'c':
                        cindex = ionindex
                    elif ion == 'c++':
                        c2index = ionindex
                    elif ion == 'y':
                        yindex = ionindex
                    elif ion == 'y++':
                        y2index = ionindex
                    elif  ion == 'z+1':
                        z_plus1_index = ionindex
                    elif ion == 'z+1++':
                        z_plus1_2index = ionindex

            if mmRE:
                tableMatch = 0
                massesMatch = 1
                ionindex = 0
                ionInfo.append(ionArray)
                ionArray = []
                continue

            if massesMatch == 1 and not re.match('\s+\<\/TR\>',line):
                tdRE = re.match('\s+<TD(.+)</TD>',line)
                if tdRE:
                    ionindex+=1
                    tdGroups = tdRE.groups()
                    tdinfo = tdGroups[0]
                    boldRedRE = re.match('.+\#FF0000.+',line)
                    if boldRedRE:
                        if ionindex == bindex:
                            ionArray.append('b')
                        elif ionindex == b2index:
                            ionArray.append('b++')
                        if ionindex == cindex:
                            ionArray.append('c')
                        elif ionindex == c2index:
                            ionArray.append('c++')
                        elif ionindex == yindex:
                            ionArray.append('y')
                        elif ionindex == y2index:
                            ionArray.append('y++')
                        elif ionindex == z_plus1_index:
                            ionArray.append('z+1')
                        elif ionindex == z_plus1_2index:
                            ionArray.append('z+1++')


            if massesMatch == 1 and re.match('\<\/TABLE\>',line):
                massesMatch = 0

        ionInfo.append(ionArray)
        ionInfo = ionInfo[1:]

        return (peptide,ionInfo,phosphoPos,oxidPos,img,text, 'Success')

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

    def close(self, out=None):
        '''Close this Mascot search session.

        'out' is a file-like object which will be called on the output buffer--useful
        for debugging the response
        '''
        #self.crl.close()

        if out and self.verbose:
            out.write(self.output.getvalue())
        elif self.verbose:
            print self.output.getvalue()

        self.output.close()
        os.unlink(self.cookie_file_name)


class FieldParser(dict):
    def __init__(self, output, version):
        '''Initializes the parser. Takes a reference to a file-like object to
        write the output to'''
        self.output = output
        self.version = version

        self.current_field = None
        self.selected = {}
        self['MODS'] = []
        self['ALLMODS'] = []
        self.regexes = [re.compile(r'SELECT (?:.* )?NAME="(\w+)"', re.I),
                        re.compile(r'OPTION>(.+)</OPTION', re.I),
                        re.compile(r'OPTION SELECTED>(.+)</OPTION', re.I),
                        re.compile(r'</SELECT>', re.I),
                        re.compile(r'(?<!all)[M,m]ods\[\d+\] = [\',"](.+)[\',"];', re.I),
                        re.compile(r'allMods\[\d+\] = [\',"](.+)[\',"];')]
        

    def __call__(self,s):
        '''The __call__ method stores the string in an internal buffer, and also
        writes it to the output buffer. It has to use an internal buffer because
        it can't guarantee that the input will end on a newline--Curl breaks the
        page into arbitrary chunks'''
        self.output.write(s)

    def parse_fields(self):
        '''Parses the buffer and builds up the dictionaries of options and defaults'''
        for line in self.output.getvalue().splitlines():
            self.parse_line(line)

        if 'IT_MODS' not in self:
            self['IT_MODS'] = self['MODS']

        assert any(self.values()), "Failed to get Mascot parameter fields.  Are your Mascot server URL and security settings correctly set?"

        return self

    def parse_line(self,s):
        '''Parses a single line of the buffer, looking for something that looks
        like a list of options'''
        if not self.current_field:
            m0 = self.regexes[0].search(s) # search for beginning of new param list
            if self.version < '2.3' and 'hiddenMods' in s:
                try:
                    self['ALLMODS'].append(unescape(s.split('"')[1]))            
                except IndexError:
                    pass
            elif m0:
                if m0.group(1) not in ('MODS', 'IT_MODS') or self.version < '2.3':
                    self.current_field = m0.group(1)
                    self[self.current_field] = []
            elif self.version >= '2.3':
                m4 = self.regexes[4].search(s) # search for Mascot 2.3-style mod
                if m4:
                    self['MODS'].append(unescape(m4.group(1)))
                elif self.regexes[5].search(s): # Mods not in the 'hidden' list.
                    self['ALLMODS'].append(unescape(self.regexes[5].search(s).group(1)))
                else:
                    pass #NotAllMods
        else:
            m1 = self.regexes[1].search(s) # search for an option
            m2 = self.regexes[2].search(s) # search for the default option
            m3 = self.regexes[3].search(s) # search for the end of the list
            if m1:
                self[self.current_field].append(unescape(m1.group(1)))
            elif m2:
                self[self.current_field].append(unescape(m2.group(1)))
                self.selected[self.current_field] = unescape(m2.group(1))
            elif m3:
                self.current_field = None
                self.parse_line(s[m3.end():])
            elif 'OPTION VALUE' in s: # Regexes are silly.
                value = s.split('>')[1].split('<')[0]
                self[self.current_field].append(value)
                



class MascotSearcher(object):
    def __init__(self, server = None, version = None, verbose = True,
                 open_tabs = False):
        
        if server:
            self.server = str(server)
            self.version = version
        else:
            self.server = settings.mascot_server
            self.version = settings.mascot_version
            
        self.verbose = verbose
        self.open_tabs = open_tabs
        
        self.cookies = {}
        self.output = StringIO()
    
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
        
        
        err = None
        login_re = re.compile(r'<B><FONT COLOR=#FF0000>Error:(.*)</FONT></B>')
        for line in r.text.splitlines():
            self.output.write(line)
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
        #return self.get_fields_better()
        search_form_url = self.server + r'/cgi/search_form.pl'
        field_parser = FieldParser(self.output, version)
        search_form = dict([('FORMVER','1.01'), ('SEARCH','MIS')])
        
        r = requests.get(search_form_url + "?" + urlencode(search_form),
                         cookies = self.cookies)
        
        field_parser(r.text)
        
        return field_parser.parse_fields()
          
    
    def get_fields_better(self):
        # Doesn't seem to return a bunch of fields that RC_MascotSearch wants.
        
        class fieldParser():
            def __init__(self):
                self.content = cStringIO.StringIO()
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
            fieldDict['PFA'] = range(0, 9)
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
            print '\n'.join(rec)
                
        return (dat_file_id,err)     
    
    def parse_response(self, s):
        self.output.write(s)
        m = re.search(r'A HREF="\.\.(/cgi/master_results.+/data/(\d+)/F(.+)\.dat)">Click here to see Search Report',s)
        m2 = re.search(r'.*\[M\d+\]',s)
        if m:
            if self.open_tabs:
                print "Opening web browser to result page."
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

        import multiplierz.msparser as ms

        self.dat_file_path = str(dat_file_path)
        self.res_file = ms.ms_mascotresfile(str(dat_file_path))

        if not self.res_file.isValid():
            print dat_file_path
            errors = []
            for i in range(1, self.res_file.getNumberOfErrors() + 1):
                errors.append(self.res_file.getErrorString(i))
                
            if any(['missing or corrupt headers' in x for x in errors]):
                errors.append("\n(Typically this occurs because Mascot couldn't find \n"
                              "the requested search and returned an invalid .DAT file.)")
            
            raise IOError, "MSParser errors: \n %s" % '\n'.join(errors)
        
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
            print "Experimental error-tolerant search results feature!"
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

    def getFILENAME_correctly(self):
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

        mascot_header.append(['Peak list data path', self.getFILENAME_correctly()])
        mascot_header.append(['Peak list format', self.params.getFORMAT()])
        mascot_header.append(['Search type', self.params.getSEARCH()])
        mascot_header.append(['Mascot version', self.res_file.getMascotVer()])
        if self.res_file.getMascotVer() >= '2.3':
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

        #mascot_header.append([' ', ' '])
        #mascot_header.append(['Mascot Decoy Search', 'True' if self.params.getDECOY() else 'False'])
        
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
        if self.res_file.getMascotVer() >= '2.3':
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
                if self.res_file.getMascotVer() >= '2.4':
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
        



def unescape(text):
    def fixup(m):
        text = m.group(0)
        if text[:2] == "&#":
            # character reference
            try:
                if text[:3] == "&#x":
                    return unichr(int(text[3:-1], 16))
                else:
                    return unichr(int(text[2:-1]))
            except ValueError:
                pass
        else:
            # named entity
            try:
                text = unichr(htmlentitydefs.name2codepoint[text[1:-1]])
            except KeyError:
                pass
        return text # leave as is
    return re.sub("&#?\w+;", fixup, text)


