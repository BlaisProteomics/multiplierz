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
import re

import multiplierz.mzReport
#import multiplierz.mzFunctions
from multiplierz.mass_biochem import mw

# This used to prefer lxml-based ElementTree; there *shouldn't* be differences?
import xml.etree.cElementTree as ET

def xTandemArray(xTandem_file, rev_mods=None, db_regex=None):
    """Parses XTandem output file and organizes the data as an array of dicts,
    where each dict contains the value for each column of the output sheet.
    """

    # Parse the xml file
    outtree = ET.parse(xTandem_file)

    # Quick scan of all the protein elements in xmltree in order to establish
    # protein rank
    prot_array = []
    name_array = []
    prot_hash = {}

    if not rev_mods:
        rev_mods = {}

    if db_regex:
        (acc_regex,desc_regex) = db_regex

    # started some regex stuff...the GUI has a regex menu, and we try that.
    # if it fails, we'll just stick everything in prot_desc

    # Make a list of tuples (protein score, protein name) organized by
    # decreasing protein score
    for protHit in outtree.iterfind('.//protein'):
        if db_regex:
            m1 = acc_regex.match(protHit.get('label'))
            m2 = desc_regex.match(protHit.get('label'))
            if m1:
                prot_acc = m1.group(1)
            else:
                prot_acc = ''
            if m2:
                prot_desc = m2.group(1)
            else:
                prot_desc = protHit.get('label')
        else:
            prot_acc = ''
            prot_desc = protHit.get('label')
        name_array.append(prot_acc or prot_desc)
        score = -float(protHit.get('expect'))
        prot_array.append((score, prot_acc or prot_desc))

    prot_array = list(set(prot_array)) # removes duplicates
    prot_array.sort(reverse=True)

    #Make a hash of protein ranks based on the position of the name in protList
    for k in range(len(prot_array)): # what is the purpose of this hash?
        prot_hash[prot_array[k][1]] = k + 1 # will have duplicates

    #Drill down layer by layer to get data needed for multiplierz
    x_array = []

    #Look at each query one at a time
    for group in outtree.iterfind('group[@type="model"]'):
        #Rank peptides in each query
        pep_array = []
        pep_hash = {}
        for dom in group.iterfind('.//domain'):
            pep_array.append((dom.get('hyperscore'), dom.get('seq')))
        pep_array = list(set(pep_array))
        pep_array.sort(reverse=True)
        pep_hash = dict(zip([p[1] for p in pep_array], range(1, len(pep_array)+1)))

        #General information for all hits in the group id
        pep_query = group.get('id')
        pep_exp_mh = float(group.get('mh'))
        pep_exp_z = int(group.get('z'))
        pep_exp_mz = str((pep_exp_mh + (pep_exp_z - 1)*1.00727638)/pep_exp_z)

        #Within each query, examine the 'fragment ion mass spectrum' group
        for subg in group.iterfind('group[@label="fragment ion mass spectrum"]'):
            StringTitle = subg.find('note[@label="Description"]').text

            #Xval = subg.find('.//*[@units="MASSTOCHARGERATIO"]').find('*').text.split()
            #Yval = subg.find('.//*[@units="UNKNOWN"]').find('*').text.split()

        subre = re.compile(r'\s+')

        # Look at each protein hit for the query
        for prot in group.iterfind('protein'):
            prot_seq = prot.find('peptide').text
            prot_seq = subre.sub('', prot_seq)
            prot_mass = int(mw(prot_seq))

            if db_regex:
                m1 = acc_regex.match(prot.get('label'))
                m2 = desc_regex.match(prot.get('label'))
                if m1:
                    prot_acc = m1.group(1)
                else:
                    prot_acc = ''
                if m2:
                    prot_desc = m2.group(1)
                else:
                    prot_desc = prot.get('label')
            else:
                prot_acc = ''
                prot_desc = prot.get('label')

            prot_score = -float(prot.get('expect'))

            prot_rank = prot_hash[prot_acc or prot_desc]
            prot_matches = name_array.count(prot_acc or prot_desc)

            #Look at the peptide that generated the protein hit
            for dom in prot.iterfind('*/domain'):
                pep_seq = dom.get('seq')
                pep_start = dom.get('start')

                # Incomprehensible list comprehension: rev_mods is a dictionary of spec:modnames
                # for the modifications available. The comprehension fetches the name if it exists
                # and creates a string listing all the mods for this hit
                # - example string: "C1: Carbamidomethyl (C); C8: -17.0265"
                var_mods = '; '.join('%s%d: %s' % (a.get('type'),
                                                   int(a.get('at')) - int(pep_start) + 1,
                                                   rev_mods.get('%.4f@%s' % (float(a.get('modified')),
                                                                             a.get('type')),
                                                               a.get('modified')))
                                     for a in dom.iterfind('aa[@at][@type][@modified]'))

                pep_calc_mr = float(dom.get('mh')) - 1.00727638
                pep_delta = dom.get('delta')
                pep_score = dom.get('hyperscore')
                pep_expect = dom.get('expect')
                pep_rank = pep_hash[pep_seq]
                pep_end = dom.get('end')
                before = dom.get('pre')[-1]
                after = dom.get('post')[0]
                try:
                    miss = dom.get('missed_cleavages')
                except KeyError:
                    miss = '0'

                #For each peptide domain, make a hash for one row in a multiplierz spreadsheet
                row = { #'Protein Rank': str(prot_rank),
                        'Accession Number': prot_acc,
                        #'Protein Description': prot_desc,
                        #'Protein Mass': str(prot_mass),
                        #'Protein Matches': prot_matches,
                        #'Protein Score': str(prot_score),
                        'Peptide Sequence': pep_seq,
                        'Variable Modifications': var_mods,
                        'Experimental mz': pep_exp_mz,
                        'Charge': str(pep_exp_z),
                        'Predicted mr': str(pep_calc_mr),
                        'Delta': pep_delta,
                        'Peptide Score': pep_score,
                        'Peptide Rank': str(pep_rank),
                        'Start Position': pep_start,
                        'End Position': pep_end,
                        'Preceding Residue': before,
                        'Following Residue': after,
                        'Missed Cleavages': miss,
                        'Spectrum Description': StringTitle,
                        'Query': pep_query,
                        'Expect': pep_expect }

                x_array.append(row)

    # we'll sort the array by rank before returning it
    x_array.sort(key = lambda x: int(x['Protein Rank']))

    return x_array, ['Expect',
                     'Peptide Rank',
                     'Peptide Score',
                     #'Protein Description',
                     #'Protein Matches',                  
                     #'Protein Rank',
                     #'Protein Score', 
                     #'Protein Mass',
                     'Experimental mz',
                     'Predicted mr',
                     'Delta',
                     'Accession Number',
                     'Missed Cleavages',
                     'Peptide Sequence',
                     'Variable Modifications',                  
                     'Charge',
                     'Spectrum Description',
                     'Preceding Residue',
                     'Following Residue',                  
                     'Start Position',
                     'End Position',                  
                     'Query']



headerSequence = ['Expect',
                  'Peptide Rank',
                  'Peptide Score',
                  #'Protein Description',
                  #'Protein Matches',                  
                  #'Protein Rank',
                  #'Protein Score', 
                  #'Protein Mass',
                  'Experimental mz',
                  'Predicted mr',
                  'Delta',
                  'Accession Number',
                  'Missed Cleavages',
                  'Peptide Sequence',
                  'Variable Modifications',                  
                  'Charge',
                  'Spectrum Description',
                  'Preceding Residue',
                  'Following Residue',                  
                  'Start Position',
                  'End Position',                  
                  'Query']

def xtandemParse(tandemfile):
    """
    PSM-based X Tandem results parser.
    """
    root = ET.parse(tandemfile).getroot()
    
    models = [x for x in root if x.get('type') == 'model']
    psms = []
    for model in models:
        spectrumEl = [x for x in model.iter('group') if x.get('label') == 'fragment ion mass spectrum'][0]
        description = spectrumEl.iter('note').next().text
        
        pep_query = model.get('id')
        pep_exp_mh = float(model.get('mh'))
        pep_exp_z = int(model.get('z'))        
        
        bestpep = None
        bestdom = None
        bestscore = -99
        for peptide in model.iter('peptide'):
            domain = peptide.iter('domain').next()
            hyperscore = float(domain.get('hyperscore'))
            if hyperscore > bestscore:
                bestpep = peptide
                bestdom = domain
        dom = bestdom
        
        sites = dom.iter('aa')
        mods = []
        for site in sites:
            loc = 1 + int(site.get('at')) - int(dom.get('start'))
            amino = site.get('type')
            mod = site.get('modified')
            mods.append('%s%d: %s' % (amino, loc, mod))
        modstring = '; '.join(mods)
        
        pep_seq = dom.get('seq')
        pep_start = dom.get('start')        
        pep_calc_mr = float(dom.get('mh')) - 1.00727638
        pep_delta = dom.get('delta')
        pep_score = dom.get('hyperscore')
        pep_expect = dom.get('expect')
        pep_rank = 1 # By construction, since we're taking the best pep for spectrum?
        pep_end = dom.get('end')
        before = dom.get('pre')[-1]
        after = dom.get('post')[0]  
        try:
            miss = dom.get('missed_cleavages')
        except KeyError:
            miss = '0'        
    
        row = { #'Protein Rank': None,
                'Accession Number': model.get('label'),
                #'Protein Description': None,
                #'Protein Mass': None,
                #'Protein Matches': None,
                #'Protein Score': None,
                'Peptide Sequence': pep_seq,
                'Variable Modifications': modstring,
                'Experimental mz': str((pep_exp_mh + (pep_exp_z - 1)*1.00727638)/pep_exp_z),
                'Charge': str(pep_exp_z),
                'Predicted mr': str(pep_calc_mr),
                'Delta': pep_delta,
                'Peptide Score': pep_score,
                'Peptide Rank': str(pep_rank),
                'Start Position': pep_start,
                'End Position': pep_end,
                'Preceding Residue': before,
                'Following Residue': after,
                'Missed Cleavages': miss,
                'Spectrum Description': description,
                'Query': pep_query,
                'Expect': pep_expect }        
        psms.append(row)
    
    return psms
    



def format_XML(xTandem_file, save_file, parameters = None,
               rev_mods=None, db_regex=None):
    '''Create a report from the xTandem output (which is an XML file)'''

    #rows, headers = xTandemArray(xTandem_file, rev_mods, db_regex)
    rows = xtandemParse(xTandem_file)

    if os.path.exists(save_file):
        try:
            os.remove(save_file)
        except WindowsError:
            pass

    report = multiplierz.mzReport.writer(save_file, columns=headerSequence)

    for row in rows:
        report.write(row)

    report.close()
    report = multiplierz.mzReport.writer(save_file, sheet_name = "XTandem_Header",
                                         columns = ["Category", "Setting", "Value"])
    if parameters:
        for category, settings in parameters.items():
            for setting, value in settings.items():
                report.write({'Category':category,
                              'Setting':setting,
                              'Value':value})
    else:
        report.write(['No', 'Parameters', 'Found'])
    report.close()

