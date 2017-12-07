import os, sys
import re
import time
import multiplierz.unimod as unimod
from multiplierz import myData
from multiplierz.mzReport import writer
from subprocess import call
from collections import defaultdict
from multiplierz.mgf import parse_mgf, extractor_name, standard_title_parse
from multiplierz.mzAPI import mzFile
from multiplierz import protonMass, __version__
import xml.etree.ElementTree as ET
import csv



__all__ = ['CometSearch']


parameterPreamble = ['# comet_version 2017.01 rev. 2\n',
                     '# Comet MS/MS search engine parameters file.\n',
                     "# Everything following the '#' symbol is treated as a comment.\n"]

# For instance, matches 'comet.2015012.win64.exe'
cometRegex = r'comet\.20[0-9]{5}\.win(64|32)\.exe'
cometExecutable = re.compile(cometRegex)

from multiplierz import SettingsFile, myData
import multiplierz.settings as settings
cometPath = settings.get_comet()


# An error will be raised on invoking a CometSearch object if this path is
# missing or invalid.


    
    
enzymeDict = {'Arg_C': ['1', 'R', 'P'],
              'Asp_N': ['0', 'D', '-'],
              'CNBr': ['1', 'M', '-'],
              'Chymotrypsin': ['1', 'FWYL', 'P'],
              'Glu_C': ['1', 'DE', 'P'],
              'Lys_C': ['1', 'K', 'P'],
              'Lys_N': ['0', 'K', '-'],
              'No_enzyme': ['0', '-', '-'],
              'PepsinA': ['1', 'FL', 'P'],
              'Trypsin': ['1', 'KR', 'P'],
              'Trypsin/P': ['1', 'KR', '-']}
unitsDict = {'amu':'0', 'ppm':'2'}


fixedModTypesForward = {'A': 'add_A_alanine',
                        'B': 'add_B_user_amino_acid',
                        'C': 'add_C_cysteine',
                        'D': 'add_D_aspartic_acid',
                        'E': 'add_E_glutamic_acid',
                        'F': 'add_F_phenylalanine',
                        'G': 'add_G_glycine',
                        'H': 'add_H_histidine',
                        'I': 'add_I_isoleucine',
                        'J': 'add_J_user_amino_acid',
                        'K': 'add_K_lysine',
                        'L': 'add_L_leucine',
                        'M': 'add_M_methionine',
                        'N': 'add_N_asparagine',
                        'O': 'add_O_ornithine',
                        'P': 'add_P_proline',
                        'Q': 'add_Q_glutamine',
                        'R': 'add_R_arginine',
                        'S': 'add_S_serine',
                        'T': 'add_T_threonine',
                        'U': 'add_U_user_amino_acid',
                        'V': 'add_V_valine',
                        'W': 'add_W_tryptophan',
                        'X': 'add_X_user_amino_acid',
                        'Y': 'add_Y_tyrosine',
                        'Z': 'add_Z_user_amino_acid'}
fixedModTypesBackward = dict([(x,y) for y, x in fixedModTypesForward.items()])


# Does not include varmods args.
defaultParameters = {'activation_method': 'ALL',
                     'add_A_alanine': '0.0000',
                     'add_B_user_amino_acid': '0.0000',
                     'add_C_cysteine': '57.021464',
                     'add_Cterm_peptide': '0.0',
                     'add_Cterm_protein': '0.0',
                     'add_D_aspartic_acid': '0.0000',
                     'add_E_glutamic_acid': '0.0000',
                     'add_F_phenylalanine': '0.0000',
                     'add_G_glycine': '0.0000',
                     'add_H_histidine': '0.0000',
                     'add_I_isoleucine': '0.0000',
                     'add_J_user_amino_acid': '0.0000',
                     'add_K_lysine': '0.0000',
                     'add_L_leucine': '0.0000',
                     'add_M_methionine': '0.0000',
                     'add_N_asparagine': '0.0000',
                     'add_Nterm_peptide': '0.0',
                     'add_Nterm_protein': '0.0',
                     'add_O_ornithine': '0.0000',
                     'add_P_proline': '0.0000',
                     'add_Q_glutamine': '0.0000',
                     'add_R_arginine': '0.0000',
                     'add_S_serine': '0.0000',
                     'add_T_threonine': '0.0000',
                     'add_U_user_amino_acid': '0.0000',
                     'add_V_valine': '0.0000',
                     'add_W_tryptophan': '0.0000',
                     'add_X_user_amino_acid': '0.0000',
                     'add_Y_tyrosine': '0.0000',
                     'add_Z_user_amino_acid': '0.0000',
                     'allowed_missed_cleavage': '2',
                     'clear_mz_range': '0.0 0.0',
                     'clip_nterm_methionine': '0',
                     'database_name': '/some/path/db.fasta',
                     'decoy_prefix': 'DECOY_',
                     'decoy_search': '0',
                     'digest_mass_range': '600.0 5000.0',
                     'fragment_bin_offset': '0.4',
                     'fragment_bin_tol': '1.0005',
                     'isotope_error': '0',
                     'mass_type_fragment': '1',
                     'mass_type_parent': '1',
                     'max_fragment_charge': '3',
                     'max_precursor_charge': '6',
                     'max_variable_mods_in_peptide': '5',
                     'minimum_intensity': '0',
                     'minimum_peaks': '10',
                     'ms_level': '2',
                     'nucleotide_reading_frame': '0',
                     'num_enzyme_termini': '2',
                     'num_output_lines': '5',
                     'num_results': '100',
                     'num_threads': '0',
                     'output_outfiles': '0',
                     'output_pepxmlfile': '1',
                     'output_percolatorfile': '0',
                     'output_sqtfile': '0',
                     'output_sqtstream': '0',
                     'output_suffix': '',
                     'output_txtfile': '0',
                     'override_charge': '0',
                     'peptide_mass_tolerance': '3.00',
                     'peptide_mass_units': '0',
                     'precursor_charge': '0 0',
                     'print_expect_score': '1',
                     'remove_precursor_peak': '0',
                     'remove_precursor_tolerance': '1.5',
                     'require_variable_mod': '0',
                     'sample_enzyme_number': '1',
                     'scan_range': '0 0',
                     'search_enzyme_number': '1',
                     'show_fragment_ions': '0',
                     'skip_researching': '1',
                     'spectrum_batch_size': '0',
                     'theoretical_fragment_ions': '1',
                     'use_A_ions': '0',
                     'use_B_ions': '1',
                     'use_C_ions': '0',
                     'use_NL_ions': '1',
                     'use_X_ions': '0',
                     'use_Y_ions': '1',
                     'use_Z_ions': '0',
                     'use_sparse_matrix': '1'}

# This needs '[COMET ENZYME INFO]' placed in front of it in the param file.
defaultEnzymes = {'0': {'active': 0, 'name': 'No_enzyme', 'p': '-', 'specificity': '-'},
                  '1': {'active': 1, 'name': 'Trypsin', 'p': 'P', 'specificity': 'KR'},
                  '2': {'active': 1, 'name': 'Trypsin/P', 'p': '-', 'specificity': 'KR'},
                  '3': {'active': 1, 'name': 'Lys_C', 'p': 'P', 'specificity': 'K'},
                  '4': {'active': 0, 'name': 'Lys_N', 'p': '-', 'specificity': 'K'},
                  '5': {'active': 1, 'name': 'Arg_C', 'p': 'P', 'specificity': 'R'},
                  '6': {'active': 0, 'name': 'Asp_N', 'p': '-', 'specificity': 'D'},
                  '7': {'active': 1, 'name': 'CNBr', 'p': '-', 'specificity': 'M'},
                  '8': {'active': 1, 'name': 'Glu_C', 'p': 'P', 'specificity': 'DE'},
                  '9': {'active': 1, 'name': 'PepsinA', 'p': 'P', 'specificity': 'FL'},
                  '10': {'active': 1, 'name': 'Chymotrypsin', 'p': 'P', 'specificity': 'FWYL'},
                  '11': {'active': 1, 'name': 'Glu_C/Tryp', 'p': 'P', 'specificity': 'DEKR'}}

defaultVarmods = [{'mass' : 15.994915,
                   'residues' : 'M',
                   'binary' : 0,
                   'max_mods_per_peptide' : 5,
                   'term_distance' : 0,
                   'N/C-term' : 0,
                   'required' : 0}]

def readVarmods():
    from multiplierz import load_mods
    modList = load_mods()
    
    varmods = {}
    for mod, site, mass in modList:
        varmods["%s(%s)" % (mod, site)] = str(mass)
    return varmods
        
            
#unimodDB = unimod.UnimodDatabase(os.path.join(myData, 'unimod.sqlite'))
from multiplierz.mass_biochem import unimod as unimodDB

modLookup = unimodDB.get_pycomet_lookup()
modLookup.update(readVarmods())
modLookup.update({'iTRAQ8plex':{'add_Nterm_peptide':'304.20536', 'add_K_lysine':'304.20536'},
                  'iTRAQ4plex':{'add_Nterm_peptide':'144.10206245', 'add_K_lysine':'144.10206245'},
                  'Carbamidomethyl(C)': {'add_C_cysteine':'57.021464'},
                  'Methylthio(C)': {'add_C_cysteine':'45.98772'}})
                            
toStandardPSMConversions = {'scan':'Query',
                            'num':'Peptide Rank',
                            'charge':'Charge',
                            'exp_neutral_mass':'Predicted mz',
                            'calc_neutral_mass':'Predicted mr',
                            'e-value':'Expectation Value',
                            'xcorr':'Cross-Correlation',
                            'delta_cn':'Delta CN',
                            'sp_score':'SP Score',
                            'ions_matched':'Ions Matched',
                            'ions_total':'Total Predicted Ions',
                            'plain_peptide':'Peptide Sequence',
                            'peptide_modifications':'Variable Modifications', # Old version.
                            'modified_peptide':'Variable Modifications', # New version.
                            'prev_aa':'Preceeding Residue',
                            'next_aa':'Following Residue',
                            'protein':'Accession Number',
                            'duplicate_protein_count':'Protein Count'}

#def convertCometOutputToMascot(cometcsv, mgffile = None):


    #from multiplierz.mzReport import reader, writer
    #psms = reader(cometcsv)
    #output = writer('.'.join(cometcsv.split('.')[:-1] + ['xlsx']),
                    #columns = ['Spectrum Description'] + 
                    #[toStandardPSMConversions[x] for x in psms.columns])
    
    #for psm in psms:
        #convpsm = {}
        #for key in psm:
            #convpsm[toStandardPSMConversions[key]] = psm[key]
        #output.write(convpsm)
    
    #psms.close()
    #output.close()
    
    
    
    
def perform_comet_search(ms2file, database, fixed_mods = None, var_mods = None):
    assert os.path.exists(ms2file), "Input %s not found!" % ms2file
    assert os.path.exists(database), "Database %s not found!" % database
    
    try:
        performSearch, saveParameters = kwargs['runType']
    except KeyError:
        performSearch, saveParameters = True, False 
    if not (performSearch or saveParameters):
        print "Neither 'Perform Search' nor 'Save Parameter File' selected.  Aborting."
        return
    
    searchObj = CometSearch()
    searchObj.database_name = database
    searchObj.allowed_missed_cleavage = missed_cleavages
    searchObj.peptide_mass_tolerance = prectol
    searchObj.peptide_mass_units = unitsDict[precunit]
    searchObj.fragment_bin_tol = fragbintol
    searchObj.fragment_bin_offset = fragbinoffset
    searchObj.num_threads = kwargs['NumThreads']
    searchObj.spectrum_batch_size = kwargs['NumSpectra']
        
    if fixed_mods[:4] == '-FM=':
        fixed_mods = fixed_mods[4:]
    if var_mods[:4] == '-VM=':
        var_mods = var_mods[4:] 
    
    fixed_mod_list = [x.strip() for x in fixed_mods.split(',') if x]
    var_mod_list = [x.strip() for x in var_mods.split(',') if x]    
    
    
    
    
class MS2Writer(object):
    """
    Writer class for .MS2 files, primarily for use in mzFile- or MGF-to-MS2 scripts.
    """
    
    def __init__(self, outputfile, date = None, extractor = None, 
                 extractorVersion = None, extractorOptions = None):
        self.file = open(outputfile, 'w')
        
        self.file.write('H\tCreationDate\t%s\n' % date)
        self.file.write('H\tExtractor\t%s\n' % extractor)
        self.file.write('H\tExtractorVersion\t%s\n' % extractorVersion)
        self.file.write('H\tExtractorOptions\t%s\n' % extractorOptions)
        
    def write(self, spectrum, scanNum, charge, precursorMZ = None, precursorMass = None):
        assert precursorMZ or precursorMass
        if not precursorMZ:
            precursorMZ = ((precursorMass + (charge * protonMass))/charge)
        if not precursorMass:
            precursorMass = (precursorMZ * charge) - (charge * protonMass)
        
        if isinstance(scanNum, tuple):
            self.file.write('S\t%s\t%s\t%s\n' % (scanNum[0], scanNum[1], precursorMZ))
        else:
            try:
                scanNum = int(scanNum)
            except ValueError:
                scanNum = abs(int(hash(scanNum))) % 100000 # Perhaps that'll work?
            self.file.write('S\t%s\t%s\t%s\n' % (scanNum, scanNum, precursorMZ))
        
        self.file.write('Z\t%s\t%s\n' % (charge, precursorMass))
        for pt in spectrum:
            self.file.write('%s\t%s\n' % (pt[0], pt[1]))
        
    def close(self):
        self.file.close()
        
    
    
def mgf_to_ms2(mgfFile, outputfile = None):
    """
    Converts the commonly-used MGF format to an MS2 file.
    """
    
    if not outputfile:
        outputfile = mgfFile[:-4]+'.ms2'
    
    mgf = parse_mgf(mgfFile, header = False)
    ms2 = MS2Writer(outputfile, date = time.localtime(),
                    extractor = 'Converted from MGF by Multiplierz',
                    extractorVersion = __version__,
                    extractorOptions = "None")
    
    sampleTitle = mgf.values()[0]['title']
    if 'MultiplierzMGF' in sampleTitle:
        titleMode = 'standard'
    elif sampleTitle.count('.') > 2 and sampleTitle.split('.')[1] == sampleTitle.split('.')[2] and '|' in sampleTitle:
        titleMode = 'classic' # From that old, barbaric format.
    else:
        titleMode = 'mysterious!'
        print "WARNING: Spectrum title format not recognized; original spectrum numbering will be lost!"
        scanNum = 0
        
    for entry in mgf.values():
        if titleMode == 'standard':
            info = standard_title_parse(entry['title'])
            scanNum = info['scan']
        elif titleMode == 'classic':
            scanNum = entry['title'].split('.')[1]
        else:
            scanNum += 1
            
        ms2.write(entry['spectrum'], 
                  scanNum = scanNum,
                  charge = entry.get('charge', 0), # TODO Add default charge to PyComet?
                  precursorMZ = entry['pepmass'])
    
    ms2.close()
    
    return outputfile
    
    
        
        
findmods = re.compile('\\[(\\d+.\\d+)\\]')
def convertVarmods(psm):
    peptide = psm['Peptide Sequence']
    cometvm = psm['Variable Modifications'][2:-2]
    mods = [(x.start(), x.group(1)) for x in findmods.finditer(cometvm)]
    assert all([cometvm[i-1].isalpha() for i, _ in mods])
    modstrs = ['%s%d: %s' % (cometvm[i-1], i, m) for i, m in mods]
    psm['Variable Modifications'] = '; '.join(modstrs)
    return psm
        
def format_report(reportfile, outputfile = None, mgffile = None, parameters = None,
                  most_rank = None, most_exp = None):
    """
    Renders a native Comet output .txt file into an mzReport-compatible and
    prettier .xlsx format.
    
    (Native .txt output is noncompatible mostly due to a space instead of 
    underscore in the 'modified peptide' column; hopefully that will be fixed
    soon.)
    """
    if most_rank:
        most_rank = int(most_rank)
    if most_exp:
        most_exp = float(most_exp)
    
    if mgffile:
        from multiplierz.mgf import parse_to_generator
        mgfgen = parse_to_generator(mgffile)
        queryToDesc = dict(enumerate(x['title'] for x in mgfgen), start = 1)
    else:
        queryToDesc = {}    
    
    columns = []
    rows = []
    report = open(reportfile, 'r')
    
    headeritems = report.next().split('\t')
    header = {'Program':headeritems[0],
              'Data':headeritems[1],
              'Search Run Time':headeritems[2],
              'Database':headeritems[3].strip()}
    columnline = report.next()
    
    # Fix for presumed bug; omit if this column title is changed in later Comet versions.
    columnline = columnline.replace('peptide\tmodifications', 'peptide_modifications')
    
    def tryNum(thing):
        try:
            return int(thing)
        except ValueError:
            try:
                return float(thing)
            except ValueError:
                return thing
    
    columns = [toStandardPSMConversions.get(x, x) for x in columnline.strip().split('\t')]
    for line in report:
        values = [tryNum(x.strip()) for x in line.split('\t')]
        row = dict(zip(columns, values))
        row = convertVarmods(row)
        row['Spectrum Description'] = queryToDesc.get(row['Query'], 'Unknown')
        rows.append(row)
    
    report.close()
    if not outputfile:
        outputfile = '.'.join(reportfile.split('.')[:-1] + ['xlsx'])
    
    if outputfile.lower().endswith('xlsx') or outputfile.lower().endswith('xls'):
        headerwriter = writer(outputfile, columns = ['Program', 'Data',
                                                     'Search Run Time', 'Database'],
                                          sheet_name = 'Comet_Header')
        headerwriter.write(header)
        headerwriter.write(['', '', '', ''])
        if parameters:
            for setting, value in sorted(parameters.items()):
                headerwriter.write({'Program':setting, 'Data':value,
                                    'Search Run Time':'', 'Database':''})
        headerwriter.close()
    mainwriter = writer(outputfile, columns = ['Spectrum Description'] + columns, sheet_name = 'Data')
    for row in rows:
        if most_rank and row['Peptide Rank'] > most_rank:
            continue
        if most_exp and row['Expectation Value'] > most_exp:
            continue
        mainwriter.write(row)
    mainwriter.close()
    
    return outputfile
    
    
    
        
     
    
    
    
    
    
    
    
def bestType(value):
    try:
        try:
            return int(value)
        except ValueError:
            return float(value)
    except ValueError:
        return str(value)
    
class CometSearch(dict):
    """
    Represents a comet parameters file, so that they can be manipulated via
    Python commands easily, and executes a Comet search using said parameters.
    
    If an existant parameters file is given, this is read in; if a non-existant
    parameters file name is given, defaults are used and saved to the given
    file name; if no file name is given, parameters are saved in a temp file
    only for use with Comet.
    
    .fields gives the available basic parameters; .enzymes contains the list
    of currently given enzymes; .varmods contains the current list of variable
    modifications.  Any of these three lists can be modified, and changes will
    be reflected in the resulting file.
    
    Any field can be accessed as a property of the object; for instance
    
    >settings = CometSearch()
    >settings.database_name
    'some/path/db.fasta'
    >settings.database_name = 'C:/real/path/to/fasta/database.fasta'
    """
    
    def __init__(self, file_name = None, save_parameter_file = False):
        if not (cometPath and os.path.exists(cometPath)):
            raise RuntimeError, ("Comet executable not found at %s; update the "
                                 "multiplierz settings file to indicate your Comet "
                                 "installation." % cometPath)
        
        self.file_name = file_name
        self.fields = []
        self.enzymes = {}
        self.varmods = []
        
        # Parameters for the run itself, used by run_comet_search().
        self.enzyme_selection = None
        
        if file_name:
            with open(file_name, 'r') as parameters:
                for line in parameters:
                    line = line.split('#')[0].strip()
                    
                    if not line:
                        continue
                    
                    try:
                        field, value = line.split('=')
                    except ValueError:
                        if line.strip() == '[COMET_ENZYME_INFO]':
                        #self.enzymes.append(line)                  
                            while True:
                                try:
                                    line = parameters.next()
                                    num, name, active, specificity, p = line.split()
                                except (IOError, ValueError, StopIteration):
                                    break
                                enzyme = {'name' : name,
                                          'active' : int(active),
                                          'specificity' : specificity,
                                          'p' : p
                                          }
                                self.enzymes[num.strip('.')] = enzyme
                            continue # Should be done; could be a return.
                        
                    if field[:12] == 'variable_mod' :
                        words = value.split()
                        if len(words) == 4:
                            mass, res, binary, maxMods = words
                            termDistance = -1
                            ncTerm = 0
                            req = 0
                        else:
                            mass, res, binary, maxMods, termDistance, ncTerm, req = words
                            
                        if float(mass) == 0.0:
                            continue
                        varmod = {'mass' : bestType(mass),
                                  'residues' : res,
                                  'binary' : bestType(binary),
                                  'max_mods_per_peptide' : bestType(maxMods),
                                  'term_distance' : bestType(termDistance),
                                  'N/C-term' : ncTerm,
                                  'required' : req}
                        self.varmods.append(varmod)
                    else:
                        self.fields.append(field.strip())
                        self[field.strip()] = bestType(value.strip())
                        

                        
        else:
            #self.update(defaultParameters)
            for field, value in defaultParameters.items():
                self[field] = value
                self.fields.append(field)
            self.enzymes = defaultEnzymes
            self.varmods = defaultVarmods
        
    def write(self, file_name = None):
        assert file_name or self.file_name, "No output file specified!"
        
        assert len(self.varmods) <= 9, "Comet only allows up to 9 variable modifications."
        
        if (not file_name):
            file_name = self.file_name
        
        
        main_pars = []
        for field in self.fields:
            value = self[field]
            main_pars.append(('%s = %s\n' % (field, value)))
        main_pars.sort()
        with open(file_name, 'w') as parfile:
            for line in parameterPreamble:
                parfile.write(line)
                
            for line in main_pars:
                parfile.write(line)
            
            for num, mod in enumerate(self.varmods, start = 1):
                residues = mod['residues'].upper().replace('N-TERM', 'n').replace('C-TERM', 'c')
                residues = ''.join(x for x in residues if x.isalpha())
                parfile.write('variable_mod0%s = %.4f %s %s %s %s %s %s\n' % (num, mod['mass'], residues,
                                                                              mod['binary'], mod['max_mods_per_peptide'],
                                                                              mod['term_distance'], mod['N/C-term'],
                                                                              mod['required']))
                
            parfile.write('[COMET_ENZYME_INFO]\n')
            for num, enz in self.enzymes.items():
                parfile.write('%s.\t%s\t%s\t%s\t%s\n' % (num, enz['name'], enz['active'],
                                                         enz['specificity'], enz['p']))
            

        
    def run_search(self, data_file, outputfile = None, most_rank = None, most_exp = None, verbose = False):
        #assert self.enzyme_selection != None, "Must specify enzyme selection (attribute .enzyme_selection)!" 
        
        self['digest_mass_range'] = '450.0 6000.0' # Remove this if this parameter gets added to the GUI!
        self['output_txtfile'] = '1'
        
        ext = data_file.split('.')[-1]
        if ext.lower() in ['raw', 'wiff', 'd', 'mzml']:
            from multiplierz.mgf import extract
            if verbose:
                print "Extracting to MGF..."
            data_file = extract(data_file)
            if verbose:
                print "Extracted %s" % datafile
        
        #if not data_file.lower().endswith('ms2'):
            #ms2_file = mgf_to_ms2(data_file)
        
        parfile = os.path.join(myData, 'COMET.par.temp')
        
        while os.path.exists(parfile): # Paranoia; avoid collisions.
            parfile += '.temp'
        self.write(parfile)
        
        if not outputfile:
            outputfile = data_file + '.xlsx'
        
        try:
            expectedResultFile = data_file[:-3] + 'txt'
            
            print('Initiating Comet search...')
            result = call([cometPath, '-P' + parfile, data_file])
            print('Comet search completed with return value %s' % result)    
            assert os.path.exists(expectedResultFile), "Comet failed to produce expected result file."
            
            if outputfile.split('.')[-1].lower() in ['xlsx', 'xls', 'mzd']:
                resultfile = format_report(expectedResultFile, outputfile, data_file,
                                           parameters = dict(self),
                                           most_rank = most_rank,
                                           most_exp = most_exp)
            else:
                resultfile = expectedResultFile
            
            if outputfile:
                os.rename(resultfile, outputfile)
            else:
                outputfile = resultfile
                
            return outputfile
        
        finally:
            os.remove(parfile)

