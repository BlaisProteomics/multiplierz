import os, sys
import re
import time
import multiplierz.unimod as unimod
from multiplierz import myData
from subprocess import call
from collections import defaultdict
from multiplierz.mgf import parse_mgf, extractor_name, standard_title_parse
from multiplierz.mzAPI import mzFile
from multiplierz import protonMass, __version__
import xml.etree.ElementTree as ET
import csv



__all__ = ['CometSearch']


parameterPreamble = ['# comet_version 2015.01 rev. 2\n',
                     '# Comet MS/MS search engine parameters file.\n',
                     "# Everything following the '#' symbol is treated as a comment.\n"]

# For instance, matches 'comet.2015012.win64.exe'
cometRegex = r'comet\.20[0-9]{5}\.win(64|32)\.exe'
cometExecutable = re.compile(cometRegex)

cometPath = None
if not sys.executable.lower().endswith('python.exe'):
    # mzDesktop mode; use packaged-with copy of Comet.
    sourcePath = os.path.dirname(sys.executable)
    for path, _, files in os.walk(sourcePath):
        executables = [x for x in files if cometExecutable.match(x)]
        if executables:
            cometPath = os.path.join(path, executables[0])
            break
    assert cometPath, "No Comet executable matching pattern %s found in %s" % (cometRegex,
                                                                               sourcePath)

if not cometPath:
    # mzPy mode (also mzDesktop failback.)
    from multiplierz import SettingsFile, myData
    import multiplierz.settings as settings
    cometPath = settings.comet
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
                     'use_sparse_matrix': '1',
                     'variable_mod01': '15.9949 M 0 3 -1 0 0',
                     'variable_mod02': '0.0 X 0 3 -1 0 0',
                     'variable_mod03': '0.0 X 0 3 -1 0 0',
                     'variable_mod04': '0.0 X 0 3 -1 0 0',
                     'variable_mod05': '0.0 X 0 3 -1 0 0',
                     'variable_mod06': '0.0 X 0 3 -1 0 0',
                     'variable_mod07': '0.0 X 0 3 -1 0 0',
                     'variable_mod08': '0.0 X 0 3 -1 0 0',
                     'variable_mod09': '0.0 X 0 3 -1 0 0'}

# This needs '[COMET ENZYME INFO]' placed in front of it in the param file.
defaultEnzymes = ['0.  No_enzyme              0      -           -',
                  '1.  Trypsin                1      KR          P',
                  '2.  Trypsin/P              1      KR          -',
                  '3.  Lys_C                  1      K           P',
                  '4.  Lys_N                  0      K           -',
                  '5.  Arg_C                  1      R           P',
                  '6.  Asp_N                  0      D           -',
                  '7.  CNBr                   1      M           -',
                  '8.  Glu_C                  1      DE          P',
                  '9.  PepsinA                1      FL          P',
                  '10. Chymotrypsin           1      FWYL        P']


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




toStandardPSMConversions = {'assumed_charge':'Charge',
                            'calc_neutral_pep_mass':'Predicted mr',
                            #'deltacn':'Delta',
                            'hit_rank':'Peptide Rank',
                            'index':'Query',
                            'massdiff':'Delta',
                            'modified_peptide':'Variable Modifications',
                            'peptide':'Peptide Sequence',
                            'peptide_next_aa':'Following Residue',
                            'peptide_prev_aa':'Preceding Residue',
                            'precursor_neutral_mass':'Experimental mz', # But its actually not the mz, is it?
                            'protein':'Accession Number',
                            'spectrum':'Spectrum Description',
                            # Thats all the easily Mascot-equivalents.  Here's the rest:
                            'deltacn':'DeltaCN',
                            'deltacnstar':'DeltaCN*',
                            'end_scan':'End Scan',
                            'expect':'Expect',
                            'num_matched_ions':'Matched Ions #',
                            'num_matched_peptides':'Matched Peptides #',
                            'num_missed_cleavages':'Missed Cleavages',
                            'num_tol_term':'num_tol_term', # I DON'T KNOW.
                            'num_tot_proteins':'Total Proteins #',
                            'retention_time_sec':'Retention Time (sec)',
                            'sprank':'SP Rank',
                            'spscore':'SP Score',
                            'start_scan':'Start Scan',
                            'tot_num_ions':'Total Ion #',
                            'xcorr':'Cross-Correlation'
                            }

def convertCometOutputToMascot(cometcsv):
    from multiplierz.mzReport import reader, writer
    psms = reader(cometcsv)
    output = writer('.'.join(cometcsv.split('.')[:-1] + ['xlsx']),
                    columns = [toStandardPSMConversions[x] for x in psms.columns])
    
    for psm in psms:
        convpsm = {}
        for key in psm:
            convpsm[toStandardPSMConversions[key]] = psm[key]
        output.write(convpsm)
    
    psms.close()
    output.close()
    




def renderFixedModSettings(fixmodstrs):
    fixSettings = defaultdict(float)
    fixSettings['add_C_cysteine'] = 0.0
    for fixstr in fixmodstrs:
        value = modLookup[fixstr]
        if not isinstance(value, dict):
            value = float(value)
            specificities = list(fixstr.split('(')[1].split(')')[0])
            for aa in specificities:
                aminoStr = fixModTypes[aa]
                fixSettings[aminoStr] += value
        
        else:
            for aminoStr, mass in value.items():
                fixSettings[aminoStr] += float(mass)
    
    return dict(fixSettings)

def renderVarModSettings(varmodstrs):
    varSettings = {'variable_mod01': '0.0 X 0 3 -1 0 0'}
    modBase = 'variable_mod0'
    
    for index, varstr in enumerate(varmodstrs, start = 1):
        field = modBase + str(index)
        value = modLookup[varstr]
        
        if not isinstance(value, dict):
            value = float(value)
            specificities = varstr.split('(')[1].split(')')[0]
            
            varValue = "%.4f %s 0 3 -1 0 0" % (value, specificities)
            varSettings[field] = varValue
        
        else:
            aas = ''
            for aminoStr, mass in value.items():
                aas += fixedModTypesBackward[aminoStr]
            varValue = "%.4f %s 0 3 -1 0 0" % (float(mass), aas)
            varSettings[field] = varValue
    
    return dict(varSettings)
            
            
            
            
parameterSequence = None

def readParameters(parameterFile):
    global parameterSequence
    parameterSequence = []
    
    pars = {}
    with open(parameterFile) as parameters:
        for line in parameters:
            if '#' in line:
                commentFrom = line.index('#')
                line = line[:commentFrom]
            line = line.strip()
            try:
                field, value = [x.strip() for x in line.split('=')]
                pars[field] = value
                parameterSequence.append(field)
            except ValueError:
                pass
    
    return pars

def writeParameters(parameters, parameterFile, enzyme):
    with open(parameterFile, 'w') as parOut:
        for line in parameterPreamble:
            parOut.write(line)
            
        #for field, value in sorted(parameters.items()):
        for field in parameterSequence:
            value = parameters[field]
            parOut.write('%s = %s\n' % (field, value))
            
        parOut.write('[COMET_ENZYME_INFO]\n')
        enzymeProperties = enzymeDict[enzyme]
        parOut.write('0.\t%s\t%s\t%s\t%s\n' % tuple([enzyme] + enzymeProperties))


def perform_comet_search(ms2file, database, fixed_mods=None, var_mods=None,
                         fdr="-PYes", enzyme = 'Trypsin', missed_cleavages='2',
                         prectol='10', precunit='ppm', fragbintol='1.0005',
                         fragbinoffset='0.4', **kwargs):
    # kwargs may contain: runType, parent, NumThreads, NumSpectra, ExpDim, fdrSort
    
    assert os.path.exists(ms2file), "Input %s not found!" % ms2file
    assert os.path.exists(database), "Database %s not found!" % database
    
    parentGUI = kwargs.get('parent', None)
    if parentGUI:
        parentGUI.post_message("Parsing MS2 spectra. " + parentGUI.make_time(time.localtime()))
        parentGUI.post_message(ms2file)         
    
    try:
        performSearch, saveParameters = kwargs['runType']
    except KeyError:
        performSearch, saveParameters = True, False 
    if not (performSearch or saveParameters):
        print "Neither 'Perform Search' nor 'Save Parameter File' selected.  Aborting."
        return    
    
    baseParameterFile = os.path.join(os.path.dirname(cometPath), 'comet.params.new')
    if not os.path.exists(baseParameterFile):
        os.chdir(os.path.dirname(cometPath))
        call([cometPath, '-p'])
        assert os.path.exists(baseParameterFile)
        os.chdir(myData)
    
    parameters = readParameters(baseParameterFile)
    
    parameters['database_name'] = database
    parameters['search_enzyme_number'] = '0' # COMET_ENZYME_INFO constructed accordingly.
    parameters['allowed_missed_cleavage'] = missed_cleavages
    parameters['peptide_mass_tolerance'] = prectol
    parameters['peptide_mass_units'] = unitsDict[precunit]
    parameters['fragment_bin_tol'] = fragbintol
    parameters['fragment_bin_offset'] = fragbinoffset
    parameters['num_threads'] = kwargs['NumThreads']
    parameters['spectrum_batch_size'] = kwargs['NumSpectra']
    
    exp = kwargs['ExpDim'] == '1D'
    fdrSort = kwargs['fdrSort']
    
    
    # Remove once pyComet.py is fixed.
    if fixed_mods[:4] == '-FM=':
        fixed_mods = fixed_mods[4:]
    if var_mods[:4] == '-VM=':
        var_mods = var_mods[4:]     
        
    fixed_mod_list = [x.strip() for x in fixed_mods.split(',') if x]
    var_mod_list = [x.strip() for x in var_mods.split(',') if x]
    
    fixedModFields = renderFixedModSettings(fixed_mod_list)
    varModFields = renderVarModSettings(var_mod_list)
    
    parameters.update(fixedModFields)
    parameters.update(varModFields)
    
    writtenParameters = ms2file + '.comet.params'
    writeParameters(parameters, writtenParameters, enzyme)
    
    if performSearch:
        print('Initiating Comet search...')
        result = call([cometPath, '-P' + writtenParameters, ms2file])
        print('Comet search completed with return value %s' % result)
        
        expectedResultFile = ms2file[:-3] + 'pep.xml'
        assert os.path.exists(expectedResultFile)
    
        print('Processing XML...')
        #csvResultFile = convertToSpreadsheet(expectedResultFile)
        csvResultFile = expectedResultFile[:-3] + 'csv'
        process_file(expectedResultFile)
        assert os.path.exists(csvResultFile)
        
        if fdr:
            #csvResultFile = cometFDR(csvResultFile)
            fdrCSVFile = csvResultFile[:-4] + '.fdr.csv'
            calc_fdr(csvResultFile, fdr_filter = 0.01,
                                       sort_column=fdrSort,
                                       rev_txt=kwargs['RevDbId'])
            assert os.path.exists(fdrCSVFile)
            try:
                os.remove(csvResultFile)
            except (WindowsError, IOError) as err:
                print err
            os.rename(fdrCSVFile, csvResultFile)
            
        
        os.remove(expectedResultFile)
        print "Result file written to %s" % csvResultFile
    
    if not saveParameters:
        os.remove(writtenParameters)
        
    print "Done!"
    
    
    
    
def perform_comet_search2(ms2file, database, fixed_mods = None, var_mods = None):
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
    if extractorName in sampleTitle:
        titleMode = 'standard'
    elif sampleTitle.split('.')[1] == sampleTitle.split('.')[2] and '|' in sampleTitle:
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
                  charge = entry['charge'],
                  precursorMZ = entry['pepmass'])
    
    ms2.close()
    
    return outputfile
    
    
        
        
        
        
        
        
        
        
        
        
        
def calc_fdr(filename, fdr_filter = None, sort_column='expect', rev_txt='rev_gi'):
    '''
    
    Calculates fdr for COMET Search data
    
    INPUT: COMET output as csv (generated by process_file.py)
    
    If fdr filter is given a value, rows containing reverse hits or that have an fdr > the given value
    are not written to the result file.
    
    OUTPUT: csv file
    
    '''
    sort_reverse=False
    if sort_column in ['xcorr']:
        sort_reverse=True
    rdr = csv.reader(open(filename, 'r'))
    data = []
    headers = None
    counter = 0
    print "Reading file in memory..."
    for row in rdr:
        if counter:
            data.append(row)    
            data[counter-1][headers.index(sort_column)] = float(data[counter-1][headers.index(sort_column)])                   
        else:
            headers = [x.strip() for x in row] + ['fdr']
        counter += 1
        
    print "Sorting..."
    
    data.sort(key=lambda t:t[headers.index(sort_column)], reverse=sort_reverse)
    
    fwd = 0
    rev = 0
    
    print "Calc fdr..."
    for i, member in enumerate(data):
        if member[headers.index('protein')].find(rev_txt) == -1:
            fwd += 1
        else:
            rev += 1
        try:
            fdr = float(rev)/float(fwd)
        except:
            fdr = 'NA'
        data[i].append(fdr)
    
    print "Writing..."
    wtr = csv.writer(open(filename[:-4] + '.fdr.csv', 'wb'), delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    wtr.writerow(headers)
    for row in data:
        if not fdr_filter:
            wtr.writerow(row)
        else:
            if row[headers.index('fdr')] < fdr_filter and row[headers.index('protein')].find(rev_txt) == -1:
                wtr.writerow(row)
    del wtr
    del rdr
    del data
            

def process_file(filename, first_rank_only=True):
    '''
    
    GENERATES A CSV File from a pepxml file
    
    USAGE:
    
    filename = r'C:\Data\COMET\2012-08-23-Enolase-100fmol_1_RECAL.pep.xml'
    
    process_file(filename)
    
    '''    
    tree = ET.parse(filename)
    root = tree.getroot()
    for child in root:
        for s in child.findall('.'): # {http://regis-web.systemsbiology.net/pepXML}spectrum_query
            a = s.getchildren()
    hit_list = []
    master_keys = []
    for member in a:
        if member.tag == '{http://regis-web.systemsbiology.net/pepXML}spectrum_query':
            query_data = {}
            x = member.getchildren()
            stuff = member.attrib
            query_data.update(stuff)
            for query in x:
                    z = query.getchildren() # this is search_hit
                    for y in z: # y = search_hit
                        hit_data = {}
                        for key in y.attrib.keys():
                            hit_data[key]=y.attrib[key] 
                        u = y.getchildren()
                        score_dict = {}
                        for i in u:
                            if 'name' in i.keys():
                                try:
                                    score_dict[i.attrib['name']]=i.attrib['value']
                                except:
                                    pass
                            else:
                                score_dict.update(i.attrib)
                            
                        #CONSTRUCT "HIT"
                        hit_dict = {}
                        hit_dict.update(query_data)
                        hit_dict.update(hit_data)
                        hit_dict.update(score_dict)
                        if 'xcorr' in hit_dict.keys():
                            hit_list.append(hit_dict)
                            if len(hit_dict.keys()) > len(master_keys):
                                master_keys = hit_dict.keys()
    wtr = csv.writer(open(filename[:-4] + '.csv', 'wb'), delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    x = master_keys
    x.sort()
    #print x
    wtr.writerow(x)
    for row in hit_list:
        proceed = False
        if first_rank_only:
            if row['hit_rank'] == '1':
                proceed = True
        else:
            proceed = True
        if proceed:
            elem = []
            x = master_keys
            x.sort()
            for key in x:
                if key in row.keys():
                    elem.append(row[key])
                else:
                    elem.append('')
            wtr.writerow(elem)
    del wtr
    del hit_list
    del master_keys
    del root
    del tree
    
    
    
    
    
    
    
    
    
    
    
def bestType(value):
    try:
        try:
            return int(value)
        except ValueError:
            return float(value)
    except ValueError:
        return str(value)
    
class CometSearch(object):
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
                                 "installation.")
        
        self.file_name = file_name
        self.fields = []
        self.enzymes = []
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
                                self.enzymes.append(enzyme)
                    
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
                        self.__dict__[field.strip()] = bestType(value.strip())
                        

                        
        else:
            #self.update(defaultParameters)
            for field, value in defaultParameters.items():
                self.__dict__[field] = value
                self.fields.append(field)
            self.enzymes = defaultEnzymes
            # self.varmods = defaultVarmods
        
    def write(self, file_name = None):
        assert file_name or self.file_name, "No output file specified!"
        
        assert len(self.varmods) <= 9, "Comet only allows up to 9 variable modifications."
        
        if (not file_name):
            file_name = self.file_name
        
        
        with open(file_name, 'w') as parfile:
            for line in parameterPreamble:
                parfile.write(line)
                
            for field in self.fields:
                value = self.__dict__[field]
                parfile.write('%s = %s\n' % (field, value))
            
            for num, mod in enumerate(self.varmods, start = 1):
                parfile.write('variable_mod0%s = %.4f %s %s %s %s %s %s\n' % (num, mod['mass'], mod['residues'],
                                                                              mod['binary'], mod['max_mods_per_peptide'],
                                                                              mod['term_distance'], mod['N/C-term'],
                                                                              mod['required']))
                
            parfile.write('[COMET_ENZYME_INFO]\n')
            for num, enz in enumerate(self.enzymes):
                parfile.write('%s.\t%s\t%s\t%s\t%s\n' % (num, enz['name'], enz['active'],
                                                         enz['specificity'], enz['p']))
            

        
    def run_search(self, data_file):
        #assert self.enzyme_selection != None, "Must specify enzyme selection (attribute .enzyme_selection)!" 
        
        ext = data_file.split('.')[-1]
        if ext.lower() in ['raw', 'wiff', 'd', 'mzml']:
            from multiplierz.mgf import extract
            data_file = extract(data_file)
        
        if not data_file.lower().endswith('ms2'):
            ms2_file = mgf_to_ms2(data_file)
        
        parfile = os.path.join(myData, 'COMET.par.temp')
        
        while os.path.exists(parfile): # Paranoia; avoid collisions.
            parfile += '.temp'
        self.write(parfile)
        
        try:
            expectedResultFile = ms2_file[:-3] + 'pep.xml'
            
            print('Initiating Comet search...')
            result = call([cometPath, '-P' + parfile, ms2_file])
            print('Comet search completed with return value %s' % result)    
            assert os.path.exists(expectedResultFile), "Comet failed to produce expected result file."
            
            os.remove(parfile)
            
            csvResultFile = expectedResultFile[:-3] + 'csv'
            process_file(expectedResultFile)
            assert os.path.exists(csvResultFile)    
            
            return csvResultFile
        
        except Exception as err:
            os.remove(parfile)
            raise err
        
    
if __name__ == '__main__':
    print 'TEST MODE'
    foo = CometSearch(r'C:\Users\Max\Downloads\comet_binaries_2015012\comet.params.new')
    foo.database_name = r'C:\Users\Max\Desktop\Dev\coverageThings\HUMAN.fasta'
    result = foo.run_search(r'C:\Users\Max\Desktop\SpectrometerData\2014-07-23-Mito-ISD-SDS.CID_ITMS.mgf')
    print "FOO"
    

#if __name__ == '__main__':
    #foo = r'C:\Users\Max\Desktop\SpectrometerData\2014-07-23-Mito-ISD-SDS.CID_ITMS.pep.csv'
    #convertCometOutputToMascot(foo)