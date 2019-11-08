import os, sys
from multiplierz import myData
from multiplierz.settings import settings
from subprocess import call
import xml.etree.ElementTree as ET
from multiplierz.mzTools.mzIdentMLAPI import mzIdentML
from multiplierz.mzReport import writer
from multiplierz.mzReport.formats.xtandem import format_XML
import xml.dom.minidom as minidom # For XML pretty-printing.

__all__ = ['TandemSearch']

xtandemExe = settings.xtandem
defaultParameterFile = os.path.join(os.path.dirname(xtandemExe), 'default_input.xml')



defaultParameters = [('list path, default parameters', 'default_input.xml'),
                     ('list path, taxonomy information', None),
                     ('spectrum, fragment monoisotopic mass error', '0.4'),
                     ('spectrum, parent monoisotopic mass error plus', '100'),
                     ('spectrum, parent monoisotopic mass error minus', '100'),
                     ('spectrum, parent monoisotopic mass isotope error', 'yes'),
                     ('spectrum, fragment monoisotopic mass error units', 'Daltons'),
                     ('spectrum, parent monoisotopic mass error units', 'ppm'),
                     ('spectrum, fragment mass type', 'monoisotopic'),
                     ('spectrum, dynamic range', '100.0'),
                     ('spectrum, total peaks', '50'),
                     ('spectrum, maximum parent charge', '4'),
                     ('spectrum, use noise suppression', 'yes'),
                     ('spectrum, minimum parent m+h', '500.0'),
                     ('spectrum, minimum fragment mz', '150.0'),
                     ('spectrum, minimum peaks', '15'),
                     ('spectrum, threads', '1'),
                     ('spectrum, sequence batch size', '1000'),
                     ('spectrum, use neutral loss window', 'yes'), #
                     ('spectrum, neutral loss window', '15.0'),    # From "fully filled out" file.  Valid?
                     ('spectrum, neutral loss mass', '466.0'),     #                    
                     ('residue, modification mass', '57.022@C'),
                     ('residue, potential modification mass', None),
                     ('residue, potential modification motif', None),
                     ('protein, taxon', None),
                     ('protein, cleavage site', '[RK]|{P}'),
                     ('protein, modified residue mass file', None),
                     ('protein, cleavage C-terminal mass change', '+17.002735'),
                     ('protein, cleavage N-terminal mass change', '+1.007825'),
                     ('protein, N-terminal residue modification mass', '0.0'),
                     ('protein, C-terminal residue modification mass', '0.0'),
                     ('protein, homolog management', 'no'),
                     ('refine', 'yes'),
                     ('refine, modification mass', None),
                     ('refine, sequence path', None),
                     ('refine, tic percent', '20'),
                     ('refine, spectrum synthesis', 'yes'),
                     ('refine, maximum valid expectation value', '0.1'),
                     ('refine, potential N-terminus modifications', '+42.010565@['),
                     ('refine, potential C-terminus modifications', None),
                     ('refine, unanticipated cleavage', 'yes'),
                     ('refine, potential modification mass', None),
                     ('refine, point mutations', 'no'),
                     ('refine, use potential modifications for full refinement', 'no'),
                     ('refine, point mutations', 'no'),
                     ('refine, potential modification motif', None),
                     ('scoring, minimum ion count', '4'),
                     ('scoring, maximum missed cleavage sites', '1'),
                     ('scoring, x ions', 'no'),
                     ('scoring, y ions', 'yes'),
                     ('scoring, z ions', 'no'),
                     ('scoring, a ions', 'no'),
                     ('scoring, b ions', 'yes'),
                     ('scoring, c ions', 'no'),
                     ('scoring, cyclic permutation', 'no'),
                     ('scoring, include reverse', 'no'),
                     ('scoring, cyclic permutation', 'no'),
                     ('scoring, include reverse', 'no'),
                     ('output, log path', None),
                     ('output, message', 'testing 1 2 3'),
                     ('output, one sequence copy', 'no'),
                     ('output, sequence path', None),
                     ('output, path', 'output.mzid'),
                     ('output, sort results by', 'protein'),
                     ('output, path hashing', 'yes'),
                     ('output, xsl path', 'tandem-style.xsl'),
                     ('output, parameters', 'yes'),
                     ('output, performance', 'yes'),
                     ('output, spectra', 'yes'),
                     ('output, histograms', 'yes'),
                     ('output, proteins', 'yes'),
                     ('output, sequences', 'yes'),
                     ('output, one sequence copy', 'no'),
                     ('output, results', 'valid'),
                     ('output, maximum valid expectation value', '0.1'),
                     ('output, histogram column width', '30')]

parameterFields = [['list path', 'default parameters'],
                   ['list path', 'taxonomy information'],
                   ['spectrum', 'fragment monoisotopic mass error'],
                   ['spectrum', 'parent monoisotopic mass error plus'],
                   ['spectrum', 'parent monoisotopic mass error minus'],
                   ['spectrum', 'parent monoisotopic mass isotope error'],
                   ['spectrum', 'fragment monoisotopic mass error units'],
                   ['spectrum', 'parent monoisotopic mass error units'],
                   ['spectrum', 'fragment mass type'],
                   ['spectrum', 'dynamic range'],
                   ['spectrum', 'total peaks'],
                   ['spectrum', 'maximum parent charge'],
                   ['spectrum', 'use noise suppression'],
                   ['spectrum', 'minimum parent m+h'],
                   ['spectrum', 'minimum fragment mz'],
                   ['spectrum', 'minimum peaks'],
                   ['spectrum', 'threads'],
                   ['spectrum', 'sequence batch size'],
                   ['spectrum', 'use neutral loss window'],
                   ['spectrum', 'neutral loss window'],
                   ['spectrum', 'neutral loss mass'],
                   ['residue', 'modification mass'],
                   ['residue', 'potential modification mass'],
                   ['residue', 'potential modification motif'],
                   ['protein', 'taxon'],
                   ['protein', 'cleavage site'],
                   ['protein', 'modified residue mass file'],
                   ['protein', 'cleavage C-terminal mass change'],
                   ['protein', 'cleavage N-terminal mass change'],
                   ['protein', 'N-terminal residue modification mass'],
                   ['protein', 'C-terminal residue modification mass'],
                   ['protein', 'homolog management'],
                   ['refine', ''],
                   ['refine', 'modification mass'],
                   ['refine', 'sequence path'],
                   ['refine', 'tic percent'],
                   ['refine', 'spectrum synthesis'],
                   ['refine', 'maximum valid expectation value'],
                   ['refine', 'potential N-terminus modifications'],
                   ['refine', 'potential C-terminus modifications'],
                   ['refine', 'unanticipated cleavage'],
                   ['refine', 'potential modification mass'],
                   ['refine', 'point mutations'],
                   ['refine', 'use potential modifications for full refinement'],
                   ['refine', 'point mutations'],
                   ['refine', 'potential modification motif'],
                   ['scoring', 'minimum ion count'],
                   ['scoring', 'maximum missed cleavage sites'],
                   ['scoring', 'x ions'],
                   ['scoring', 'y ions'],
                   ['scoring', 'z ions'],
                   ['scoring', 'a ions'],
                   ['scoring', 'b ions'],
                   ['scoring', 'c ions'],
                   ['scoring', 'cyclic permutation'],
                   ['scoring', 'include reverse'],
                   ['scoring', 'cyclic permutation'],
                   ['scoring', 'include reverse'],
                   ['output', 'log path'],
                   ['output', 'message'],
                   ['output', 'one sequence copy'],
                   ['output', 'sequence path'],
                   ['output', 'path'],
                   ['output', 'sort results by'],
                   ['output', 'path hashing'],
                   ['output', 'xsl path'],
                   ['output', 'parameters'],
                   ['output', 'performance'],
                   ['output', 'spectra'],
                   ['output', 'histograms'],
                   ['output', 'proteins'],
                   ['output', 'sequences'],
                   ['output', 'one sequence copy'],
                   ['output', 'results'],
                   ['output', 'maximum valid expectation value'],
                   ['output', 'histogram column width']]

requiredParameters = [] # Figure out what these are!

# IT is unclear to me whether XTandem actually uses the non-taxon parts of the
# XML data files.
defaultEnzymes = [('[R]|[A-Z]', 'Arg-C'),
                  ('[A-Z]|[D]', 'Asp-N'),
                  ('[KAY]|[A-Z]', 'Bromelain'),
                  ('[M]|[A-Z]', 'CNBr_HSer'),
                  ('[M]|[A-Z]', 'CNBr_HSerLac'),
                  ('[R]|[A-Z]', 'Cathepsin B'),
                  ('[LF]|{VAG}', 'Cathepsin D'),
                  ('[YWF]|[A-Z]', 'Cathepsin G'),
                  ('[YWFL]|[A-Z]', 'Chymotrypsin'),
                  ('[R]|[P]', 'Clostripain'),
                  ('[AVLIGS]|[A-Z]', 'Elastase'),
                  ('[E]|[A-Z]', 'Glu-C_Bic'),
                  ('[ED]|[A-Z]', 'Glu-C_Phos'),
                  ('[N]|[G]', 'Hydroxylamine'),
                  ('[K]|[A-Z]', 'Lys-C'),
                  ('[A-Z]|[K]', 'Lys-N'),
                  ('[RK]|[A-Z]', 'Papain'),
                  ('[LF]|{VAG}', 'Pepsin'),
                  ('[YWF]|[A-Z]', 'Proteinase K'),
                  ('{RHK}|[A-Z]', 'Subtilisin'),
                  ('[LFIVMA]|{P}', 'Thermolysin'),
                  ('[KR]|{P}', 'Trypsin')]



def writeTaxonomyForFasta(fastafiles, taxon, taxfile = None):
    """
    Writes a taxonomy file with default enzymes and parsers, set to use the
    given FASTA file (overriding notions of 'taxon'.)
    """
    if not taxfile:
        taxfile = fastafiles[0] + '.taxon.xml'
    
    bioml = ET.Element('bioml')
    bioml.set('label', 'x! taxon-to-file matching list')
    
    taxTree = ET.ElementTree(bioml)    
    
    arb = ET.Element('taxon')
    arb.set('label', taxon)
    
    for fastafile in fastafiles:
        fil = ET.Element('file')
        fil.set('format', 'peptide')
        fil.set('URL', fastafile)
        arb.append(fil)
        
    bioml.append(arb)
    taxTree.write(taxfile, xml_declaration = True)
    
    return taxfile


def adaptXtandemXML(xmlfile, inpath, outpath, fasta = None):
    xml = ET.parse(xmlfile)
    bioml = xml.getroot()
    
    for thing in bioml:
        label = thing.get('label')
        if label == 'input, path':
            thing.text = inpath
        elif label == 'output, path':
            thing.text = outpath
        elif label == 'protein, taxon' and fasta:
            thing.text = 'arbitrary'
        elif label == 'list path, taxonomy information' and fasta:
            taxfile = writeTaxonomyForFasta(fasta, 'arbitrary')
            thing.text = taxfile
    
    outputxml = xmlfile[:-4] + '.modified.xml'
    xml.write(outputxml, xml_declaration = True)
    
    return outputxml
    


class TandemSearch(dict):
    """
    Represents an XTandem parameters file, with associated (dummy) taxonomy file.
    
    If an existant parameters file is given, this is read in; if a non-existant
    parameters file name is given, defaults are used and saved to the given
    file name; if no file name is given, parameters are saved in a temp file only
    for use by the script itself.
    
    Fields are accessed through two-layer lookup, where the first layer represents
    the settings categories and the second the setting itself.  For instance:
    
    >settings = TandemParameters()
    >settings['spectrum']['dynamic range'] = 200
    >settings['residue']['modification mass']
    57.022@C
    """
    def __init__(self, parameterfile = None, save_parameter_file = False):
        self.file_name = parameterfile
        self.fields = []
        
        self.save_parameter_file = save_parameter_file
        self.fasta_files = None
        
        if parameterfile and os.path.exists(parameterfile):
            bioml = ET.parse(parameterfile).getroot()
            for thing in bioml:
                fullfield = thing.get('label')
                if not fullfield:
                    continue
                
                try:
                    category, field = [x.strip() for x in fullfield.split(',')]
                    if category not in self:
                        self[category] = {}
                    self[category][field] = thing.text
                except ValueError:
                    field = fullfield.strip()
                    if field not in self:
                        self[field] = {}
                    self[field][''] = thing.text
                        
                
        else: # No file given; use defaults.
            for fullfield, value in defaultParameters:
                try:
                    category, field = [x.strip() for x in fullfield.split(',')]
                    
                    if category not in self:
                        self[category] = {}
                    self[category][field] = value
                except ValueError:
                    assert fullfield == 'refine'
                    if fullfield not in self:
                        self[fullfield] = {}
                    self[fullfield][''] = value
    
        
    def __str__(self):
        output = []
        unset = []
        for category, field in parameterFields:
            if category in self and field in self[category]:
                output.append("\t%s, %s:\t%s\n" % (category, field, self[category][field]))
            else:
                unset.append((category, field))
        
        output.append('\n---------\n')
        output.append('Optional Parameters Not Set:\n')
        output.append('\n'.join([', '.join(x) for x in unset if x not in requiredParameters]))
        output.append('\n---------\n')
        output.append('REQUIRED Parameters Not Set:\n')
        output.append('\n'.join([', '.join(x) for x in unset if x in requiredParameters]))
        output.append('\n---------\n')
        
        return ''.join(output)
    
    def write(self, outputfile = None):
        if not outputfile:
            outputfile = os.path.join(myData, 'XTANDEMTEMP.TEMP')
            self.save_parameter_file = False
        
        bioml = ET.Element('bioml')
        
        note = ET.Element('note')
        note.text = 'XTandem parameter file written by Multiplierz.'
        bioml.append(note)
        
        for category, fields in list(self.items()):
            for field, value in list(fields.items()):
                if category == 'refine' and not field:
                    fullfield = 'refine'
                else:
                    fullfield = ', '.join([category, field])
                
                par = ET.Element('note')
                par.set('type', 'input')
                par.set('label', fullfield)
                par.text = str(value) if value != None else ''
                
                bioml.append(par)
        
        xmlform = ET.ElementTree(bioml)
        domform = minidom.parseString(ET.tostring(xmlform.getroot()))
        output = open(outputfile, 'w')
        output.write(domform.toprettyxml())
        output.close()
        #xml.write(outputfile)
        
        
        return outputfile
    
    
    
    def run_search(self, data_file, outputfile = None, fasta_files = None):
        """
        Runs an XTandem search using the XTandem instance specified in the
        multiplierz settings file.
        
        If fasta_file is specified, this is made to override any existing
        database/taxon listed in the original settings file (if any.)
        """
        assert self.fasta_files or fasta_files, "No sequence database specified!"
        if not fasta_files:
            fasta_files = self.fasta_files
        
        #for i in range(len(fasta_files)):
            #if fasta_files[i].lower().endswith('fasta'):
                #if not os.path.exists(fasta_files[i] + '.pro'):
                    #print "Converting %s..." % fasta_files[i]
                    #fastaConverter = os.path.join(os.path.dirname(xtandemExe), 'fasta_pro.exe')
                    #call([fastaConverter, fasta_files[i]])
                #fasta_files[i] = fasta_files[i] + '.pro'
        
        # This will cause calling run_search to visibly mutate the object!
        # Which may be the only reasonable way to do it- the mutated fields
        # were ignored to begin with, in order to account for the FASTA file.
        
        if not self['protein']['taxon']:
            taxon = 'arbitrary'
            self['protein']['taxon'] = taxon
        else:
            taxon = self['protein']['taxon']
        taxfile = writeTaxonomyForFasta(fasta_files, taxon)
        self['list path']['taxonomy information'] = taxfile
            
        if not self['output']['path']:
            self['output']['path'] = data_file + '.mzid'        
        
        ext = data_file.split('.')[-1]
        if ext.lower() in ['raw', 'wiff', 'd']: # mzml not included, since that can go directly!
            from multiplierz.mgf import extract
            mgf_file = extract(data_file)
        else:
            mgf_file = data_file
        
        if not outputfile:
            outputfile = data_file + '.xlsx'
        else:
            assert outputfile.split('.')[-1] in ['xls', 'xlsx', 'csv', 'mzd']
            
    
        self['spectrum']['path'] = mgf_file
        
        self['list path']['default parameters'] = defaultParameterFile
        
        self['output']['path hashing'] = 'no' # Required to be able to retrieve output.
        expectedOutput = self['output']['path']
        
        parametersFile = self.write()
        
        result = call([xtandemExe, parametersFile])
        assert not result, "XTandem failed with error code %s" % result
        
        assert os.path.exists(expectedOutput), "Output file not found!"
        
        format_XML(expectedOutput, outputfile, parameters = dict(self))
        
        return outputfile
        
    
