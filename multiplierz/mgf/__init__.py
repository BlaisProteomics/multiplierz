from multiplierz.mzAPI import mzFile
from multiplierz.internalAlgorithms import splitOnFirst
from multiplierz import protonMass, __version__, vprint
import os
from numpy import average, std

from multiplierz.mass_biochem import add_protons

__all__ = ['parse_mgf', 'write_mgf', 'MGF_Writer', 'extract', 'apply_spectral_process']
           #'charge_reduce', 'top_n_peaks', 'exclusion_radius', 'signal_noise', 'intensity_threshold']



MGFTopMatter = (['ACCESSION', 'CHARGE', 'CLE', 'COM', 'CUTOUT', 'DB', 'DECOY',
                 'ERRORTOLERANT', 'FORMAT', 'FRAMES', 'INSTRUMENT', 'IT_MODS',
                 'ITOL', 'ITOLU', 'MASS', 'MODS', 'MULTI_SIDE_MODS', 'PEP_SITE_MODS',
                 'PEP_ISOTOPE_ERROR', 'PFA', 'PRECURSOR', 'QUANTITATION', 'REPORT',
                 'REPTYPE', 'SEARCH', 'SEG', 'TAXONOMY', 'TOL', 'TOLU', 'USEREMAIL',
                 'USERNAME'] + ['USER0%s' % x for x in range(0, 9)] +
                ['USER%s' % x for x in range(10, 13)])

# Careful changing this; other multiplierz modules currently
# detect standard title format by checking for this string.
extractor_name = 'MultiplierzMGF'



def standard_title_write(filename, mz = None, scan = None, **info):
    """
    Stores relevant spectral information in a string that can be used as a conveniently
    informative spectrum title, which facilitates downstream processing.
    """

    if isinstance(filename, dict):
        # Can just give an info dict, if e.g. the script has a dict from
        # standard_title_parse and wants to write back modified version.
        info = filename.copy()
        filename = info['file']
        mz = info['mz']
        scan = info['scan']
        del info['file']
        del info['mz']
        del info['scan']
    else:
        assert filename and mz and scan, ("Arguments must either be a dict or "
                                          "include filename, scan, and mz values.")
    
    #title = os.path.basename(filename) + '|%s' % extractorName
    title = '%s|%s|SCAN:%s|MZ:%s' % (os.path.basename(filename), extractor_name,
                                     scan, mz)
    
    for field, value in sorted(info.items()):
        if any(x.isupper() for x in field):
            assert field.lower() not in info, ("Title data should not be case "
                                               "sensitive. (%s <-> %s)") % (field, field.lower())
   
        try: # For field.
            assert '|' not in field and ':' not in field, "| and : cannot appear in field strings.  E.g. %s" % field
        except TypeError:
            pass   
        if isinstance(value, float):
            value = '%.3f' % value if value else '0'
        else:
            try: # For value.
                assert '|' not in value and ':' not in value, "| and : cannot appear in value strings.  E.g. %s" % value
            except TypeError:
                pass
        
        
        title += '|%s:%s' % (field, value)
    
    return title

def standard_title_parse(title):
    """
    Retrieves information from a multiplierz-standard MGF spectrum title.
    """
    
    info = title.split('|')
    filename = info[0]
    data = {'file':filename}
    for fieldval in info[2:]:
        field, value = fieldval.split(':')
        data[field.lower()] = value
    
    return data

def parse_mgf(mgffile, labelType = (lambda x: x), header = False, rawStrings = False):
    """
    Loads a Mascot Generic Format file and returns it in dict form.

    labelType can be a callable object that transforms the 'TITLE='
    value of an MGF entry into what will be used for the corresponding
    key in the dict.  If this is unspecified, the key will be the
    whole TITLE value.
    
    If "header" is set to True the output dict has an extra entry
    of key 'header' that contains the MGF header info.
    
    If raw_strings is set to True, the charge, pepmass, etc are returned
    as 
    """


    f = open(mgffile, "r")

    topMatter = header
    if header:
        data = {'header':{}}
    else:
        data = {}
    for line in f:
        if "BEGIN IONS" in line:
            topMatter = False
            entry = {}
            key = None
            rt = charge = mass = None
            spectrum = []
            for line in f:
                if 'END IONS' in line:
                    break
                elif '=' in line:
                    field, value = splitOnFirst(line, '=')
                    
                    if (not rawStrings) and field == 'CHARGE':
                        value = int(value.strip('\n\r+ .0'))
                    elif (not rawStrings) and field == "PEPMASS":
                        value = float(value.split()[0].strip())
                    else:
                        value = value.strip()
                        
                    entry[field.strip().lower()] = value
                    if field == 'TITLE':
                        key = value.strip()    
                elif any(line):
                    # Ignores fragment charge values.
                    spectrum.append(tuple(map(float, line.split()[:2]))) 

            entry["spectrum"] = spectrum

            key = labelType(key)
            data[key] = entry
        elif topMatter and '=' in line:
            field, value = line.split('=')
            if field in MGFTopMatter:
                data['header'][field.lower()] = value
        elif 'SEARCH=' in line or 'MASS=' in line:
            continue
        else:
            vprint("Unexpected line: %s" % line)

    return data    

def parse_to_generator(mgffile, labelType = (lambda x: x), header = False, rawStrings = False):
    """
    Loads a Mascot Generic Format file and returns it in dict form.

    labelType can be a callable object that transforms the 'TITLE='
    value of an MGF entry into what will be used for the corresponding
    key in the dict.  If this is unspecified, the key will be the
    whole TITLE value.
    
    If "header" is set to True the output dict has an extra entry
    of key 'header' that contains the MGF header info.
    
    If raw_strings is set to True, the charge, pepmass, etc are returned
    as 
    """


    f = open(mgffile, "r")
    topMatter = True
    for line in f:
        if "BEGIN IONS" in line:
            topMatter = False
            entry = {}
            key = None
            rt = charge = mass = None
            spectrum = []
            for line in f:
                if 'END IONS' in line:
                    break
                elif '=' in line:
                    field, value = line.split('=')[:2]
                    
                    if (not rawStrings) and field == 'CHARGE':
                        value = int(value.strip('\n\r+ '))
                    elif (not rawStrings) and field == "PEPMASS":
                        value = float(value.split()[0].strip())
                    else:
                        value = value.strip()
                        
                    entry[field.strip().lower()] = value
                    if field == 'TITLE':
                        key = value.strip()    
                elif any(line):
                    spectrum.append(tuple(map(float, line.split()[:2])))

            entry["spectrum"] = spectrum

            key = labelType(key)
            yield entry
        elif topMatter and '=' in line:
            try:
                field, value = line.split('=')
                if field in MGFTopMatter:
                    #data['header'][field.lower()] = value
                    continue
            except ValueError:
                continue
        elif 'SEARCH=' in line or 'MASS=' in line:
            continue
        else:
            vprint("Unexpected line: %s" % line)




def write_mgf(entries, outputName, header = []):
    output = open(outputName, 'w')

    if isinstance(entries, dict):
        entries = entries.values()

    if isinstance(header, dict):
        header = header.items()

    for field, value in header:
        output.write("%s=%s\n" % (field, value))

    for entry in entries:
        output.write('BEGIN IONS\n')
        output.write('TITLE=%s\n' % entry['title'])
        for field in entry:
            if field.lower() == 'title' or field.lower() == 'spectrum': continue
            output.write('%s=%s\n' % (field.upper(), entry[field]))
        for datapt in entry['spectrum']:
            if len(datapt) == 2:
                mz, intensity = datapt
                output.write('%s\t%s\n' % (mz, intensity))
            elif len(datapt) == 3:
                mz, intensity, chg = datapt
                output.write('%s\t%s\t%s\n' % (mz, intensity, chg))
            else:
                raise NotImplementedError
        output.write('END IONS\n')

    output.close()
    return outputName
    
class MGF_Writer(object):
    def __init__(self, outputfile, header = []):
        self.file = open(outputfile, 'w')
        
        if isinstance(header, dict):
            header = header.items()
        if header:
            for field, value in header:
                self.file.write("%s=%s\n" % (field, value))
    
    def __enter__(self, *args):
        return self
    
    def write(self, scan, title, mass, charge = None):
        self.file.write('BEGIN IONS\n')
        self.file.write('TITLE=%s\n' % title)
        if charge: # Don't even bother with 0 charge, either.
            charge = str(int(charge))
            self.file.write('CHARGE=%s\n' % charge)
        if mass: # 0 mass gives errors, and is usually due to null values.
            self.file.write('PEPMASS=%s\n' % mass)
        
        for pt in scan:
            if len(pt) == 2:
                self.file.write("%s\t%s\n" % (pt[0], pt[1]))
            elif len(pt) >= 3:
                self.file.write("%s\t%s\t%s\n" % (pt[0], pt[1], pt[2]))
            else:
                raise ValueError, "Scan datapoints must have both MZ and intensity!"
        
        self.file.write("END IONS\n")
    
    def add(self, entry):
        self.write(entry['spectrum'], entry['title'],
                   entry['pepmass'], entry.get('charge', None))
    
    def close(self):
        self.file.close()
    
    def __exit__(self, *etc):
        self.close()
        

# 6plex from classic extractor.
example_6plex_dict = {'126':{'127':8.9, '128':0.4},
                      '127':{'126':0.6, '128':8.0, '129':0.4},
                      '128':{'127':1.1, '129':6.5, '130':0.2},
                      '129':{'128':1.6, '130':5.6, '131':0.2},
                      '130':{'128':0.1, '129':1.8, '131':4.8, 'other':0.1},
                      '131':{'129':0.1, '130':3.3, 'other':4.3}}   

example_10plex_dict = {'126':{'127C':5},
                       '127N':{'128N':5.8, 'other':0.2},
                       '127C':{'126':0.3, '128C':4.8},
                       '128N':{'127N':0.3, '129N':4.1},
                       '128C':{'127C':0.6, '129C':3.0},
                       '129N':{'128N':0.8, '130N':3.5},
                       '129C':{'128C':1.4, '130C':2.4},
                       '130N':{'128N':0.1, '129N':1.5, '131':2.4, 'other':3.2},
                       '130C':{'129C':1.8, 'other':2.1},
                       '131':{'130N':1.8, 'other':1.7}}




        
def extract(datafile, outputfile = None, default_charge = 2, centroid = True,
            scan_type = None, deisotope_and_reduce_charge = True,
            maximum_precursor_mass = 15999,
            long_ms1 = True, derive_precursor_via = 'All',
            deisotope_and_reduce_MS1_args = {},
            deisotope_and_reduce_MS2_args = {},
            min_mz = 140, precursor_tolerance = 0.005,
            isobaric_labels = None, label_tolerance = 0.01,
            channel_corrections = None):
    """
    Converts a mzAPI-compatible data file to MGF.
    
    Writes only MS2 spectra where these can be determined, otherwise takes
    every spectrum in the file.  Likewise writes the precursor charge
    and mass if these can be determined.
    
    deisotope_and_reduce_charge deisotopes and charge-reduces each MS2
    spectrum, which generally improves results from peptide database search
    algorithms. However, it should be disabled for very low-resolution scans.
    """
    
    for key, val in [('tolerance', 0.01),
                     ('min_peaks', 2),
                     ('enforce_isotopic_ratios', True)]:
        if key not in deisotope_and_reduce_MS1_args:
            deisotope_and_reduce_MS1_args[key] = val

    if not outputfile:
        outputfile = datafile + '.mgf'

    if os.path.exists(outputfile):
        assert outputfile.lower().endswith('mgf'), ("Overwriting a non-MGF file %s with "
                                                    "the MGF extractor is probably a mistake."
                                                    % outputfile)    
    
    data = mzFile(datafile)
    from extraction import _extractor_
    extractor = _extractor_(data, datafile, default_charge, centroid,
                            scan_type, deisotope_and_reduce_charge, derive_precursor_via,
                            maximum_precursor_mass, long_ms1,
                            deisotope_and_reduce_MS1_args, deisotope_and_reduce_MS2_args,
                            min_mz, precursor_tolerance, isobaric_labels, label_tolerance,
                            channel_corrections)
    writer = MGF_Writer(outputfile)
    
    for scan, title, mz, charge in extractor.run():
        writer.write(scan, title, mass = mz, charge = charge)
    writer.close()
    
    
    if extractor.inconsistent_precursors:
        vprint("Precursor inconsistencies: %s/%s" % (extractor.inconsistent_precursors,
                                                     extractor.scans_written))

    return outputfile    

    



def apply_spectral_process(mgfFile, functions, outputFile = None):
    """
    Takes an MGF file name and a list of functions together with their
    arguments, and applies each function; compatible with the 
    functions from spectral_process.py that act on scan spectra.
    """
    
    mgf = parse_mgf(mgfFile)
    
    if not outputFile:
        outputFile = mgfFile
        
    def applyFunctions(spectrum):
        for args, func in functions:
            if args:
                spectrum = func(spectrum, args)
            else:
                spectrum = func(spectrum)
        return spectrum
        
    try:
        header = mgf['header']
        del mgf['header']
    except KeyError:
        header = None
    with MGF_Writer(outputFile, header) as writer:
        for entry in mgf.values():
            entry['spectrum'] = applyFunctions(entry['spectrum'])
            writer.write(entry['spectrum'], entry['title'],
                         mass = entry['pepmass'], charge = entry.get('charge', None))
    
    return outputFile
    


        
if __name__ == '__main__':
    import cProfile as profile
    #extract(r'\\rc-data1\blaise\ms_data_share\Max\CSF\example_RAW\2017-04-11-CSF-BIOFIND-SET5-PROTEOME-TMT-3D-24-200.RAW', outputfile = r'C:\Users\Max\Desktop\ExampleData\foobar.mgf', long_ms1 = True)
    profile.run(r"extract(r'\\rc-data1\blaise\ms_data_share\Max\CSF\example_RAW\2017-04-11-CSF-BIOFIND-SET5-PROTEOME-TMT-3D-24-200.RAW', long_ms1 = True)",
                r'C:\Users\Max\Desktop\Projects\extractor_stats_newmode_3')