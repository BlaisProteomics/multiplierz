import multiplierz.mzAPI.raw
from multiplierz.mzAPI import mzFile
from multiplierz import protonMass, __version__
import os
from numpy import average, std



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
            
        try:
            assert '|' not in field and ':' not in field, "| and : cannot appear in field strings.  E.g. %s" % field
        except TypeError:
            pass
        try:
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
                    field, value = line.split('=')
                    
                    if (not rawStrings) and field == 'CHARGE':
                        value = int(value.strip('\n+ '))
                    elif (not rawStrings) and field == "PEPMASS":
                        value = float(value.split()[0].strip())
                    else:
                        value = value.strip()
                        
                    entry[field.strip().lower()] = value
                    if field == 'TITLE':
                        key = value.strip()    
                elif any(line):
                    spectrum.append(tuple(map(float, line.split())))

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
            print "Unexpected line: %s" % line

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
                    field, value = line.split('=')
                    
                    if (not rawStrings) and field == 'CHARGE':
                        value = int(value.strip('\n+ '))
                    elif (not rawStrings) and field == "PEPMASS":
                        value = float(value.split()[0].strip())
                    else:
                        value = value.strip()
                        
                    entry[field.strip().lower()] = value
                    if field == 'TITLE':
                        key = value.strip()    
                elif any(line):
                    spectrum.append(tuple(map(float, line.split())))

            entry["spectrum"] = spectrum

            key = labelType(key)
            yield entry
        elif topMatter and '=' in line:
            field, value = line.split('=')
            if field in MGFTopMatter:
                #data['header'][field.lower()] = value
                continue
        elif 'SEARCH=' in line or 'MASS=' in line:
            continue
        else:
            print "Unexpected line: %s" % line  




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
    
    def close(self):
        self.file.close()
    
    def __exit__(self, *etc):
        self.close()
        
    

#ScanFeatureData = {} # Isotopic feature data per scan; should be refreshed per file processed!
    
#def _get_precursor(pMZ, scan, scanNum, charge, tolerance):
    #"""
    #Used in extract(); attempts to find the charge of a given precursor peak
    #in an MS1 scan.   
    #"""
    
    #scan = [x[:2] for x in scan]
    
    #global ScanFeatureData
    #if scanNum not in ScanFeatureData:
        #c12Table = {}
        #for chg in range(1, 20):
            #c12Table[chg] = [] # To avoid gaps.
            
        #featureData, etc = peak_pick(scan, max_charge = 10, tolerance = tolerance)
        #for chg, isotopeSeqs in featureData.items():
            #for iso in isotopeSeqs:
                ##c12Table.append((iso, chg))
                #c12Table[chg].append((iso, chg))
        
        #ScanFeatureData[scanNum] = c12Table
    #else:
        #c12Table = ScanFeatureData[scanNum]
    
    #if charge:
        #nearSequences = [iso for iso, chg in c12Table[charge]
                         #if iso[-1][0] > pMZ - 0.05 and iso[0][0] < pMZ + 0.05]
        ##                Matching charge, and given mz falling (given tolerances)
        ##                somewhere within the range of the isotopic envelope.
    #else:
        ## If not given a charge, we want to recover it from the isotopes table.
        #nearSequences, charges = zip(*[(iso, chg) for iso, chg in sum(c12Table.values(), []) if
                                       #iso[-1][0] > pMZ - 0.05 and iso[0][0] < pMZ + 0.05])
        
    #if not nearSequences:
        #precMZ = None
    #else:
        #if len(nearSequences) == 1: # Efficiency; assume its right if its the only option.
            #seq = nearSequences[0]
            #if not charge:
                #charge = charges[0]
                
        #else: # Multiple possible isotopic sequences.
            #seqIndex, seq = min(enumerate(nearSequences),
                                #key = lambda x: min([abs(mz - pMZ) for mz, int in x[1]]))
            ##                   Take the isotopic sequence with the nearest peak.
            
            #if not charge:
                #charge = charges[seqIndex]
                ## ...and the corresponding charge, if necessary.
        
        ##                                 Most intense peak in the isotopic sequence.
        #highIndex, highPt = max(enumerate(seq), key = lambda x: x[1][1])
        
        ## highPt is likely the C12 if its the leftmost peak, otherwise its likely
        ## the C13.  In the latter case, if its not the next-to-leftmost, the
        ## one beside it is considered the C12 and any peaks further left than
        ## that are ???mysteries???
        #if highIndex > 0:
            #precMZ = seq[highIndex - 1][0]
        #else:
            #precMZ = highPt[0]
        ## Complicated!
        
    #return precMZ, charge
        
    
                
                    
                
                
    
    
def extract(datafile, outputfile = None, default_charge = 2, centroid = True,
            scan_type = None, deisotope_and_reduce_charge = True,
            min_mz = 140, precursor_tolerance = 0.005,
            isobaric_labels = None, label_tolerance = 0.01):
    """
    Converts a mzAPI-compatible data file to MGF.
    
    Writes only MS2 spectra where these can be determined, otherwise takes
    every spectrum in the file.  Likewise writes the precursor charge
    and mass if these can be determined.
    
    deisotope_and_reduce_charge deisotopes and charge-reduces each MS2
    spectrum, which generally improves results from peptide database search
    algorithms. However, it should be disabled for very low-resolution scans.
    """
    # Currently doesn't compensate for injection time! Would be required in
    # order to deal with iTRAQ/TMT labels.
    
    from multiplierz.spectral_process import deisotope_reduce_scan, peak_pick, centroid
    
    def _get_precursor(mz, possible_prec, charge):
        try:
            return min([x for x in possible_prec
                        if (charge == None or x[1] == charge)],
                       key = lambda x: abs(x[0] - mz))
        except ValueError:
            return None, None

    if not outputfile:
        outputfile = datafile + '.mgf'
        
    if os.path.exists(outputfile):
        assert outputfile.lower().endswith('mgf'), ("Overwriting a non-MGF file %s with "
                                                    "the MGF extractor is probably a mistake."
                                                    % outputfile)
    
    writer = MGF_Writer(outputfile)    
    
    data = mzFile(datafile)
    scanInfo = data.scan_info()
    
    # Coerce that scanInfo be in order of time, so that for .WIFF files
    # we can still use the previous-MS1 method to look up precursor charges.
    scanInfo.sort(key = lambda x: x[0])
    
    if datafile.lower().endswith('.raw'): # May also exist for WIFF?
        filters = dict(data.filters())
        
        # For RAW files only, there's the option to filter by a given
        # scan type.  (It would be more efficient in many cases to
        # actually split files in a single run, though.)
        if scan_type:    
            scan_type = scan_type
            assert (scan_type.lower() in
                    ['cid', 'hcd', 'etd', 'etdsa']), ("Invalid scan type %s, must be one"
                                                      "of (CID, HCD, ETD, ETDSA).") % scan_type
            typestr = "@%s" % scan_type.lower()
            
            scanInfo = [x for x in scanInfo if x[3] == 'MS1' or
                        typestr in filters[x[0]]]
    else:
        filters = None
        assert not scan_type, "Scan type filtering only enabled with .RAW format files."
    
    
    if isobaric_labels:
        assert centroid, "Isobaric tags can only be read from centroided data; set 'centroid' to True."
    
    if not isobaric_labels:
        labels = []
    elif isobaric_labels == 4 or isobaric_labels == '4plex':
        labels = zip([114, 115, 116, 117], [114.11, 115.11, 116.11, 117.12])
    elif isobaric_labels == 6 or isobaric_labels == '6plex':
        labels = zip([126, 127, 128, 129, 130, 131],
                     [126.127, 127.131, 128.134, 129.138, 130.141, 131.138])
    elif isobaric_labels == 8 or isobaric_labels == '8plex':
        labels = zip([113, 114, 115, 116, 117, 118, 119, 121],
                     [113.11, 114.11, 115.11, 116.11, 117.12, 118.12, 119.12, 121.12])
    elif isobaric_labels == 10 or isobaric_labels == '10plex':
        labels = zip(['126', '127N', '127C', '128N', '128C', 
                      '129N', '129C', '130N', '130C', '131'],
                     [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                      129.131471, 129.137790, 130.134825, 130.141145, 131.138180])
        
        assert label_tolerance < 0.005, ("label_tolerance must be lower "
                                         "than 0.005 for 10-plex experiments!" % label_tolerance)
    else:
        raise NotImplementedError, ("Labels of type %s not recognized.\n"
                                    "Should be one of [4,6,8,10] or None.")
            
    def read_labels(scan):
        partscan = [x for x in scan if x[0] < labels[-1][1] + 3]
        if not partscan:
            return dict([(str(l), '0') for l in zip(*labels)[0]])
        
        scan_values = {}
        for label, mz in labels:
            nearpt = min(partscan, key = lambda x: abs(x[0] - mz))
            if abs(nearpt[0] - mz) < label_tolerance:
                scan_values[str(label)] = '%.3f' % nearpt[1]
            else:
                scan_values[str(label)] = '0' # Report noise value?
                    
        return scan_values
    
    
    inconsistent_precursors = 0
    scans_written = 0
    
    lastMS1 = None
    for time, mz, scanNum, scanLevel, scanMode in scanInfo:
        scanName = scanNum if isinstance(scanNum, int) else time
        
        if scanLevel == 'MS1':
            lastMS1ScanName = scanName
            
            possible_precursors = None
            def calculate_precursors():
                if isinstance(data, multiplierz.mzAPI.raw.mzFile):
                    lastMS1 = data.lscan(lastMS1ScanName)
                else:
                    try:
                        lastMS1 = data.scan(lastMS1ScanName,
                                            centroid = True)
                    except NotImplementedError:
                        lastMS1 = centroid(data.scan(lastMS1ScanName))      
                        
                envelopes = peak_pick(lastMS1, tolerance = 0.01, min_peaks = 2,
                                      enforce_isotopic_ratios=False)[0]
                return sum([[(x[0][0], c) for x in xs]
                            for c, xs in envelopes.items()], [])
            
            continue
        
        try:
            scan = data.scan(scanName, centroid = True)
        except NotImplementedError as err:
            # Currently calls to mzWiff with centroid enabled return
            # an error.
            if centroid:
                scan = centroid(data.scan(scanName))
            else:
                raise err
        
        if filters and not mz:
            mz = float(filters[time].split('@')[0].split(' ')[-1])            
        
        mzP = None
        chargeP = None
        if "scanPrecursor" in dir(data):
            assert isinstance(scanName, int)
            mzP, chargeP = data.scanPrecursor(scanName)
        
        if not mzP: # .scanPrecursor sometimes returns charge and not mzP.
            if possible_precursors == None:
                possible_precursors = calculate_precursors()
                
            mzP, chargeP = _get_precursor(mz, possible_precursors, chargeP)
            if not mzP:
                # Release presumed charge possibly obtained from scanPrecursor.
                mzP, chargeP = _get_precursor(mz, possible_precursors, None)
                if mz and chargeP:
                    inconsistent_precursors += 1
        
            

            
        if mzP and (abs(mz - mzP) < 2 or not mz): 
            mz = mzP
            charge = chargeP
        else:
            charge = default_charge
        
        if not charge:
            charge = default_charge
        


        
        if not mz:
            import warnings
            errmgf = os.path.abspath(datafile)
            warnings.warn('Unable to recover all precursor masses from %s' % errmgf)
        else:    
            if labels:
                scan_labels = read_labels(scan)
            else:
                scan_labels = {}
            
            title = standard_title_write(datafile, rt = time, mz = mz,
                                         mode = scanMode, scan = scanNum,
                                         **scan_labels)
        
            # Should expand extract() call to include arguments to this.
            scan = deisotope_reduce_scan(scan, tolerance = 0.0001)  
            scan = [x for x in scan if x[0] > min_mz]
            assert charge, title
            writer.write(scan, title, mass = mz, charge = charge)
            scans_written += 1
        
    writer.close
    
    if inconsistent_precursors:
        print "Precursor inconsistencies: %s/%s" % (inconsistent_precursors,
                                                    scans_written)

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
    print "TEST MODE"
    #extract(r'C:\Users\Max\Desktop\Projects\rcprot_testing\2015-09-15-K562-10ng_1E5-TARG_50ms-mft_25-UFILL-1.raw',
            #r'C:\Users\Max\Desktop\Projects\rcprot_testing\2015-09-15-K562-10ng_1E5-TARG_50ms-mft_25-UFILL-1.raw.newExtractor.mgf')
    extract(r'\\rc-data1\blaise\ms_data_share\Max\CDK14\2017-03-01-CDK14-FMF-1.raw',
            r'\\rc-data1\blaise\ms_data_share\Max\CDK14\2017-03-01-CDK14-FMF-1.raw.newtest.mgf')