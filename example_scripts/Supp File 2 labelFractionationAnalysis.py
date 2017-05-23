from multiplierz.mzReport import reader, writer
from multiplierz.mzAPI import mzFile
from multiplierz.mgf import extract
from multiplierz.mzSearch import MascotSearch
from multiplierz.post_process import combine_accessions, calculate_FDR
from multiplierz.post_process import fractionation_plot, multimode_fractionation_plot
import matplotlib.pyplot as pyt
import os
from collections import defaultdict
from multiplierz.internalAlgorithms import collectByCriterion
from numpy import average
import time



# If running this script with a copy of Mascot that has
# security enabled, write in the username and password
# of the account to submit the search here.
MASCOT_USERNAME = '...'
MASCOT_PASSWORD = '...'

# Specify the source directory containing the MS data files here.
# Results will be written into the same directory.
directory = r'C:\my_data_directory'

# This list indicates the (un)labelled states that will be examined.
# Each entry should be a three-tuple; the first element is a descriptive
# name of the label type, the second element is the base of the name
# of the subdirectory where search results with each labelling are stored,
# and the third element is a parameter file proper for searching
# with the given label type.
label_state_subdirs = [('Unlabelled', '0plex', 'HumanCarbOx-Generic.par'),
                       ('TMT 4-plex', '4plex', 'HumanCarbOx-4plexCold.par'),
                       ('iTRAQ 6-plex', '6plex', 'HumanCarbOx-6plexCold.par'),
                       ('TMT 8-plex', '8plex', 'HumanCarbOx-8plex.par')]



def typeInDir(directory, ext):
    import os
    return [os.path.join(directory, x) for x in os.listdir(directory) if x.lower().endswith(ext.lower())]

def parseFractionFilename(filename):
    """
    Parses the filenames of each RAW file from the run to determine
    the organic/salt fraction represented.
    """
    procname = filename.replace('Ac', '900')
    return procname.split('.')[0].replace('p','.').split('-')[-2:] + [filename]

def extractionProcess(directory):
    """
    Manages the RAW-to-MGF extraction function, which converts MS2 data from
    the RAW file into a set of spectra in Mascot Generic Format.
    """
    raws = typeInDir(directory, 'raw')
    mgffiles = []
    for datafile in raws:
        mgffiles += extract(datafile)
        
    return mgffiles

def searchProcess(targetDir, parfile, mgfs):
    """
    Manages the interface to the Mascot search server.
    
    In the call to MascotSearch.run_search(), "user" and "password" must match
    a valid user for the server specified in the parameter file (or omitted
    if security is not enabled on that server.)
    """
    searcher = MascotSearch(parfile)
    outputfiles = []
    for mgffile in mgfs:
        try:
            outputfile = searcher.run_search(mgffile,
                                             user = MASCOT_USERNAME,
                                             password = MASCOT_PASSWORD)
        except Exception as err:
            print "---------------- Failed on %s" % mgffile
            raise err
        targetoutput = os.path.join(targetDir, os.path.basename(outputfile))
        os.rename(outputfile, targetoutput)
        outputfiles.append(targetoutput)
    
    return outputfiles
    
    

def postProcessingSteps(targetDir):
    """
    For every search result file, two important post-processing steps are
    performed;
    - combine_accessions combines redundant peptide hits produced by Mascot
    - calculate_FDR derives the False Discovery Rate based on the rate of matches
    reversed sequences embedded in the sequence database, and filters out
    PSMs that exceed 1% FDR.
    """    
    
    for dirpath, dirnames, filenames in os.walk(targetDir):
        for filename in filenames:
            filename = os.path.join(dirpath, filename)
            if filename.lower().endswith('.xls') or filename.lower().endswith('.xlsx'):
                if 'FDR' not in filename and 'COM' not in filename:
                    combinedfilename = '.'.join(filename.split('.')[:-1] + ['COM', 'xlsx'])
                    fdrfilename = '.'.join(combinedfilename.split('.')[:-1] + ['FDR', 'xlsx'])
                    
                    combine_accessions(filename, combinedfilename)
                    calculate_FDR(filename, fdrfilename)

  
  
    
    

def peptideKey(psm):
    """
    Utility function used in psm_intersection(), used to help group PSMs by
    the peptide and modification state they represent.
    """
    return (psm['Peptide Sequence'],
            frozenset(map(lambda y: y.strip(),
                          psm['Variable Modifications'].split('; '))))

def psm_intersection(directory, mode_subdirs):
    """
    To give a more accurate depiction of the relative elution profile of each
    label state, the final results will only consider peptides that appear in
    the results for all four states. This determines the overlapping peptide
    repertoire detected across all four experiments, and produces subset
    result files that only include these peptides.
    """
    
    psmByCondition = defaultdict(list)
    for mode, subdir, par in mode_subdirs:
        files = typeInDir(os.path.join(directory, subdir), 'xlsx')
        conditionPSMs = []
        for resultfile in files:
            if not 'FDR' in resultfile:
                continue
            conditionPSMs += list(reader(resultfile))
        psmByCondition[subdir] = collectByCriterion(conditionPSMs, peptideKey)
    
    consistentPSMs = reduce(set.intersection, [set(x.keys()) for x in psmByCondition.values()],
                            set(psmByCondition.values()[0].keys()))
    
    newSubdirs = []
    for mode, subdir, par in mode_subdirs:
        newSubdir = subdir + '_intersection_sheets'
        newSubdirs.append((mode, newSubdir))
        try:
            os.mkdir(os.path.join(directory, newSubdir))
        except:
            pass
        files = typeInDir(os.path.join(directory, subdir), 'xlsx')
        for filename in files:
            alreadySeenPeptides = set()
            if not 'FDR' in filename:
                continue            
            psms = reader(filename)
            filterfile = writer(os.path.join(directory, newSubdir, os.path.basename(filename)),
                                columns = psms.columns)
            for psm in psms:
                pepKey = peptideKey(psm)
                if pepKey in consistentPSMs and pepKey not in alreadySeenPeptides:
                    alreadySeenPeptides.add(pepKey)
                    filterfile.write(psm)
            filterfile.close()
            
    return newSubdirs
        
        
    
        
    
def psm_XIC_localized(directory, subdirs):  
    """
    A peptide may appear in multiple fractions due various factors, but for
    the purpose of this analysis it is useful to consider a peptide as
    "belonging" only to the fraction in which the main bulk of the elution
    occurred. For each fraction in which a given peptide appeared, we take
    XICs over the m/z values for a set of possible charge and compare their
    total intensity; the fraction with the most intense XIC(s) is assigned
    that peptide for the final count.
    """
    
    tolerance = 0.1
    time_tolerance = 15
    
    rawfiles = dict([(x.split('.')[0], mzFile(os.path.join(directory, x)))
                     for x in os.listdir(directory) if x.lower().endswith('raw')])
    columns = None
    
    start = time.clock()
    for subdir in subdirs:
        resultfiles = typeInDir(os.path.join(directory, subdir), 'xlsx')
        resultfiles = [x for x in resultfiles if 'XIC_localized' not in x]
        
        peptidesForFile = defaultdict(dict)
        for resultfile in resultfiles:
            rdr = reader(resultfile)
            columns = rdr.columns
            psmsByPeptide = collectByCriterion(list(rdr),
                                               lambda x: (x['Peptide Sequence'],
                                                          x['Variable Modifications']))
            for peptide, psms in psmsByPeptide.items():
                peptidesForFile[peptide][resultfile] = psms
        
        outputByFile = defaultdict(list)
        for peptide, psmsByFile in peptidesForFile.items():
            xicsByFile = []
            
            allPSMs = sum(psmsByFile.values(), [])
            mass = allPSMs[0]['Predicted mr']
            assert len(set(x['Predicted mr'] for x in allPSMs)) == 1
            
            charges = set(x['Charge'] for x in allPSMs)
            allScans = set([tuple(x['Spectrum Description'].split('.')[:2]) for x in allPSMs])
            allRTs = set(rawfiles[x[0]].scan_time_from_scan_name(int(x[1]))
                          for x in allScans)
            minRT, maxRT = min(allRTs), max(allRTs)
            
            for resultfile, psms in psmsByFile.items():
                rawfile = rawfiles[os.path.basename(resultfile.split('.')[0])]
                xicInt = 0
                for charge in charges:
                    mz = (mass + (1.0072764 * charge)) / charge
                    xic = rawfile.xic(minRT - time_tolerance, maxRT + time_tolerance,
                                      mz - tolerance, mz + tolerance)
                    xicInt += sum(zip(*xic)[1])
                
                xicsByFile.append((xicInt, resultfile))
            
            highIntFile = max(xicsByFile, key = lambda x: x[0])[1]
            outputByFile[highIntFile].append(psmsByFile[highIntFile][0])
        
        for resultfile, psms in outputByFile.items():
            outputfile = resultfile[:-5] + '.XIC_localized.xlsx'
            output = writer(outputfile, columns = columns)
            for psm in psms:
                output.write(psm)
            output.close()
        
        
        
if __name__ == '__main__':
    # See docstrings for each function for an explanation of each step.
    
    time.clock()
    
    mgfs = extractionProcess(directory)
    print "MGFs extracted: %s" % time.clock()
    
    
    for mode, subdir, parfile in label_state_subdirs:
        state_dir = os.path.join(directory, subdir)
        parfile = os.path.join(directory, parfile)
        searchProcess(state_dir, parfile, mgfs)
    print "Searches completed: %s" % time.clock()
    
    postProcessingSteps(directory)
    print "Post-processing completed: %s" % time.clock()
    
    intersection_subdirs = psm_intersection(directory, label_state_subdirs)
    print "PSM Intersection report completed: %s" % time.clock()
    
    psm_XIC_localized(directory, zip(*intersection_subdirs)[1])
    
    all_mode_fractions = []
    for mode, subdir in intersection_subdirs:
        mode_dir = os.path.join(directory, subdir)
        files = typeInDir(mode_dir, '.XIC_localized.xlsx')
        filesWithFraction = map(parseFractionFilename, files)
        
        fractionation_plot(filesWithFraction,
                           os.path.join(mode_dir, '%s_fractionation_plot.svg' % mode))
        all_mode_fractions.append((mode, filesWithFraction))
    
    multimode_fractionation_plot(all_mode_fractions,
                                 os.path.join(directory,
                                              'all_mode_fractionation_plot.svg'))
    
        
    
    print "Analysis completed: %s" % time.clock()    
    
    print "Done!"
