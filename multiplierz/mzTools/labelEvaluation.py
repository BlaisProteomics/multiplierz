from multiplierz.mzReport import reader, writer
from multiplierz.mzTools.featureDetector import detectFeatures, Feature
from multiplierz.mzTools.featureUtilities import FeatureInterface
import pickle
from collections import defaultdict
from multiplierz.mzAPI import mzFile




__all__ = ['evaluateMascotFile']

def evaluateSILAC(outputfile, columns, results, featureIntMap):
    peptides = defaultdict(list)
    for psm in results:
        if psm['Variable Modifications']:
            varmods = [x for x in psm['Variable Modifications'].split(';') if 'Label' not in x]
        else:
            varmods = ''
        peptides[psm['Peptide Sequence'], '; '.join(varmods)].append(psm)
    
    
    output = writer(outputfile, sheet_name = 'Data',
                    columns = columns + ['Lysines', 'Arginines', 'Lysine Labels',
                                         'Arginine Labels', 'Fully Labelled',
                                         'Intensity'])    
    
    evaluation = {'total peptides':0,
                  'total lysines':0,
                  'total arginines':0,
                  'fully labelled':0,
                  'lysine labelled':0,
                  'arginine labelled':0}    
    
    labelSummary = []
    for key, pepSet in list(peptides.items()):
        fullLabelledInt = 0
        partLabelledInt = 0
        
        seq = key[0]
        lysCount = seq.count('K')
        argCount = seq.count('R')        
        
        if not (lysCount + argCount):
            for psm in pepSet:
                psm['Lysines'] = 0
                psm['Arginines'] = 0
                psm['Lysine Labels'] = 0
                psm['Arginine Labels'] = 0
                psm['Fully Labelled'] = 'N/A'
                psm['Intensity'] = ''
                output.write(psm)
            
            continue
        
        for psm in pepSet:
            assert seq == psm['Peptide Sequence']
            if psm['Variable Modifications']:
                varmods = psm['Variable Modifications'].split(';')
            else:
                varmods = ''
            
            lysLabelCount = len([x for x in varmods if x[0] == 'K' and 'Label' in x])
            argLabelCount = len([x for x in varmods if x[0] == 'R' and 'Label' in x])
            
            psm['Lysines'] = lysCount
            psm['Arginines'] = argCount
            psm['Lysine Labels'] = lysLabelCount
            psm['Arginine Labels'] = argLabelCount
            
            try:
                intensity = featureIntMap[int(psm['Spectrum Description'].split('.')[1])]
            except KeyError:
                intensity = 0

            if lysLabelCount == lysCount and argCount == argLabelCount:
                fullLabelledInt += intensity
                psm['Fully Labelled'] = True
            else:
                partLabelledInt += intensity
                psm['Fully Labelled'] = False
            
            psm['Intensity'] = intensity
            output.write(psm)
            
            evaluation['total peptides'] += 1
            evaluation['total lysines'] += lysCount
            evaluation['total arginines'] += argCount
            evaluation['fully labelled'] += int(psm['Fully Labelled'])
            evaluation['lysine labelled'] += lysLabelCount
            evaluation['arginine labelled'] += argLabelCount
        
        labelSummary.append((key, fullLabelledInt, partLabelledInt))
    output.close()
    
    summaryOutput = writer(outputfile, sheet_name = 'Label Report',
                           columns = ['Peptide Sequence', 'Non-Label Mods',
                                      '% Labelling'])
    for (seq, varmods), fullLabelledInt, partLabelledInt in labelSummary:
        row = {}
        row['Peptide Sequence'] = seq
        row['Non-Label Mods'] = varmods
        if partLabelledInt + fullLabelledInt:
            row['% Labelling'] = fullLabelledInt / (partLabelledInt + fullLabelledInt)
        else:
            row['% Labelling'] = ''
        summaryOutput.write(row)
    
    summaryOutput.close()
    
    return evaluation
    
    
def evaluateTMTiTRAQ(outputfile, columns, results, resultIntMap):
    peptides = defaultdict(list)
    for psm in results:
        try:
            varmods = [x for x in psm['Variable Modifications'].split(';') if 'plex' not in x]
        except AttributeError:
            varmods = ''
        peptides[psm['Peptide Sequence'], '; '.join(varmods)].append(psm)
    
    output = writer(outputfile, sheet_name = 'Data',
                    columns = columns + ['Lysines','Lysine Labels','N-term Label',
                                         'Fully Labelled', 'Intensity'])
    
    evaluation = {'total peptides':0,
                  'total lysines':0,
                  'fully labelled':0,
                  'nterm labelled':0,
                  'lysine labelled':0}
    
    labelSummary = []
    for key, pepSet in list(peptides.items()):
        fullLabelledInt = 0
        partLabelledInt = 0
        
        seq = key[0]
        lysCount = seq.count('K')
        
        for psm in pepSet:
            assert seq == psm['Peptide Sequence']
            
            # Should work now?
            lysLabelCount = len([x for x in psm['Variable Modifications'].split('; ')
                                 if x and x[0] == 'K' and 'plex' in x])
            ntermLabel = any([(x[:len('N-term')] == 'N-term' and 'plex' in x) for
                              x in psm['Variable Modifications'].split('; ')])
            
            psm['Lysines'] = lysCount
            psm['Lysine Labels'] = lysCount
            psm['N-term Label'] = ntermLabel
            
            try:
                intensity = resultIntMap[int(psm['Spectrum Description'].split('.')[1])]
            except KeyError:
                intensity = 0
            
            if (lysLabelCount == lysCount) and ntermLabel:
                fullLabelledInt += intensity
                psm['Fully Labelled'] = True
            else:
                partLabelledInt += intensity
                psm['Fully Labelled'] = False
            
            psm['Intensity'] = intensity
            output.write(psm)
            
            evaluation['total peptides'] += 1
            evaluation['total lysines'] += lysCount
            evaluation['fully labelled'] += 1 if not partLabelledInt else 0
            evaluation['nterm labelled'] += int(ntermLabel)
            evaluation['lysine labelled'] += lysLabelCount
        
        labelSummary.append((key, fullLabelledInt, partLabelledInt))
    output.close()
    
    summaryOutput = writer(outputfile, sheet_name = 'Label Report',
                           columns = ['Peptide Sequence', 'Non-Label Mods',
                                      '% Labelling'])
    for (seq, varmods), fullLabelledInt, partLabelledInt in labelSummary:
        row = {}
        row['Peptide Sequence'] = seq
        row['Non-Label Mods'] = varmods
        if partLabelledInt + fullLabelledInt:
            row['% Labelling'] = fullLabelledInt / (partLabelledInt + fullLabelledInt)
        else:
            row['% Labelling'] = ''
        summaryOutput.write(row)
    
    summaryOutput.close()         
    
    return evaluation
        

def evaluateMascotFile(resultfile, datafile = None, featurefile = None, outputfile = None):
    #assert datafile or featurefile, "Either raw data or feature data must be given!"
    
    header = [list(x.values()) for x in list(reader(resultfile, sheet_name = 'Mascot_Header'))]
    
    def retrieveHeaderValue(key):
        try:
            return [[x for x in xs if x != key] for xs in header if key in xs][0][0]
        except IndexError:
            return ''
    quant = retrieveHeaderValue('Quantitation method')
    varmods = retrieveHeaderValue('Variable modifications')
    
    assert ('SILAC' in quant) or ('plex' in varmods), "Label method not recognized!"
    
    if not featurefile:
        featurefile = detectFeatures(datafile, signalToNoiseThreshold = 15)
        features = FeatureInterface(featurefile)
    else:
        features = FeatureInterface(featurefile)
    
    print("Matching features to PSMs...")
    results = reader(resultfile)
    columns = results.columns
    results = list(results)
    
    
    data = mzFile(datafile)
    ms1map = {}
    ms2s = []
    ms1 = None
    for _, _, scan, level, _ in data.scan_info(0, 999999):
        if level == 'MS1':
            for ms2 in ms2s:
                ms1map[ms2] = ms1
            ms1 = scan
            ms2s = []
        elif level == 'MS2':
            ms2s.append(scan)
    ms1map[ms1] = ms2s
    data.close()
    
    featureIntMap = {}
    for psm in results:
        mz = psm['Experimental mz']
        scan = int(psm['Spectrum Description'].split('.')[1])
        charge = int(psm['Charge'])
        for index, feature in features.mz_range(mz - 1, mz + 1):
            if feature.containsPoint(mz, ms1map[scan], charge):
                featureIntMap[scan] = feature.c12Intensity()
                break
    
    del features    
    
    if not outputfile:
        outputfile = '.'.join(resultfile.split('.')[:-1]) + '_LABEL_EVALUATION.xlsx'
    
    if 'SILAC' in quant:
        return evaluateSILAC(outputfile, columns, results, featureIntMap), outputfile
    elif 'plex' in varmods:
        return evaluateTMTiTRAQ(outputfile, columns, results, featureIntMap), outputfile




if __name__ == '__main__':
    evaluateMascotFile(r'\\rc-data1\blaise\ms_data_share\Max\LabelEvaluator\2017-10-11-Buhrlage-DUB-Set1_HCD_RECAL.mgf.xlsx',
                       r'\\rc-data1\blaise\ms_data_share\Max\LabelEvaluator\2017-10-11-Buhrlage-DUB-Set1.raw')