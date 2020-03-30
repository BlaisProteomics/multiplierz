from collections import defaultdict
from multiplierz.internalAlgorithms import insert_tag, collectByCriterion
import os

def calculate_FDR(reportfile, outputfile = None, threshold = 0.01,
                  decoyString = 'rev_', includeStatisticsSheet = True,
                  includeDuplicates = True, separateDuplicateSheet = True,
                  includeFailedSheet = True, includeReverseSheet = True,
                  single_cutoff = True):
    """
    Performs Forward/Reverse database filtering on the target file, giving back
    the true PSMs over the specified statistical threshold as well as removed decoy
    and below-threshold PSMs, in respective sheets.

    All entries in the decoy (reverse) database must have accessions that begin with
    some uniform prefix; by default, "rev_" (so that gi|198292342|X7823_EXTRA becomes
    rev_gi|198292342|X7823_EXTRA.)
    
    outputfile may be safely specified to be the same as the input file, in
    order to overwrite the original file.
    """

    from multiplierz.mzReport import reader, writer

    reportReader = reader(reportfile)
    reportRows = list(reportReader)
    columns = reportReader.columns + ['FDR']
    reportReader.close()

    reportRows.sort(key = lambda x: x['Peptide Score'], reverse = True)

    seenSpectra = {}

    passedRows = []
    failedRows = []
    duplicateRows = []
    reverseRows = []

    reverses = 0.0
    forwards = 0.0
    duplicates = 0
    passed = 0
    failed = 0
    lowPass = 999999999
    highRev = 0
    for row in reportRows:
        specDesc = row['Spectrum Description']
        if specDesc in seenSpectra:
            duplicates += 1
            fdr = seenSpectra[specDesc]
            row['FDR'] = fdr
            if includeDuplicates and not separateDuplicateSheet:
                if fdr < threshold:
                    passedRows.append(row)
                else:
                    failedRows.append(row)
            else:
                duplicateRows.append(row)
            continue

        #if decoyString in row['Accession Number'].lower():
        # Turns out that produced awful results, since high-scoring peptides
        # could just happen to be duplicated in the reverse database.
        # So instead:
        if all([decoyString in x.lower() for x in row['Accession Number'].split(';')]):
            reverses += 1
            if forwards:
                fdr = reverses / forwards
            else:
                fdr = 100
            row['FDR'] = fdr      

            if float(row['Peptide Score']) > highRev:
                highRev = float(row['Peptide Score'])

            seenSpectra[specDesc] = fdr
            reverseRows.append(row)
        else:
            forwards += 1
            fdr = reverses / forwards
            row['FDR'] = fdr            

            seenSpectra[specDesc] = fdr
            if fdr < threshold:
                passed += 1
                passedRows.append(row)
                if float(row['Peptide Score']) < lowPass:
                    lowPass = float(row['Peptide Score'])
            else:
                failed += 1
                failedRows.append(row)

    if single_cutoff:
        recovered = [x for x in failedRows if x['Peptide Score'] > lowPass]
        failedRows = [x for x in failedRows if x['Peptide Score'] <= lowPass]
        passedRows += recovered

    if not outputfile: 
        # Output format must support sheets.
        if reportfile.lower().endswith('xlsx') or reportfile.lower().endswith('xls'):
            outputfile = insert_tag(reportfile, 'FDR_filtered')
        else:
            outputfile = '.'.join(reportfile.split('.')[:-1] + ['FDR_filtered.xlsx'])
        is_excel_output = True
    else:
        ext = outputfile.lower().split('.')[-1]
        is_excel_output = ext in {'xls', 'xlsx'}
        
    percentage = round(threshold * 100)


    if includeStatisticsSheet and is_excel_output:
        statOutput = writer(outputfile, columns = ['FDR Calculation Statistics', '--------------'],
                            sheet_name = "FDR Statistics")
        statOutput.write(['', ''])
        statOutput.write(['Total Spectra', str(len(reportRows))])
        statOutput.write(['Passed %s%% FDR' % percentage, str(passed)])
        statOutput.write(['Lowest Passing Score', str(lowPass)])
        statOutput.write(['Reverse Hits', str(reverses)])
        statOutput.write(['Highest Scoring Reverse Hit', str(highRev)])
        statOutput.write(['Number of Duplicates', str(duplicates)])
        statOutput.close()

    if includeFailedSheet and is_excel_output:
        failedOutput = writer(outputfile, columns = columns,
                              sheet_name = "Failed %s%% FDR" % percentage)
        for row in failedRows:
            failedOutput.write(row)
        failedOutput.close()

    if separateDuplicateSheet and is_excel_output:
        duplicateOutput = writer(outputfile, columns = columns,
                                 sheet_name = "Duplicate Rows")
        for row in duplicateRows:
            duplicateOutput.write(row)
        duplicateOutput.close()

    if includeReverseSheet and is_excel_output:
        reverseOutput = writer(outputfile, columns = columns,
                               sheet_name = 'Reverse Hits')
        for row in reverseRows:
            reverseOutput.write(row)
        reverseOutput.close()

    passedOutput = writer(outputfile, columns = columns, 
                          sheet_name = "Data" if is_excel_output else None)
    for row in passedRows:
        passedOutput.write(row)
    passedOutput.close()   

    return outputfile

def combine_accessions(reportfile, outputfile = None):
    """
    Given a Mascot-style PSM report, this combines all protein hypotheses for a given
    MS2 spectrum into a single PSM.
    
    outputfile may be safely specified to be the same as the input file, in
    order to overwrite the original file.
    """
    
    
    from multiplierz.mzReport import reader, writer
    
    report = reader(reportfile)
    columns = report.columns

    molecules = defaultdict(list)
    for row in report:
        molecules[row['Spectrum Description']].append(row)
        
    
    outputData = []
    for rows in list(molecules.values()):
        accessions = [x['Accession Number'] for x in rows]
        newRow = max(rows, key = lambda x: x['Peptide Score'])
        
        if 'Accession Number' in columns:
            newRow['Accession Number'] = '; '.join([x['Accession Number'] for x in rows])
        if 'Protein Description' in columns:
            newRow['Protein Description'] = '; '.join([x['Protein Description'] for x in rows])
        if 'Protein Masses' in columns:
            newRow['Protein Masses'] = '; '.join([str(x['Protein Mass']) for x in rows])
        newRow['Protein Redundancy'] = len(rows)
        outputData.append(newRow)
    
    try:
        columns = [x for x in columns + ['Protein Masses'] if x in newRow]
    except UnboundLocalError:
        pass # Means there was no newRow, and thus no rows, so it's pretty arbitrary.
    
    if not outputfile:
        outputfile = insert_tag(reportfile, 'combined_accessions')
    
    output = writer(outputfile, columns = columns + ['Protein Redundancy'])
    report.close()
    for row in outputData:
        output.write(row)
    output.close()
    
    
    
    return outputfile


def concatenate_reports(reportfiles, outputfile,
                        include_file_column = False):
    from multiplierz.mzReport import reader, writer
    readers = list(map(reader, reportfiles))
    allcols = sorted(set.intersection(*list(map(set, [x.columns for x in readers]))),
                     key = lambda x: readers[0].columns.index(x))
    if not all(x.columns == allcols for x in readers):
        print("Warning- concatenation drops some columns!")
    output = writer(outputfile, columns = allcols)
    for filename, report in zip(reportfiles, readers):
        filename = os.path.basename(filename)
        for row in report:
            if include_file_column:
                row['FILE'] = filename
            output.write(row, ignore_extra = True)
    output.close()
    return outputfile





def combine_peptides(reportfile, isobaric = None, outputfile = None):    
    from multiplierz.mzReport import reader, writer
    from multiplierz.mgf import standard_title_parse
    
    isobaric_labels = {None: [],
                       4: ['114', '115', '116', '117'],
                       6: ['126', '127', '128', '129', '130', '131'],
                       8: ['113', '114', '115', '116', '117', '118', '119', '121'],
                       10: ['126', '127N', '127C', '128N', '128C', 
                            '129N', '129C', '130N', '130C', '131']}    
    
    def _byPeptide(row):
        # Not counting charge.
        varmodset = frozenset([x.strip() for x in 
                               row['Variable Modifications'].split(';')])
        return row['Peptide Sequence'], varmodset    
    
    def _getReporters(row):
        attrib = standard_title_parse(row['Spectrum Description'])
        return [float(attrib[x.lower()]) for x in isobaric_labels[isobaric]]        
    
    assert isobaric in isobaric_labels
    
    psms = reader(reportfile)
    rowsByPeptide = collectByCriterion(psms, _byPeptide)
    
    
    sum_cols = ['Sum%s' % x for x in isobaric_labels[isobaric]]
    top_cols = ['Max%s' % x for x in isobaric_labels[isobaric]]
    if not outputfile:
        outputfile = insert_tag(reportfile, 'peptide_combined')
    output = writer(outputfile,
                    columns = (psms.columns + sum_cols + top_cols + ['PSMs']))
    
    for pep, psms in rowsByPeptide.items():
        outrow = max(psms, key = lambda x: x['Peptide Score'])
        outrow['PSMs'] = len(psms)
        
        if isobaric:
            repsets = [_getReporters(x) for x in psms]
            toprepset = max(repsets, key = lambda x: sum(x))
            sumrepset = [sum(x) for x in zip(*repsets)]
            
            for rep, col in zip(toprepset, top_cols):
                outrow[col] = rep
            for rep, col in zip(sumrepset, sum_cols):
                outrow[col] = rep
        
        output.write(outrow)
    
    output.close()
    
    return outputfile        
    
    


def fractionation_plot(fractions, outputfile = None, fig_size = None, **kwargs):
    """
    Takes a list of 3-tuples (<organic fraction>, <salt fraction>, <filename>)
    describing the PSM output of a multi-fraction MS experiment; draws a
    plot where each fraction is represented as an appropriately scaled point
    (according to PSM count) in organic/salt space.
    """
    from multiplierz.mzReport import reader
    import matplotlib.pyplot as pyt
    pyt.cla()
    
    fractions = [(float(o), float(s), f) for o, s, f in fractions]
    
    organics = sorted(set(zip(*fractions)[0]))
    salts = sorted(set(zip(*fractions)[1]))
    
    orgCoords = dict([(o, i) for i, o in enumerate(organics, start = 1)])
    saltCoords = dict([(s, i) for i, s in enumerate(salts, start = 1)])
    
    
    if fig_size:
        fig = pyt.gcf()
        fig.set_size_inches(*fig_size)
    elif len(orgCoords) > 8 or len(saltCoords) > 8:
        fig = pyt.gcf()
        cursize = fig.get_size_inches()
        newsize = [cursize[0], cursize[1]]
        if len(orgCoords) > 8:
            cursize[0] = max(cursize[0], len(orgCoords) * 0.9)
        if len(saltCoords) > 8:
            cursize[1] = max(cursize[1], len(saltCoords) * 0.9)
        fig.set_size_inches(*cursize)
        
        
        
    
    scatterPts = []
    for organic, salt, psms in fractions:
        orgcoord = orgCoords[organic]
        saltcoord = saltCoords[salt]
        if isinstance(psms, int): # Can just pass the count.
            count = psms
        elif isinstance(psms, str): # Else the file
            rdr = reader(psms)
            try:
                count = rdr.get_row_count()
            except (IOError, AttributeError):
                count = len(list(rdr))
            rdr.close()
        else:
            raise Exception("Must specify PSM count or file.")
        
        scatterPts.append((orgcoord, saltcoord, count))
        pyt.text(orgcoord, saltcoord, str(count),
                 verticalalignment = 'center',
                 horizontalalignment = 'center')
        
    orgRange = max(orgCoords.values()) - min(orgCoords.values())
    saltRange = max(saltCoords.values()) - min(saltCoords.values())
    orgMargin = orgRange / 15.0
    saltMargin = saltRange / 15.0
    overallRange = min(orgRange, saltRange)
        
    ax = pyt.axes()
    ax.set_xlim(min(orgCoords.values()) - orgMargin, max(orgCoords.values()) + orgMargin)
    ax.set_ylim(min(saltCoords.values()) - saltMargin, max(saltCoords.values()) + saltMargin)
    
        
    #def count_to_size(counts):
        #counts / 
        
    orgpts, saltpts, counts = list(zip(*scatterPts))
    pyt.scatter(orgpts, saltpts, counts,
                alpha = 0.2,
                **kwargs)
    
    orgTicks = [(v, k) for k, v in list(orgCoords.items())]
    saltTicks = [(v, k) for k, v in list(saltCoords.items())]
    pyt.xticks(*list(zip(*orgTicks)))
    pyt.yticks(*list(zip(*saltTicks)))
    pyt.xlabel('Organic')
    pyt.ylabel('Salt')
    
    
    
    #pyt.xlim(min(orgCoords.values()) - 0.5, max(orgCoords.values()) + 0.5)
    #pyt.ylim(min(saltCoords.values()) - 0.5, max(saltCoords.values()) + 0.5)
    ##pyt.set_aspect('equal')    
    #pyt.autoscale(enable = False)
    
    #print pyt.xlim(), pyt.ylim()
    if not outputfile:
        pyt.show()
    else:
        pyt.savefig(outputfile)
    
    #print pyt.xlim(), pyt.ylim()
    pyt.cla()
    



def multimode_fractionation_plot(mode_fractions, outputfile = None,
                                 count_to_size = (lambda x: x/100),
                                 fig_size = None,
                                 color_sequence = None):
    """
    mode_fractions must be a list of (<mode title>, <fractions>) tuples where
    <fractions> is a list of (<organic fraction>, <salt fraction>,
    <filename>) tuples as in the input for fractionation_plot.  Output is a plot
    of each fraction in organic/salt space, where each point is a pie chart
    showing absolute magnitude of all modes in the space (via chart size) and
    relative magnitude of each mode (via chart slice sizes.)
    """


    from multiplierz.mzReport import reader
    import matplotlib.pyplot as pyt    
    
    pyt.cla()
    if fig_size:
        fig = pyt.gcf()
        fig.set_size_inches(*fig_size)    
    
    for i in range(len(mode_fractions)):
        mode_fractions[i] = mode_fractions[i][0], [(float(o), float(s), f) for o, s, f 
                                                   in mode_fractions[i][1]]
    
    modes = list(zip(*mode_fractions))[0]
    organics = sorted(set(sum([list(zip(*x[1])[0]) for x in mode_fractions], [])))
    salts = sorted(set(sum([list(zip(*x[1])[1]) for x in mode_fractions], [])))	
    
    orgCoords = dict([(o, i) for i, o in enumerate(organics, start = 1)])
    saltCoords = dict([(s, i) for i, s in enumerate(salts, start = 1)])    
    

    
    grid = defaultdict(dict)
    largest = 0
    for mode, fractions in mode_fractions:
        for organic, salt, filename in fractions:
            orgCoord = orgCoords[organic]
            saltCoord = saltCoords[salt]
            rdr = reader(filename)
            count = len(list(rdr))
            rdr.close()
            grid[orgCoord, saltCoord][mode] = count
            largest = max([largest, count])
    
    xScale = (max(list(zip(*list(grid.keys())))[0]) - min(list(zip(*list(grid.keys())))[0])) / float(len(orgCoords))
    yScale = (max(list(zip(*list(grid.keys())))[1]) - min(list(zip(*list(grid.keys())))[1])) / float(len(saltCoords))
    overallscale = min([xScale, yScale])    
    
    #def countConvert(count):
        #return (float(count) / (largest*3)) * float(overallscale)
    def countConvert(x): return 1
    
    fig  = pyt.figure()
    ax = fig.gca()
    for (orgCoord, saltCoord), modecounts in list(grid.items()):
        countForAllModes = [modecounts.get(m, 0) for m in modes]
        total = sum(modecounts.values())
        ax.pie(countForAllModes, center = (orgCoord, saltCoord),
               radius = countConvert(total), frame = True,
               colors = color_sequence)
        

    
        
    ax.set_xlim(min(orgCoords.values()) - 0.5, max(orgCoords.values()) + 0.5)
    ax.set_ylim(min(saltCoords.values()) - 0.5, max(saltCoords.values()) + 0.5)
    ax.set_aspect('equal')
    
    ax.legend(modes, loc = 'upper left')
    
    if not outputfile:
        pyt.show()
    else:
        pyt.savefig(outputfile)
    pyt.cla()
        
        
    
