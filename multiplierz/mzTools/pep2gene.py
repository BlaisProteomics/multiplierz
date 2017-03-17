from multiplierz.mzReport import reader, writer
from multiplierz.fasta import parse_to_dict
from multiplierz.settings import settings

from multiplierz.internalAlgorithms import print_progress
import cPickle as pickle
import time
import os
import re
from collections import defaultdict
import urllib, urllib2

K = 4


def convertAccessionsViaUniprot(accessions):
    uniprotReg = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    ncbiReg = '[A-Z]{2}_[0-9]*.*[0-9]*'

    uniprotAccs = all([re.match(uniprotReg, x) for x in accessions[:10]])
    ncbiAccs = all([re.match(ncbiReg, x) for x in accessions[:10]])

    assert uniprotAccs or ncbiAccs, "Accession format not recognized! %s" % accessions[0]

    if uniprotAccs:
        fromKey = 'ACC+ID'
        accessions = [x.split('-')[0] for x in accessions]
    else:
        #accessions = [x.split('.')[0] for x in accessions if 'P_' in x]
        accessions = [x.split('.')[0] for x in accessions]
        fromKey = 'P_REFSEQ_AC'    

    url = 'http://www.uniprot.org/uploadlists/'
    
    geneNames = defaultdict(list)
    for i in range(0, len(accessions), 1000):
        subAccList = accessions[i:i+1000]
        
        geneNameReq = {
            'from':fromKey,
            'to':'ACC',
            'format':'tab',
            'query':' '.join(subAccList),
        }    
        
        contact = settings.user_email    
    
        data = urllib.urlencode(geneNameReq)
        request = urllib2.Request(url, data)
        request.add_header('User-Agent', 'Multiplierz: %s' % contact)
        response = urllib2.urlopen(request)
        geneNamePage = response.read()  
        
        table = [x.split('\t') for x in geneNamePage.split('\n')[1:] if x]
        for row in table:
            acc = row[0]
            # The gene field seems to be a list of synonymous names,
            # which should be considered as a unit.
            geneNames[acc].append(row[6])
            #uniprotIDs[acc] = row[3]
    
    normedGeneNames = dict()
    for acc, genes in geneNames.items():
        normedGeneNames[acc.split('.')[0]] = (', '.join(genes), '')
    
    return normedGeneNames#, uniprotIDs


def create_fasta_index(fasta_file, outputfile, labelParser = (lambda x: x),
                       distinguish_leucine = True):
    starttime = time.clock()
    
    if not outputfile:
        outputfile = fasta_file + '.pep2gene'
    
    if isinstance(labelParser, basestring):
        # If called from the GUI, a function can't be passed for
        # the sake of async.  But presumably the string is a regex.
        parser = re.compile(labelParser)
        def labelParser(label):
            parsed = parser.search(label)
            if not parsed:
                raise IOError, "Could not parse FASTA label: %s" % label
            return parsed.group(0)  
    
    prevtime = time.clock()
    fasta = parse_to_dict(fasta_file, labelConverter = labelParser,
                          filter_string = 'rev_')
    print "FASTA file read: %.2f" % (time.clock() - prevtime)
    prevtime = time.clock()
    
    prot_to_gene = convertAccessionsViaUniprot(fasta.keys())
    print "Gene lookup acquired: %.2f" % (time.clock() - prevtime)
    prevtime = time.clock()
    
    if distinguish_leucine:
        isofasta = {}
        for prot, seq in fasta.items():
            isofasta[prot] = seq.replace('I', 'L')
    else:
        isofasta = None
        
    fmerLookup = defaultdict(set)
    isomerLookup = defaultdict(set)
    for protein, seq in fasta.items():
        for i in range(len(seq) - 3):
            fmer = seq[i:i+K]
            fmerLookup[fmer].add(protein)
            if distinguish_leucine:
                isomerLookup[fmer.replace('I', 'L')].add(protein)
    
    print "Lookup tables created: %.2f" % (time.clock() - prevtime)
    prevtime = time.clock()
    
    out = open(outputfile, 'wb')
    data = (K, fasta, fmerLookup, prot_to_gene, isofasta, isomerLookup)
    pickle.dump(data, out, protocol = 2)
    out.close()
    
    print "Data written: %.2f" % (time.clock() - prevtime)
    print "Overall database generation time: %.2f" % (time.clock() - starttime)
    
    return outputfile
    
    
def add_gene_ids(target_file, p2g_database,
                 target_sheet = None,
                 outputfile = None, inPlace = False, distinguish_leucine = True,
                 legacy_columns = True):
    starttime = time.clock()
    
    dataRdr = open(p2g_database, 'rb')
    data = pickle.load(dataRdr)
    if isinstance(data, tuple) and len(data) == 6:
        k_len, seqLookup, fmerLookup, geneLookup, isoSeqLookup, isoFmerLookup = data
    elif isinstance(data, tuple) and not len(data) == 6:
        raise Exception, str(len(data))
    else:
        print 'Legacy mode P2G database detected!'
        seqLookup = data
        fmerLookup = pickle.load(dataRdr)
        geneLookup = pickle.load(dataRdr)
        isoSeqLookup = pickle.load(dataRdr)
        isoFmerLookup = pickle.load(dataRdr)
    dataRdr.close()
    
    if distinguish_leucine:
        assert isoFmerLookup, ("Pep2Gene database does not contain leucine-isoleucine " 
                               "ambiguity data; re-compile database or "
                               "select distinguish_leucine = False .")
    assert k_len == K, "Pep2Gene database created with kmers of length %s, not %s" % (k_len, K)
    
    print "P2G database loaded: %.2f\n\n" % (time.clock() - starttime)
    prevtime = time.clock()
    
    rdr = reader(target_file, sheet_name = target_sheet)
    
    
    add_legacy_cols = ["pro_count","pro_list",
                       "gene_count","gene_symbols","entrez_gene_ids"]

    add_cols = ["Protein Count", "Proteins",
                "Gene Count", "Gene Symbols", "Gene IDs"]
    if legacy_columns:
        new_cols = add_legacy_cols
        colname = dict(zip(add_cols, add_legacy_cols))
    else:
        new_cols = add_legacy_cols
        colname = dict(zip(add_cols, add_cols))
        
    iso_legacy_cols = ['IL Ambiguity pro_count', 'IL Ambiguity pro_list', 
                       "IL Ambiguity gene_count", "IL Ambiguity gene_symbols",
                       "IL Ambiguity entrez_gene_ids"]
    iso_cols = ['I<->L Protein Count', 'I<->L Proteins', 'I<->L Gene Count',
                'I<->L Gene Symbols', 'I<->L Gene IDs']
    if legacy_columns and distinguish_leucine:
        new_cols += iso_legacy_cols
        colname.update(dict(zip(iso_cols, iso_legacy_cols)))
    elif distinguish_leucine:
        new_cols += iso_cols
        colname.update(dict(zip(iso_cols, iso_cols)))
    
        
    if not outputfile:
        ext = outputfile.split('.')[-1]
        outputfile = '.'.join(outputfile.split('.')[:-1] + ['GENES', ext])
    output = writer(outputfile,
                    columns = rdr.columns + new_cols)
    
    pepToProts = {}
    isoPepToProts = {}
    for counter, row in enumerate(rdr):
        if counter%1000 == 0:
            print_progress(counter)
        try:
            pep = row['Peptide Sequence'].upper()
        except KeyError:
            pep = row['peptide sequence'].upper()
            
        isoPep = pep.replace('I', 'L')
        if pep not in pepToProts:
            candidate_prots = reduce(set.intersection,
                                     (fmerLookup[pep[x:x+K]] for x in range(len(pep)-K)))
            pep_find = re.compile('((^M?)|[KR](?=[^P]))%s(((?<=[KR])[^P])|$)' % pep)
            pepToProts[pep] = set(prot for prot in candidate_prots
                                  if pep_find.search(seqLookup[prot]))
            
            if distinguish_leucine:
                assert isoPep not in isoPepToProts
                iso_candidate_prots = reduce(set.intersection,
                                             (isoFmerLookup[isoPep[x:x+K]] for x
                                              in range(len(isoPep)-K)))
                pep_find = re.compile('((^M?)|[KR](?=[^P]))%s(((?<=[KR])[^P])|$)' % isoPep)
                isoPepToProts[isoPep] = set(prot for prot in iso_candidate_prots
                                            if pep_find.search(isoSeqLookup[prot]))
                
                
        
        proteins = '; '.join(pepToProts[pep])
        proteinCount = len(pepToProts[pep])
        
        geneList = set(geneLookup[x] for x in pepToProts[pep] if x in geneLookup)
        geneIds = '; '.join(set(g[0] for g in geneList))
        #geneSymbols = '; '.join(set(s for _, s in geneList))
        geneCount = len(geneList)
        
        row[colname['Protein Count']] = proteinCount
        row[colname['Proteins']] = proteins
        row[colname['Gene Count']] = geneCount
        row[colname['Gene Symbols']] = geneIds
        #row[colname['Gene IDs']] = 
        
        if distinguish_leucine:
            isoProteins = '; '.join(isoPepToProts[isoPep])
            isoProteinCount = len(isoPepToProts[isoPep])
            
            isoGeneList = set(geneLookup[x] for x in isoPepToProts[isoPep] if x in geneLookup)
            isoGeneIds = '; '.join(set(g[0] for g in isoGeneList))
            #isoGeneSymbols = '; '.join(set(s for _, s in isoGeneList))
            isoGeneCount = len(isoGeneList)
            
            row[colname['I<->L Protein Count']] = isoProteinCount
            row[colname['I<->L Proteins']] = isoProteins
            row[colname['I<->L Gene Count']] = isoGeneCount
            row[colname['I<->L Gene Symbols']] = isoGeneIds
            #row[colname['I<->L Gene IDs']] = 
        
        output.write(row)
    
    print "\nGene lookup completed: %.2f" % (time.clock() - prevtime)
    prevtime = time.clock()
    rdr.close()
    output.close()
    print "Output written: %.2f" % (time.clock() - prevtime)
    return outputfile
    
        
            
        
        

        
if __name__ == '__main__':
    create_fasta_index(r'C:\Users\Max\Downloads\UP000005640_9606.fasta\Human_Uniprot_Full_8-4-16.fasta', 
                     r'C:\Users\Max\Downloads\UP000005640_9606.fasta\Human_Uniprot_Full_8-4-16_newMode.pickle',
                     labelParser = lambda x: x.split('|')[1])
                     
if __name__ == '__main__':
    add_gene_ids(r'C:\Users\Max\Desktop\Projects\guillaume_gene_analysis\somethingElse\20160710_Eif4g1-1-2.xlsx.SILAC_annotated_FDR.xlsx',
                 r'C:\Users\Max\Downloads\UP000005640_9606.fasta\Human_Uniprot_Full_8-4-16_newMode.pickle',
                 outputfile = r'C:\Users\Max\Desktop\Projects\guillaume_gene_analysis\somethingElse\timetest.xlsx',
                 target_sheet = None)    
    