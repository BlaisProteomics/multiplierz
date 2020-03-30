from multiplierz.mzReport import reader, writer
from multiplierz.fasta import parse_to_dict, parse_to_generator
from multiplierz.settings import settings
from multiplierz.internalAlgorithms import print_progress

import pickle
import time
import os
import re
from collections import defaultdict
import requests
from functools import reduce

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

    url = 'https://www.uniprot.org/uploadlists/'
    
    geneNames = defaultdict(list)
    for i in range(0, len(accessions), 1000):
        subAccList = accessions[i:i+1000]
        
        geneNameReq = {
            'from':fromKey,
            'to':'GENENAME',
            'format':'tab',
            'query':' '.join(subAccList),
        }    
        
        contact = settings.user_email    
    
        res = requests.post(url, data = geneNameReq)
        geneNamePage = res.text
        
        table = [x.split('\t') for x in geneNamePage.split('\n')[1:] if x]
        for row in table:
            acc = row[0]
            # The gene field seems to be a list of synonymous names,
            # which should be considered as a unit.
            try:
                geneNames[acc].append(row[1])
            except IndexError as err:
                print(row)
                raise err
    
    normedGeneNames = dict()
    for acc, genes in list(geneNames.items()):
        normedGeneNames[acc.split('.')[0]] = ','.join(genes)
    
    return normedGeneNames#, uniprotIDs



def readLocalMapForFasta(local_map_file, accessions, accession_parser = None):
    # Local mapping file, for now, is just assumed to be a pickled dict. 
    # Only including relevant accessions purely for the sake of output file size.
    import pickle
    local_map = pickle.load(open(local_map_file, 'rb'))
    output_map = {}
    valid = 0
    accession_re = re.compile('(NP_[0-9]{6,15})|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')
    
    for accession in accessions:
        match = accession_re.search(accession)
        if match:
            acc = match.group()
            try:
                output_map[acc] = local_map[acc]
                valid += 1
            except KeyError:
                continue
        
    if valid < float(len(accessions))/10:
        raise ValueError("Only %s matching accessions out of %s found in %s" % 
                           (valid, len(accessions), local_map_file))
    print("Local map file match: %s/%s" % (valid, len(accessions)))
    return output_map
    
    


def create_fasta_index(fasta_file, outputfile = None, labelParser = (lambda x: x),
                       distinguish_leucine = True,
                       local_mapping_file = None):
    starttime = time.clock()
    
    if not outputfile:
        outputfile = fasta_file + '.pep2gene'
    
    if isinstance(labelParser, str):
        # If called from the GUI, a function can't be passed for
        # the sake of async.  But presumably the string is a regex.
        parser = re.compile(labelParser)
        def labelParser(label):
            parsed = parser.search(label)
            if not parsed:
                raise RuntimeError(("Could not parse FASTA label %s "
                                      "with provided regular expression %s") % 
                                     (label, parser.pattern))
            return parsed.group(0)  
    
    prevtime = time.clock()
    fasta = parse_to_dict(fasta_file, labelConverter = labelParser,
                          filter_string = 'rev_')
    print("FASTA file read: %.2f" % (time.clock() - prevtime))
    prevtime = time.clock()
    
    if not local_mapping_file:
        prot_to_gene = convertAccessionsViaUniprot(list(fasta.keys()))
    else:
        prot_to_gene = readLocalMapForFasta(local_mapping_file, list(fasta.keys()))
    assert any(prot_to_gene.values()), 'WARNING: Gene mapping acquisition failed to find any genes.  Check your header parser?'
    print("Gene lookup acquired: %.2f" % (time.clock() - prevtime))
    prevtime = time.clock()
    
    if distinguish_leucine:
        isofasta = {}
        for prot, seq in list(fasta.items()):
            isofasta[prot] = seq.replace('I', 'L')
    else:
        isofasta = None
        
    fmerLookup = defaultdict(set)
    isomerLookup = defaultdict(set)
    for protein, seq in list(fasta.items()):
        for i in range(len(seq) - 3):
            fmer = seq[i:i+K]
            fmerLookup[fmer].add(protein)
            if distinguish_leucine:
                isomerLookup[fmer.replace('I', 'L')].add(protein)
    
    print("Lookup tables created: %.2f" % (time.clock() - prevtime))
    prevtime = time.clock()
    
    out = open(outputfile, 'wb')
    data = (K, fasta, fmerLookup, prot_to_gene, isofasta, isomerLookup)
    pickle.dump(data, out, protocol = 2)
    out.close()
    
    print("Data written: %.2f" % (time.clock() - prevtime))
    print("Overall database generation time: %.2f" % (time.clock() - starttime))
    
    return outputfile
    
    
def add_gene_ids(target_files, p2g_database,
                 target_sheet = None,
                 outputfile = None, inPlace = False, leucine_equals_isoleucine = True,
                 legacy_columns = True):
    starttime = time.clock()
    
    if isinstance(target_files, str):
        return_list = False
        target_files = [target_files]
    else:
        return_list = True
    
    dataRdr = open(p2g_database, 'rb')
    data = pickle.load(dataRdr)
    k_len = None
    if isinstance(data, tuple) and len(data) == 6:
        k_len, seqLookup, fmerLookup, geneLookup, isoSeqLookup, isoFmerLookup = data
    elif isinstance(data, tuple) and not len(data) == 6:
        raise Exception(str(len(data)))
    else:
        print('Legacy mode P2G database detected!')
        seqLookup = data
        fmerLookup = pickle.load(dataRdr)
        geneLookup = pickle.load(dataRdr)
        try:
            isoSeqLookup = pickle.load(dataRdr)
            isoFmerLookup = pickle.load(dataRdr)
        except EOFError:
            distinguish_leucine = False
            isoSeqLookup = None
            isoFmerLookup = None
    dataRdr.close()
    
    if isinstance(list(geneLookup.values())[0], tuple):
        print("Legacy mode gene names detected.")
        oldTupleInstance = list(geneLookup.values())[0]
        nameIndex = 0 if oldTupleInstance[0] and any(x.isalpha() for x in oldTupleInstance[0]) else 1
        for k, v in list(geneLookup.items()):
            geneLookup[k] = v[nameIndex]
    
    if leucine_equals_isoleucine:
        assert isoFmerLookup, ("Pep2Gene database does not contain leucine-isoleucine " 
                               "ambiguity data; re-compile database or "
                               "select leucine_equals_isoleucine = False .")
    if k_len:
        assert k_len == K, "Pep2Gene database created with kmers of length %s, not %s" % (k_len, K)
    
    print("P2G database loaded: %.2f\n\n" % (time.clock() - starttime))
    prevtime = time.clock()
    
    outputfiles = []
    for target_file in target_files:
        try:
            rdr = reader(target_file, sheet_name = target_sheet)
        except TypeError:
            rdr = reader(target_file) # Not an Excel file.
        
        
        add_legacy_cols = ["pro_count","pro_list",
                           "gene_count","gene_symbols",]
    
        add_cols = ["Protein Count", "Proteins",
                    "Gene Count", "Gene Symbols"]
        if legacy_columns:
            new_cols = add_legacy_cols
            colname = dict(list(zip(add_cols, add_legacy_cols)))
        else:
            new_cols = add_cols
            colname = dict(list(zip(add_cols, add_cols)))
            
        iso_legacy_cols = ['IL Ambiguity pro_count', 'IL Ambiguity pro_list', 
                           "IL Ambiguity gene_count", "IL Ambiguity gene_symbols"]
        iso_cols = ['I<->L Protein Count', 'I<->L Proteins', 'I<->L Gene Count',
                    'I<->L Gene Symbols']
        if legacy_columns and leucine_equals_isoleucine:
            new_cols += iso_legacy_cols
            colname.update(dict(list(zip(iso_cols, iso_legacy_cols))))
        elif leucine_equals_isoleucine:
            new_cols += iso_cols
            colname.update(dict(list(zip(iso_cols, iso_cols))))
        
            
        if (not outputfile) or return_list:
            ext = target_file.split('.')[-1]
            outputfile = '.'.join(target_file.split('.')[:-1] + ['GENES', ext])
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
                pep = row['Peptide'].upper()
            
            pep = ''.join([x for x in pep if x.isalpha()])
                
            if len(pep) <= K:
                continue # No 4-mers in a 3-mer!
                
            isoPep = pep.replace('I', 'L')
            if pep not in pepToProts:
                candidate_prots = reduce(set.intersection,
                                         (fmerLookup[pep[x:x+K]] for x in range(len(pep)-K)))
                # pep_find could be replaced by giving the p2g database a pre-made set of
                # hashes of all tryptic peptides in a protein, and seeing if the hash of the
                # pep is present in the set.
                pep_find = re.compile('((^M?)|[KR](?=[^P]))%s(((?<=[KR])[^P])|$)' % pep)
                pepToProts[pep] = set(prot for prot in candidate_prots
                                      if pep_find.search(seqLookup[prot]))
                
                if leucine_equals_isoleucine and isoPep not in isoPepToProts:
                    iso_candidate_prots = reduce(set.intersection,
                                                 (isoFmerLookup[isoPep[x:x+K]] for x
                                                  in range(len(isoPep)-K)))
                    pep_find = re.compile('((^M?)|[KR](?=[^P]))%s(((?<=[KR])[^P])|$)' % isoPep)
                    isoPepToProts[isoPep] = set(prot for prot in iso_candidate_prots
                                                if pep_find.search(isoSeqLookup[prot]))
                    
                    
            
            proteins = '; '.join(pepToProts[pep])
            proteinCount = len(pepToProts[pep])
            
            geneList = set(geneLookup[x] for x in pepToProts[pep] if x in geneLookup)
            geneIds = '; '.join(set(g for g in geneList))
            #geneSymbols = '; '.join(set(s for _, s in geneList))
            geneCount = len(geneList)
            
            row[colname['Protein Count']] = proteinCount
            row[colname['Proteins']] = proteins
            row[colname['Gene Count']] = geneCount
            row[colname['Gene Symbols']] = geneIds
            #row[colname['Gene IDs']] = 
            
            if leucine_equals_isoleucine:
                isoProteins = '; '.join(isoPepToProts[isoPep])
                isoProteinCount = len(isoPepToProts[isoPep])
                
                isoGeneList = set(geneLookup[x] for x in isoPepToProts[isoPep] if x in geneLookup)
                isoGeneIds = '; '.join(set(g for g in isoGeneList))
                #isoGeneSymbols = '; '.join(set(s for _, s in isoGeneList))
                isoGeneCount = len(isoGeneList)
                
                row[colname['I<->L Protein Count']] = isoProteinCount
                row[colname['I<->L Proteins']] = isoProteins
                row[colname['I<->L Gene Count']] = isoGeneCount
                row[colname['I<->L Gene Symbols']] = isoGeneIds
                #row[colname['I<->L Gene IDs']] = 
            
            output.write(row)
        
        print("\nGene lookup completed: %.2f" % (time.clock() - prevtime))
        prevtime = time.clock()
        rdr.close()
        output.close()
        print("Output written: %.2f" % (time.clock() - prevtime))
        outputfiles.append(outputfile)
    if return_list:
        return outputfiles
    else:
        return outputfile
