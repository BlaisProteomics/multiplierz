import gzip
from multiplierz.mass_biochem import digest
from multiplierz.internalAlgorithms import gzOptOpen
import re

__all__ = ['Writer', 'parse_to_dict', 'parse_to_generator', 'write_fasta',
           'partial_database', 'reverse_database', 'combine', 'pseudo_reverse']






def parse_to_dict(fastaFile, labelConverter = (lambda x: x), 
                  nonredundant_labels = False,
                  filter_string = None):
    '''
    Parses a FASTA file and provides a dict object from headers to corresponding
    sequences.  Takes longer to initialize than the generator, but random-access
    lookup of sequences can be performed much faster.
    '''

    fasta = gzOptOpen(fastaFile, mode = 'r')

    index = {}
    nextKey = None
    nextSeq = ''
    for line in fasta:
        if line[0] == ">":
            if nextSeq:
                if not filter_string or filter_string not in nextKey:
                    assert nextKey, "Missing or invalid identification line?"
                    index[labelConverter(nextKey.strip())] = nextSeq.replace('\n', '')
                nextSeq = ''

            nextKey = line[1:]

        else:
            nextSeq += line
    
    if not nonredundant_labels:
        label = labelConverter(nextKey.strip())
        index[label] = nextSeq.replace('\n', '')
    else:
        seq = nextSeq.replace('\n', '')
        keys = nextKey.split('\x01')
        for key in keys:
            index[key] = seq
    

    fasta.close()
    return index

def parse_to_generator(fasta_file, labelConverter = (lambda x: x),
                       filter_string = None):
    '''
    Parses a FASTA file and provides a generator of (header,sequence) tuples.
    Can be used efficiently even on very large FASTA files.

    Example usage:
    for header,sequence in parse_fasta(file_name):
        print header
        process(sequence)
    '''

    fasta = gzOptOpen(fasta_file, mode = 'r')

    header = ''
    sequence = []

    for line in fasta:
        if line[0] == '>':
            if header and sequence:
                if (not filter_string) or filter_string not in header:
                    # It would be slightly more efficient to check for
                    # filter_string when the header is first encountered,
                    # and not aggregate the sequence in that case?
                    yield (labelConverter(header), ''.join(sequence))
            sequence = []
            header = line[1:].strip()
        else:
            sequence.append(line.strip())	

    yield (labelConverter(header), ''.join(sequence))

    fasta.close()

# For the sake of symmetry with Writer(); since generators have a .close()
# function inherently, the methods even match!
Reader = parse_to_generator

class Writer():
    """
    Easily writes properly formatted FASTA entries.
    """
    def __init__(self, filename, line_length = 80):
        self.file = gzOptOpen(filename, mode = 'w')
        self.length = line_length
    
    def write(self, header, sequence):
        if not header[0] == '>':
            header = '>' + header
        
        sequence = '\n'.join([sequence[i : i + self.length]
                              for i in range(0, len(sequence), self.length)])
        
        self.file.write(header + '\n')
        self.file.write(sequence + '\n')
    
    def close(self):
        self.file.close()
        
        
        
            


def write_fasta(fasta, save_file, write_mode='w'):
    '''
    Creates a fasta file given a dictionary with prot: sequence or an iterator of
    (prot, sequence) pairs.
    This method adds ">" at the beginning of each header if there isn't one already.
    The proteins will not be written in any particular order.
    '''

    fasta_file = gzOptOpen(save_file, mode = write_mode)

    if type(fasta) == type(dict()):
        for prot, seq in list(fasta.items()):
            if prot[0] == '>':
                fasta_file.write("%s\n" % prot.strip())
            else:
                fasta_file.write(">%s\n" % prot.strip())
            fasta_file.write("%s\n\n" % seq)
    else:
        for prot, seq in fasta:
            if prot[0] == '>':
                fasta_file.write("%s\n" % prot.strip())
            else:
                fasta_file.write(">%s\n" % prot.strip())
            fasta_file.write("%s\n\n" % seq)	

    fasta_file.close()



def partial_database(fasta, output = None, search = '',
                     use_regex = False, include_matches = True):
    """
    fasta -> A target FASTA-format file
    output -> The output file (input file is overwritten if this is not specified.)
    search -> Either a string or a list of strings.
    use_regex -> If this is set to True, strings in search will be interpreted
    as regular expressions.  Otherwise, a search string "matches" if it is contained
    by the header.
    include_matches -> If this is True (the default,) FASTA entries are included in
    the output if they match a search string; if False, FASTA entries are included
    if they match NONE of the search strings.
    
    Creates a new fasta database by copying each entry from the original, according
    to the rules outlined above.
    """

    if not output:
        output = fasta + 'partial_database.fasta'

    if isinstance(search, str):
        search = [search]
    if use_regex:
        search = list(map(re.compile, search))
        
    
    fastaGen = parse_to_generator(fasta)

    out = gzOptOpen(output, mode = 'w')

    if use_regex:
        for header, sequence in fastaGen:
            if any([x.search(header) for x in search]) == include_matches:
                out.write('>' + header + '\n')
                out.write(sequence + '\n')                
    else:
        for header, sequence in fastaGen:
            if any([x in header for x in search]) == include_matches:
                out.write('>' + header + '\n')
                out.write(sequence + '\n')

    out.close()
    
    return output


def reverse_database(fasta, outputFile = None, include_forward = False, tag = 'rev_'):
    '''
    Writes a reverse-sequence database, where each header is annotated
    with 'rev_' and each sequence is reversed from the original.

    If includeForward is False (default), only the reverse sequences are
    included.  If True, the original database is included in the new database.
    (Combined forward-reverse databases are useful in early versions of Mascot,
    which only search one database at a time.)
    '''

    if not outputFile:
        if fasta.lower().endswith('.fasta'):
            inputName = fasta[:-6]
            
        if include_forward:
            outputFile = inputName + '-forward_reverse_db.fasta'
        else:
            outputFile = inputName + '-reverse_db.fasta'

    output = Writer(outputFile)

    inputGen = Reader(fasta)

    for protein, sequence in inputGen:
        revProtein = tag + protein
        revSequence = ''.join(list(reversed(sequence)))

        if include_forward:
            output.write(protein, sequence)
        output.write(revProtein, revSequence)

    output.close()	
    
    return outputFile
    
    
    
def combine(fasta_files, output):
    """
    Simply combines a given list of FASTA files together.  Aside from gzip support,
    this is just file concatenation.
    """
    
    outputFile = gzOptOpen(output, mode = 'w')
    for fasta in fasta_files:
        fastaFile = gzOptOpen(fasta, mode = 'r')
        outputFile.write(fastaFile.read())
        fastaFile.close()
    outputFile.close()
    
    return output
    
def pseudo_reverse(fasta, output = None, enzyme = 'Trypsin', tag = 'rev_', include_forward = False):
    """
    Creates a psuedo-reverse database, in which each individual fragment sequence
    specified by the given enzyme is reversed.  The cleavage sites themselves are
    preserved and not reversed.
    """
    
    if not output:
        output = fasta + '.pseudo_reversed.fasta'
    
    fastaFile = parse_to_generator(fasta)
    outputFile = Writer(output)
    
    for header, sequence in fastaFile:
        subsequences = [x[0] for x in digest(sequence, enzyme)]
        
        reconstructed = []
        for subseq in subsequences:
            reconstructed.append(subseq[:-1][::-1] + subseq[-1])
        reconstructed = ''.join(reconstructed)
        
        if include_forward:
            outputFile.write(header, sequence)
        outputFile.write(tag+header, reconstructed)
    
    outputFile.close()

    return output

