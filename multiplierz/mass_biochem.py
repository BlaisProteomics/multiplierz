# Copyright 2008 Dana-Farber Cancer Institute
# multiplierz is distributed under the terms of the GNU Lesser General Public License
#
# This file is part of multiplierz.
#
# multiplierz is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiplierz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiplierz.  If not, see <http://www.gnu.org/licenses/>.


import os
import re

from random import randint, sample
from functools import reduce
try:
    from collections import defaultdict, Counter
except ImportError:
    from collections import defaultdict
    print("Some functions may not currently work in Python 2.6.")
    
from itertools import permutations

from multiplierz import myData, logger_message, protonMass
from multiplierz.unimod import UnimodDatabase



#__all__ = ['digest', 'fragment', 'mz', 'mw']

unimodpath = os.path.join(myData, 'unimod.sqlite')
try:
    # load the unimod database
    unimod = UnimodDatabase(unimodpath)
except IOError:
    print ("Unimod database not found at %s!\n"
           "This may be due to an incomplete multiplierz installation.  Please\n"
           "reinstall multiplierz or copy the unimod.sqlite file from an existing\n"
           "multiplierz repository.")
    

# From Unimod.
AtomicMass = {'13C': 13.00335483,
              '15N': 15.00010897,
              '18O': 17.9991603,
              '2H': 2.014101779,
              'Ag': 106.905092,
              'As': 74.9215942,
              'Au': 196.966543,
              'B': 11.0093055,
              'Br': 78.9183361,
              'C': 12.0,
              'Ca': 39.9625906,
              'Cd': 113.903357,
              'Cl': 34.96885272,
              'Co': 58.9331976,
              'Cr': 51.9405098,
              'Cu': 62.9295989,
              'F': 18.99840322,
              'Fe': 55.9349393,
              'H': 1.007825035,
              'Hg': 201.970617,
              'I': 126.904473,
              'K': 38.9637074,
              'Li': 7.016003,
              'Mg': 23.9850423,
              'Mn': 54.9380471,
              'Mo': 97.9054073,
              'N': 14.003074,
              'Na': 22.9897677,
              'Ni': 57.9353462,
              'O': 15.99491463,
              'P': 30.973762,
              'Pd': 105.903478,
              'S': 31.9720707,
              'Se': 79.9165196,
              'Zn': 63.9291448,
              'e': 0.000549}
AW = AtomicMass # For legacy reasons.

#AverageMasses = {'Br': 78.9183361,
                 #'C': 12.01104,
                 #'Cl': 34.968852721,
                 #'D': 2.014,
                 #'F': 18.9984046,
                 #'Fe': 55.845,
                 #'H': 1.007976,
                 #'N': 14.00666,
                 #'O': 15.99932,
                 #'P': 30.973762,
                 #'S': 32.06469,
                 #'e': 0.000549}
AverageMasses = {'Al': 26.98154,
                 'Ar': 39.948,
                 'As': 74.92159,
                 'B': 10.811,
                 'Be': 9.012182,
                 'Br': 79.904,
                 'C': 12.011,
                 'Ca': 40.078,
                 'Cl': 35.4527,
                 'Co': 58.9332,
                 'Cr': 51.9961,
                 'Cu': 63.546,
                 'F': 18.9984,
                 'Fe': 55.847,
                 'Ga': 69.723,
                 'Ge': 72.61,
                 'H': 1.00794,
                 'He': 4.002602,
                 'K': 39.0983,
                 'Li': 6.941,
                 'Mg': 24.305,
                 'Mn': 54.93805,
                 'N': 14.00674,
                 'Na': 22.98977,
                 'Ne': 20.1797,
                 'Ni': 58.6934,
                 'O': 15.9994,
                 'P': 30.97376,
                 'S': 32.066,
                 'Sc': 44.95591,
                 'Se': 78.96,
                 'Si': 28.0855,
                 'Ti': 47.88,
                 'V': 50.9415,
                 'Zn': 65.39}                 
Avg_AW = AverageMasses

# dictionary of modification shorthand
mod_dictionary = {'p': 'Phospho',
                  'i': 'iTRAQ4plex',
                  'i8': 'iTRAQ8plex',
                  'o': 'Oxidation',
                  'c': 'Carbamidomethyl',
                  'd': 'Deamidated',
                  's4': 'Label:13C(3)15N(1)',
                  's6': 'Label:13C(6)',
                  's10': 'Label:13C(6)15N(4)'}

mod_masses = {'Phospho' : 79.966331,
              'Carbamidomethyl' : 57.021464,
              'Deamidated' : 0.984016,
              'iTRAQ4plex' : 144.102063,
              'iTRAQ8plex': 304.205360,
              'TMT6plex':229.162932,
              'Label:13C(6)15N(4)' : 10.008269,
              'Label:13C(6)' : 6.020129,
              'Label:13C(3)15N(1)' : 4.007099,
              'Oxidation' : 15.994915,
              'Methyl' : 14.01565,
              }

ModificationFormulae = {'Phospho':{'H':1, 'O':3, 'P':1},
                        'Carbamidomethyl':{'H':3, 'C':2, 'N':1, 'O':1},
                        'Deamidated' : {'H':-1, 'N':-1, 'O':1},
                        'Oxidation' : {'O':1},
                        'Methyl' : {'C':1, 'H':2},
                        'iTRAQ4plex' : {'H':12, 'C':4, '13C':3, 'N':1, '15N':1, 'O':1},
                        'iTRAQ8plex' : {'H':24, 'C':7, '13C':7, 'N':3, '15N':1, 'O':3},
                        'TMT6plex' : {'H':20,'C':8, '13C':4, 'N':1, 'O':2, '15N':1},
                        'Label:13C(6)': {'13C':6, 'C':-6},
                        'Label:13C(6)15N(4)':{'13C':6, 'C':-6, '15N':4, 'N':-4}}


EnzymeSpecification = {'Arg-C':          '[R][A-Z]',
                       'Asp-N':          '[A-Z][D]',
                       'Bromelain':      '[KAY][A-Z]',
                       'CNBr_HSer':      '[M][A-Z]',
                       'CNBr_HSerLac':   '[M][A-Z]',
                       'Cathepsin B':    '[R][A-Z]',
                       'Cathepsin D':    '[LF][^VAG]',
                       'Cathepsin G':    '[YWF][A-Z]',
                       'Chymotrypsin':   '[YWFL][A-Z]',
                       'Clostripain':    '[R][P]',
                       'Elastase':       '[AVLIGS][A-Z]',
                       'Glu-C_Bic':      '[E][A-Z]',
                       'Glu-C_Phos':     '[ED][A-Z]',
                       'Hydroxylamine':  '[N][G]',
                       'Lys-C':          '[K][A-Z]',
                       'Lys-N':          '[A-Z][K]',
                       'Papain':         '[RK][A-Z]',
                       'Pepsin':         '[LF][^VAG]',
                       'Proteinase K':   '[YWF][A-Z]',
                       'Subtilisin':     '[^RHK][A-Z]',
                       'Thermolysin':    '[LFIVMA][^P]',
                       'Trypsin':        '[KR][^P]'
                       }



# Values are (monoisotopic mass, average mass)
# Courtesy of ExPASy: http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
AminoAcidMasses = {'A': (71.03711, 71.0788),
                   'C': (103.00919, 103.1388),
                   'D': (115.02694, 115.0886),
                   'E': (129.04259, 129.1155),
                   'F': (147.06841, 147.1766),
                   'G': (57.02146, 57.0519),
                   'H': (137.05891, 137.1411),
                   'I': (113.08406, 113.1594),
                   'K': (128.09496, 128.1741),
                   'L': (113.08406, 113.1594),
                   'M': (131.04049, 131.1926),
                   'N': (114.04293, 114.1038),
                   'P': (97.05276, 97.1167),
                   'Q': (128.05858, 128.1307),
                   'R': (156.10111, 156.1875),
                   'S': (87.03203, 87.0782),
                   'T': (101.04768, 101.1051),
                   'V': (99.06841, 99.1326),
                   'W': (186.07931, 186.2132),
                   'Y': (163.06333, 163.176)}

# These all include a bound H2O.
AminoAcidFormulas = {'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2},
                     'C': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1},
                     'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4},
                     'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4},
                     'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2},
                     'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2},
                     'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2},
                     'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
                     'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2},
                     'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
                     'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 1},
                     'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3},
                     'P': {'C': 5, 'H': 9, 'N': 1, 'O': 2},
                     'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3},
                     'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2},
                     'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3},
                     'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3},
                     'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2},
                     'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2},
                     'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3}}

H2Omass = 18.010565
NH3mass = 17.026549


formula_form = re.compile('([A-Z][0-9]*)*$')


def sum_formula(formulas):
    totals = defaultdict(int)
    for form in formulas:
        for atom, count in list(form.items()):
            totals[atom] += count
    return totals

#formula_parser = re.compile(r'([A-Z][0-9]{0,3})')
def parse_chemical_formula(formstring):
    formstring = formstring.strip()
    assert re.match(r'(([A-Z][a-z]?)([0-9]{0,3}))*$', formstring), '%s is not a valid chemical formula.' % formstring
    formdict = {}
    i = 0
    for atom, count in re.findall(r'([A-Z][a-z]?)([0-9]{0,3})', formstring):
        if count:
            count = int(count)
        else:
            count = 1
        formdict[atom] = count
    return formdict
        
    
    #while i < len(formstring):
        #x = formstring[i]
        #if x.isalpha():
            #try:
                #y = formstring[i+1]
                #if y.isdigit():
                    #formdict[x] = int(y)
                    #i += 1
                #else:
                    #formdict[x] = 1                    
            #except IndexError:
                #formdict[x] = 1
            #i += 1
        #else:
            #raise NotImplementedError, "Unable to parse chemical formula %s" % formstring    
    
    #return formdict

chemicalMassMemo = {}
def chemicalFormulaMass(formula):    
    assert isinstance(formula, str)
    
    if formula in chemicalMassMemo:
        return chemicalMassMemo[formula]        
    
    formdict = parse_chemical_formula(formula)
        
    mass = chemicalDataMass(formdict)
    chemicalMassMemo[formula] = mass
    
    return mass

def chemicalDataMass(formula):
    assert isinstance(formula, dict)
    mass = 0
    for element, count in list(formula.items()):
        mass += AW[element] * count
    return mass

def digest(protein, enzyme="Trypsin", missed_cleavages=0):
    '''Produces digested peptide sequences for a given protein sequence and an enzyme.

    Takes a protein sequence and digests it with the chosen enzyme.
    Optionally select # of missed cleavages
    Returns a list of tuples: [(peptide, (start,end), num_missed_cleavages)]

    The following enzymes are available:
        Arg-C
        Asp-N
        Bromelain
        CNBr_HSer
        CNBr_HSerLac
        Cathepsin B
        Cathepsin D
        Cathepsin G
        Chymotrypsin
        Clostripain
        Elastase
        Glu-C_Bic
        Glu-C_Phos
        Hydroxylamine
        Lys-C
        Lys-N
        Papain
        Pepsin
        Proteinase K
        Subtilisin
        Thermolysin
        Trypsin

    If the enzyme name starts with [ , the enzyme is used as the regular expression

    Example:
        >>> protein = 'MWKASAGHAVSIAQDDAGADDWETDPDFVNDVSEKEQRWGAKTVQGSGHQEHINIHKLRENVFQEHQTLKEKELETGPKASHGYGGKF'
        >>> enzyme = 'Trypsin'
        >>> digests = mzFunctions.digest(protein, enzyme)
        >>> for digest in digests:
        ... 	print digest[0],digest[1]
        ...
        MWK [1, 3]
        ASAGHAVSIAQDDAGADDWETDPDFVNDVSEK [4, 35]
        EQR [36, 38]
        WGAK [39, 42]
        TVQGSGHQEHINIHK [43, 57]
        LR [58, 59]
        ENVFQEHQTLK [60, 70]
        EK [71, 72]
        ELETGPK [73, 79]
        ASHGYGGK [80, 87]
        F [88, 88]
    '''

    protein = re.sub("\s+","",protein)
    protein = re.sub("\n+","",protein)

    #Enzyme Specificities
    # list of enzymes
    enzSpec = EnzymeSpecification

    digest_array = []

    if enzyme[0] == "[":
        reg_exp = enzyme
    else:
        reg_exp = enzSpec[enzyme]

    expr = re.compile(reg_exp)

    start = 0
    while True:
        m = expr.search(protein, start)
        if m is None:
            break
        else:
            digest_array.append((protein[start:(m.end()-1)], (start, m.end()-1), 0))
            start = 1 + m.start()

    if start != len(protein):
        digest_array.append((protein[start:], (start, len(protein)), 0))

    #Add missed cleavages
    if missed_cleavages:
        dig_all = digest_array[:]
        for i,p in enumerate(digest_array[:-1]):
            pep = p[0]
            start = p[1][0]
            for j in range(min((len(digest_array) - i - 1), missed_cleavages)):
                pep += digest_array[i+j+1][0]
                end = digest_array[i+j+1][1][1]
                dig_all.append((pep, (start,end),j+1))
        return dig_all

    return digest_array


def fragment_legacy(peptide, ions=('b', 'b++', 'y', 'y++'), labels=True):
    '''Calculates the ion masses for a given peptide, ion types, and residue modifications.

    Takes:
    - a peptide sequence (with possible modifications)
    - a list of ion types (defaults to b and y)
    - 'labels' flag: if True, function returns two lists: ion masses and corresponding labels

    If 'labels' is not set, the return is a list of lists of ion masses, in the same order as
    provided to the function. The length of the lists will vary--most ions types will not return
    the final ion (e.g. the b7 ion of 'PEPTIDE') because it is almost never observed.

    The function will calculate neutral losses

    The following ions are available:
        imm (immonium)
        a, a0 (H20 loss), a* (NH3 loss)
        b, b0 (H20 loss), b* (NH3 loss)
        c
        x
        y, y0 (H20 loss), y* (NH3 loss)
        z, z+1, z+2
        intya (internal ya), intyb (internal yb)

        Add a '++' after any of the ions to make it doubly charged

    Also, passing 'MH+', 'MH++', etc will calculate the peptide m/z values

    The following modifications are available:
        p: phospho = H2PO3 -H
        i: iTRAQ = H12C4C*3NN*O
        o: Oxidation = O
        c: Carbamidomethyl = CH2CONH2 -H
        d: Deamidated = O -HN

    Alternatively, modifications can be specified using brackets with either a mass
    or Unimod name inside.

    Example:
        >>> peptide = '[iTRAQ4plex]-PEPpTIDE'
        >>> ions = ['b','y','b++','y++','b0','y0']
        >>> frags = mass_biochem.fragment_legacy(peptide, ions)
        >>> b_ions = frags[0]
        >>> print b_ions
        [242.16210299999997, 371.20469600000001, 468.25746000000004, 649.27147000000002, 762.35553400000003, 877.38247699999999]
    '''

    calc_masses = []
    label_dict = {}

    ion_list = ions
    include = set(ion_list)

    sequence, mods, vm_masses, neutral_loss = mz_pep_decode(peptide)

    neutral_loss_list = []
    neutral_loss_C_term = 0
    neutral_loss_N_term = 0

    net_NL_list = []
    left_label_list = []
    right_label_list = []
    var_mod_string_list = []

    # polarity is + for now
    polarity = '+'
    double = polarity * 2

    # amino acid masses
    masses = dict([(k, v[0]) for k, v in list(AminoAcidMasses.items())])
    masses['-'] = 0.0
    masses['C-term'] = 17.00274
    masses['N-term'] = 1.007825

    min_internal_mass = 0.0
    max_internal_mass = 700.0

    running_sum = []

    # first calculate running sum of residues, including variable mods
    # ie calculate the approximate mass at each fragmentation point

    running_sum.append(masses[sequence[0]])
    if mods[1]:
        running_sum[0] += sum(vm_masses[mod] for mod in mods[1])
        neutral_loss_list.extend(neutral_loss[mod][0] for mod in mods[1])
    else:
        neutral_loss_list.append(0.0)

    for i,aa in enumerate(sequence[1:]):
        running_sum.append(running_sum[-1] + masses[aa])
        if mods[i+2]:
            running_sum[-1] += sum(vm_masses[mod] for mod in mods[i+2])
            for mod in mods[i+2]:
                neutral_loss_list.append(neutral_loss_list[-1] + neutral_loss[mod][0])
        else:
            neutral_loss_list.append(neutral_loss_list[-1])

    # n-term and c-term modifications:
    if mods[0]:
        masses['N-term'] += sum(vm_masses[mod] for mod in mods[0])
        neutral_loss_N_term = sum(neutral_loss[mod][0] for mod in mods[0])
    else:
        neutral_loss_N_term = 0.0

    if mods[-1]:
        masses['C-term'] += sum(vm_masses[mod] for mod in mods[-1])
        neutral_loss_C_term = sum(neutral_loss[mod][0] for mod in mods[-1])
    else:
        neutral_loss_C_term = 0.0

    calc_MR = running_sum[-1] + masses['N-term'] + masses['C-term']

    hydrogen = AW['H']
    carbon = AW['C']
    nitrogen = AW['N']
    oxygen = AW['O']
    electron = AW['e']

    CO = carbon + oxygen
    NH3 = nitrogen + 3 * hydrogen
    H2O = 2 * hydrogen + oxygen
    if (polarity == "+"):
        charge_mass = hydrogen - electron
    else:
        charge_mass = - hydrogen + electron

    # nl_label = the net neutral as a text string for a peak label
    # frag_seq = the sequence of a fragment
    # mod_string = substring of mods corresponding to the fragment
    #   used for figuring out where neutral losses need to be permuted
    #   for high-energy fragment, string excludes new terminal residue
    # net_nl = total neutral loss for a fragment

    ions = defaultdict(list)

    RKNQ_re = re.compile(r'[RKNQ]')
    STED_re = re.compile(r'[STED]')

    len_seq = len(sequence)

    # n-term fragments
    for i in range(len_seq - 1):
        frag_seq = sequence[:i+1]

        mod_string = mods[:i+2]

        net_nl = neutral_loss_list[i] + neutral_loss_N_term

        # label is blank if the loss is ~zero
        nl_label = '%s' % (int(round(-net_nl)) or '')

        ions['a'].append(running_sum[i] + masses['N-term'] - CO
                         - net_nl - hydrogen + charge_mass)

        ions['b'].append(running_sum[i] + masses['N-term']
                         - net_nl - hydrogen + charge_mass)

        ions['c'].append(ions['b'][-1] + NH3)

        for k in ('a', 'b', 'c'):
            if k in include:
                calc_masses.append(ions[k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

        # NH3-loss if fragment includes [RKNQ]
        #if RKNQ_re.search(frag_seq):
        if True:
            ions['a*'].append(ions['a'][-1] - NH3)
            ions['b*'].append(ions['b'][-1] - NH3)

            for k in ('a*', 'b*'):
                if k in include:
                    calc_masses.append(ions[k][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('%s(%d)' % (k, i+1))
                    right_label_list.append('')
                    var_mod_string_list.append(mod_string)
                    label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

                    if ('%s++' % k) in include:
                        ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                        calc_masses.append(ions['%s++' % k][-1])
                        net_NL_list.append(net_nl)
                        left_label_list.append('%s(%d)' % (k, i+1))
                        right_label_list.append(double)
                        var_mod_string_list.append(mod_string)
                        label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

        # water-loss if fragment includes [STED]
        #if STED_re.search(frag_seq):
        if True:
            ions['a0'].append(ions['a'][-1] - H2O)
            ions['b0'].append(ions['b'][-1] - H2O)

            for k in ('a0', 'b0'):
                if k in include:
                    calc_masses.append(ions[k][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('%s(%d)' % (k, i+1))
                    right_label_list.append('')
                    var_mod_string_list.append(mod_string)
                    label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

                    if ('%s++' % k) in include:
                        ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                        calc_masses.append(ions['%s++' % k][-1])
                        net_NL_list.append(net_nl)
                        left_label_list.append('%s(%d)' % (k, i+1))
                        right_label_list.append(double)
                        var_mod_string_list.append(mod_string)
                        label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

        # doubly-charged ions
        for k in ('a', 'b', 'c'):
            if ('%s++' % k) in include:
                ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                calc_masses.append(ions['%s++' % k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append(double)
                var_mod_string_list.append(mod_string)
                label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

    # c-term fragments
    for i in range(len_seq - 1):
        frag_seq = sequence[-i-1:]

        mod_string = mods[-i-2:]

        net_nl = (neutral_loss_list[len_seq - 1]
                  - neutral_loss_list[len_seq - 2 - i]
                  + neutral_loss_C_term)

        nl_label = '%s' % (int(round(-net_nl)) or '')

        ions['y'].append(running_sum[len_seq - 1]
                         - running_sum[len_seq - 2 - i]
                         + masses['C-term'] + hydrogen
                         - net_nl + charge_mass)

        ions['x'].append(ions['y'][-1] - 2 * hydrogen + CO)

        ions['z'].append(ions['y'][-1] - NH3)

        ions['z+1'].append(ions['z'][-1] + hydrogen)

        ions['z+2'].append(ions['z+1'][-1] + hydrogen)

        for k in ('y', 'x', 'z', 'z+1', 'z+2'):
            if k in include:
                calc_masses.append(ions[k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

        #if RKNQ_re.search(frag_seq):
        if True:
            ions['y*'].append(ions['y'][-1] - NH3)

            if 'y*' in include:
                calc_masses.append(ions['y*'][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('y*(%d)' % (i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions['y*'][-1]] = 'y*(%d)%s' % (i+1, nl_label)

                if 'y*++' in include:
                    ions['y*++'].append((ions['y*'][-1] + charge_mass) / 2)

                    calc_masses.append(ions['y*++'][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('y*(%d)' % (i+1))
                    right_label_list.append(double)
                    var_mod_string_list.append(mod_string)
                    label_dict[ions['y*++'][-1]] = 'y*(%d)%s%s' % (i+1, nl_label, double)

        #if STED_re.search(frag_seq):
        if True:
            ions['y0'].append(ions['y'][-1] - H2O)

            if 'y0' in include:
                calc_masses.append(ions['y0'][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('y0(%d)' % (i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions['y0'][-1]] = 'y0(%d)%s' % (i+1, nl_label)

                if 'y0++' in include:
                    ions['y0++'].append((ions['y0'][-1] + charge_mass) / 2)

                    calc_masses.append(ions['y0++'][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('y0(%d)' % (i+1))
                    right_label_list.append(double)
                    var_mod_string_list.append(mod_string)
                    label_dict[ions['y0++'][-1]] = 'y0(%d)%s%s' % (i+1, nl_label, double)

        # doubly-charged ions
        for k in ('y', 'x', 'z', 'z+1', 'z+2'):
            if ('%s++' % k) in include:
                ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                calc_masses.append(ions['%s++' % k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append(double)
                var_mod_string_list.append(mod_string)
                label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

    # immonium
    if 'imm' in include:
        for i,aa in enumerate(sequence):
            ions['imm'].append(masses[aa] - CO + charge_mass)
            net_nl = 0
            if mods[i+1]:
                #ions['imm'][-1] += vm_masses[mods[i+1]]
                ions['imm'][-1] += sum(vm_masses[mod] for mod in mods[i+1])
                #net_nl = neutral_loss[mods[i+1]][0]
                net_nl = sum(neutral_loss[mod][0] for mod in mods[i+1])
            #elif neutral_loss[aa]:
                #net_nl = neutral_loss[aa][0]
            ions['imm'][-1] -= net_nl

            nl_label = '%s' % (int(round(-net_nl)) or '')

            calc_masses.append(ions['imm'][-1])
            net_NL_list.append(net_nl)
            left_label_list.append(aa)
            right_label_list.append('')
            var_mod_string_list.append(mods[i+1])
            label_dict[ions['imm'][-1]] = '%s%s' % (aa, nl_label)
            
            if 'imm++' in include:
                    ions['imm++'].append((ions['imm'][-1] + charge_mass)/2)
                    calc_masses.append(ions['imm++'][-1])
                    
        

    # internals
    if 'intya' in include or 'intyb' in include:
        for i in range(len_seq - 3):
            for j in range(i+2, len_seq - 1):
                frag_seq = sequence[i+1:j+1]

                mod_string = mods[i+2:j+2]

                net_nl = neutral_loss_list[j] - neutral_loss_list[i]

                nl_label = '%s' % (int(round(-net_nl)) or '')

                ions['intyb'].append(running_sum[j] - running_sum[i]
                                     - net_nl + charge_mass)

                ions['intya'].append(ions['intyb'][-1] - CO)

                for k,lbl in (('intyb',''), ('intya','-CO')):
                    if k in include:
                        if min_internal_mass <= ions[k][-1] <= max_internal_mass:
                            calc_masses.append(ions[k][-1])
                            net_NL_list.append(net_nl)
                            left_label_list.append('%s%s' % (frag_seq, lbl))
                            right_label_list.append('')
                            var_mod_string_list.append(mod_string)
                            label_dict[ions[k][-1]] = '%s%s%s' % (frag_seq, lbl, nl_label)

    # going to ignore d/d', v, w/w' ions
    #if 'd' in include:
        #if sequence[0] == 'R':
            #ions['d'].append(masses['N-term'] - neutral_loss_N_term + charge_mass
                             #+ carbon * 2 + hydrogen * 4 + nitrogen)
            #calc_masses.append(ions['d'][-1])
            #net_NL_list.append(neutral_loss_N_term)
            #left_label_list.append('d(1)')
            #right_label_list.append('')
            #var_mod_string_list.append(mods[0]

    # if 'v' in include:
    # ...
    # if 'w' in include:
    # ...

    # If any variable mods have multiple neutral losses, we now need to
    # permute out additional calculated values
    for mod in range(1, len(vm_masses) + 1):
        if len(neutral_loss[mod]) > 1:
            for i in range(len(calc_masses)):
                if mod in reduce(set.union, var_mod_string_list[i], set()):
                    count = [m for ms in var_mod_string_list[i] for m in ms].count(mod)

                    for nl_mod in neutral_loss[mod][1:]:
                        delta = count * (nl_mod - neutral_loss[mod][0])
                        if delta:
                            net_NL_list.append(net_NL_list[i] + delta)
                            left_label_list.append(left_label_list[i])
                            right_label_list.append(right_label_list[i])
                            var_mod_string_list.append(var_mod_string_list[i])

                            nl_label = '%s' % (int(round(-net_NL_list[-1])) or '')

                            if right_label_list[i]:
                                # allow for charge
                                calc_masses.append(calc_masses[i] - delta / 2)
                            else:
                                calc_masses.append(calc_masses[i] - delta)

                            label_dict[calc_masses[-1]] = (left_label_list[i]
                                                           + nl_label
                                                           + right_label_list[i])


    # this code never worked and never ran...
    #for pep_nl_list in pep_neutral_loss.values():
        #for nl in pep_nl_list:
            #calc_masses.append((calc_MR - nl + getCharge() * charge_mass)
                               #/ pep.getCharge())
            #label_dict[calc_masses[-1]] = "M%.0f%s" % (-nl, polarity * pep.getCharge())

    #for pep_nl_list in req_pep_neutral_loss.values():
        #for nl in pep_nl_list:
            #calc_masses.append((calc_MR - nl + pep.getCharge() * charge_mass)
                               #/ pep.getCharge())
            #label_dict[calc_masses[-1]] = "M%.0f%s" % (-nl, polarity * pep.getCharge())

    mh_list = [mh for mh in include if mh.startswith('MH')]
    for mh in mh_list:
        charge = mh.count(polarity)
        if charge > 0:
            calc_masses.append((calc_MR + charge_mass * charge) / charge)
            label_dict[calc_masses[-1]] = mh

    if not labels:
        return [ions[i] for i in ion_list]
    else:
        return calc_masses, label_dict

    
    
def n15_label(pepFormula, modFormula, spMass):
    # We assume modifications are not themselves labelled.
    pepFormula['15N'] = pepFormula.pop('N', 0)
    return pepFormula, modFormula, spMass

SpecialModifications = {'N15-Labelled' : n15_label}    
mascotVarModPattern = re.compile(r'([A-Z][0-9]{1,3})|N-term|C-term: .*')
mascotFixModPattern = re.compile(r'.* \(.*\)')

def peptide_mass(peptide, modifications = None, use_monoisotopic = True,
                 special_aa_formulae = {}):
    massType = 0 if use_monoisotopic else 1
    
    if not modifications: # Due to that odd Python default-argument quirk.
        modifications = []
    elif isinstance(modifications, str):
        modifications = [x.strip() for x in modifications.split(';')]
    else:
        modifications = list(modifications)
        
    pepFormula = defaultdict(int)
    for aa in peptide:
        if aa in special_aa_formulae:
            for el, num in list(special_aa_formulae[aa].items()):
                pepFormula[el] += num
        else:
            for el, num in list(AminoAcidFormulas[aa].items()):
                pepFormula[el] += num
    
    aaChainLinks = len(peptide) - 1
    pepFormula['H'] -= aaChainLinks * 2
    pepFormula['O'] -= aaChainLinks
    
    specialMass = 0
    modFormula = defaultdict(int)
    modOps = []
    for mod in modifications:
        if mod in SpecialModifications: # Mod as entry in SpecialModifications
            modOps.append(SpecialModifications[mod])
        elif isinstance(mod, str) and mascotVarModPattern.match(mod):
            # Mod as Mascot variable mod substring.
            # Presumably this is being given once per modification, e.g.
            # 'S1: Phospho; S12: Phospho'.
            submod = mod.split()[1]
            assert submod in ModificationFormulae, ("Could not recognize variable "
                                                    "modfication %s" % mod)
            for atom, num in list(ModificationFormulae[submod].items()):
                modFormula[atom] += num
        elif isinstance(mod, str) and mascotFixModPattern.match(mod):
            # Mod as Mascot fixed mod substring.
            # This is taken as being given once for however many instances
            # of the site may be present (which may be none.)
            submod, sites = mod.split()
            assert submod in ModificationFormulae, ("Could not recognize variable "
                                                    "modfication %s" % mod)
            sites = sites.strip('()')
            instances = 0
            for site in sites:
                instances += peptide.count(site)
            if instances:
                for atom, num in list(ModificationFormulae[submod].items()):
                    modFormula[atom] += num * instances
        elif mod in ModificationFormulae:
            # Mod as plain mod name.
            for atom, num in list(ModificationFormulae[mod].items()):
                modFormula[atom] += num
        elif isinstance(mod, str) and formula_form.match(mod):
            # Mod as written-out chemical formula.
            moddict = parse_chemical_formula(mod)
            for atom, num in moddict:
                modFormula[atom] += num
        elif isinstance(mod, dict):
            # Mod as chemical formula dict.
            for atom, num in mod:
                modFormula[atom] += num
        elif isinstance(mod, float):
            # Mod as plain mass
            specialMass += mod
        else:
            try:
                # Mod as mass but in string format.
                specialMass += float(mod)
            except ValueError:
                raise IOError("Unrecognized modification type: %s" % str(mod))
    
    
    for op in modOps:
        pepFormula, modFormula, specialMass = op(pepFormula, modFormula, specialMass)
    
    combinedFormula = [(el, pepFormula[el] + modFormula[el])
                       for el in set(pepFormula) | set(modFormula)]
    assert all([x[1] >= 0 for x in combinedFormula]), "Impossible molecule; %s" % combinedFormula
    
    mass = specialMass
    if use_monoisotopic:        
        for el, num in combinedFormula:
            mass += AW[el] * num
    else:       
        for el, num in combinedFormula:
            mass += Avg_AW[el] * num  
                
    return mass
    
    
mw = peptide_mass # Legacy.

def peptide_mz(peptide, mods, charge):
    mass = peptide_mass(peptide, mods)
    return ((mass + (protonMass * charge)) / charge)   

mz = peptide_mz # Legacy.





def _placeCharge(mass, chg): # Assume uncharged to start.
    return (mass + (chg*protonMass)) / chg

knownNeutralLosses = {'Phospho':chemicalFormulaMass('H3PO4')}

def fragment(peptide, mods = [], charges = [1], 
             ions = ['b', 'y'], 
             neutralPhosLoss = False,
             neutralLossDynamics = None,
             waterLoss = False,
             ammoniaLoss = False,
             use_monoisotopic = True,
             special_AAs = {}):
    # Currently only records one neutral loss even if there's more than
    # one loss-inducing mod.  Simplicity is a virtue?

   

    massType = 0 if use_monoisotopic else 1
    aminoMasses = dict((aa,AminoAcidMasses[aa][massType]) for aa in list(AminoAcidMasses.keys()))
    aminoMasses.update(special_AAs)
    
    if neutralLossDynamics is None:
        neutralLossDynamics = {}
    assert not (neutralPhosLoss and neutralLossDynamics)
    if neutralPhosLoss:
        neutralLossDynamics[mod_masses['Phospho']] = 97.9769
    elif neutralLossDynamics:
        for key, value in list(neutralLossDynamics.items()):
            if isinstance(key, str):
                neutralLossDynamics[mod_masses[key]] = value
    
    #assert not neutralPhosLoss, "Not currently set up for phos loss!"
    #czLoss = 15.010899035
    cLoss = 17.026000420104097 # NH3+
    zLoss = 16.018724 # NH2
    
    nterminusMass = 0
    cterminusMass = H2Omass
    
    modBySite = defaultdict(list)
    if isinstance(mods, str):
        mods = mods.split('; ')
    for modstr in [x for x in mods if x]:
        if ':' in modstr:
            mod = modstr.split(': ')[1].strip()
            try:
                modmass = float(mod)
            except ValueError:
                modmass = mod_masses[mod]
            if modstr[:6].lower() == 'n-term':
                nterminusMass = modmass 
            elif modstr[:6].lower() == 'c-term':
                cterminusMass = modmass + protonMass
            else:
                site = int(modstr.split(': ')[0][1:])
                modBySite[site].append(modmass)
        else:
            # Fixed mod processing.
            assert re.match('[A-Za-z0-9\-\. ]* \([A-Z]\)', modstr), modstr
            mod, loc = modstr.split()
            try:
                modmass = float(mod)
            except ValueError:
                modmass = mod_masses[mod]
            loc = loc.strip('()')
            for site, aa in enumerate(peptide, start = 1):
                if aa == loc:
                    modBySite[site].append(modmass)
        
    
    fragmentSetsByIonType = {}
    if 'b' in ions or 'c' in ions:
        bfrags = []
        cfrags = []
        bmass = nterminusMass
        presentmods = []
        for site in range(1, len(peptide)): # 1-indexed.
            aa = peptide[site-1]
            mods = modBySite[site]
            presentmods += mods
            bmass += AminoAcidMasses[aa][massType] + sum(mod_masses.get(x, x) for x in mods)
                    
            neutralLoss = 0
            for mod in presentmods:
                if mod in neutralLossDynamics:
                    neutralLoss -= neutralLossDynamics[mod]            

            if 'b' in ions:
                bfrags.append(('b'+str(site), bmass))
            if 'c' in ions:
                cfrags.append(('c'+str(site), bmass + cLoss))

            if neutralLoss:
                if 'b' in ions:
                    bfrags.append(('b%d-%.0f'%(site, abs(neutralLoss)), bmass + neutralLoss))                
                if 'c' in ions:
                    cfrags.append(('c%d-%.0f'%(site, abs(neutralLoss)), bmass + neutralLoss + cLoss))
        if 'b' in ions:
            fragmentSetsByIonType['b'] = bfrags
        if 'c' in ions:
            fragmentSetsByIonType['c'] = cfrags
    
    if 'y' in ions or 'z' in ions:
        yfrags = []
        zfrags = []
        ymass = cterminusMass
        presentmods = []
        for site in range(len(peptide)-1, 0, -1): # 0-indexed.
            aa = peptide[site]
            mods = modBySite[site+1]
            presentmods += mods
            ymass += AminoAcidMasses[aa][massType] + sum(mod_masses.get(x, x) for x in mods)
            

            neutralLoss = 0
            for mod in presentmods:
                if mod in neutralLossDynamics:
                    neutralLoss -= neutralLossDynamics[mod]
            
            realsite = len(peptide)-site
            if 'y' in ions:
                yfrags.append(('y'+str(realsite), ymass))
            if 'z' in ions:
                zfrags.append(('z'+str(realsite), ymass - zLoss))

            if neutralLoss:
                yfrags.append(('y%d-%.0f'%(realsite, abs(neutralLoss)), ymass + neutralLoss))
                if 'z' in ions:
                    zfrags.append(('z%d-%.0f'%(realsite, abs(neutralLoss)), ymass + neutralLoss - zLoss))
        
        if 'y' in ions:
            fragmentSetsByIonType['y'] = yfrags
        if 'z' in ions:
            fragmentSetsByIonType['z'] = zfrags
                
    # Similar for z and w and whatever, I guess.
    
    if 'by-int' in ions:
        int_frags = []
        for bsite in range(1, len(peptide)):
            sub_intfrags = []
            bMass = b
    
    
    chargedFragmentSets = {}
    # Replicate sequences for multiply-charged states.
    for chg in charges:
        chg = int(chg)
        assert chg >= 1, "Positive charge states only!"
                
        for iontype in ions:
            if chg == 1:
                chgion = iontype
            else:
                chgion = iontype + '+'*chg
            chargedFragmentSets[chgion] = []
            for site, (prelabel, preion) in enumerate(fragmentSetsByIonType[iontype], start = 1):
                newion = _placeCharge(preion, chg)
                if chg == 1:
                    newlabel = prelabel
                else:
                    newlabel = prelabel + '+'*chg
                chargedFragmentSets[chgion].append((newlabel, newion))
                
                
    
    # Water-loss duplicates of ALL ions!
    if waterLoss or ammoniaLoss:
        for fragtype, labelions in list(chargedFragmentSets.items()):
            waterlosses = []
            ammonialosses = []
            for label, ion in labelions:
                chg = label.count('+') if '+' in label else 1
                mass = (ion * chg) - (chg * protonMass)                
                if waterLoss:
                    newlabel = list(label)
                    newlabel.insert(1, '0')
                    newlabel = ''.join(newlabel)
                    newmass = mass - H2Omass
                    newion = (newmass + (chg * protonMass)) / chg
                    waterlosses.append((newlabel, newion))
                if ammoniaLoss:
                    newlabel = list(label)
                    newlabel.insert(1, '*')
                    newlabel = ''.join(newlabel)
                    newmass = mass - NH3mass
                    newion = (newmass + (chg * protonMass)) / chg
                    ammonialosses.append((newlabel, newion))
            chargedFragmentSets[fragtype] += waterlosses
            chargedFragmentSets[fragtype] += ammonialosses
    
    #if 1 not in charges:
        #for iontype in ions:
            #del fragmentSetsByIonType[iontype]
    
    return chargedFragmentSets



def proline_fragments(peptide, *etc, **etcetc):
    if 'P' not in peptide:
        return []
    
    pIndices = [x[0] for x in enumerate(peptide) if x[1] == 'P']
    
    fragments = fragment(peptide, *etc, **etcetc)
    pMZs = []
    for ion, ionmz in fragments['b']:
        site = int(ion[1:].split('-')[0])
        if site in pIndices or site-1 in pIndices:
            pMZs.append(ionmz)
    for ion, ionmz in fragments['y']:
        site = len(peptide) - int(ion[1:].split('-')[0])
        if site in pIndices or site-1 in pIndices:
            pMZs.append(ionmz)  
            
    return pMZs
    
    

       


def generate_labels(scan, peptide, ions, charge=None, tolerance=0.6, **settings):
    '''Takes an MS2 scan, and a peptide (in modification format as described
    in the fragment function above) and generates a list of mass-label pairs
    by matching theoretical ion masses to experimental values within a
    tolerance.

    - scan should be a list (or other sequence) of mz,intensity pairs to be matched.
      This should probably NOT be an entire scan--that will result in a lot of
      likely-erroneous matches (to noise). One option is to use the top 50 most
      intense ions
    - peptide should be in fragment-ready format (with modification masses)
    - ions should be a list or set of ions to look for
    - tolerance is in Daltons, 0.6 is a common value
    - settings is a dictionary of optional arguments (default in brackets):
      - show_theor_mz = [ True ] | False - Display theoretical ion mass
      - ms2_mz_figs = integer >= 0, default 2 - Significant figures for mass error
      - show_mass_error = True | [ False ] - Display mass error of ions
      - mass_error_figs = integer >= 0, default 2 - Significant figures for mass error
      - mass_error_units = ['ppm'] | 'Da' - Units to display mass error

    Returns a tuple of mass-label pairs--the experimental peaks that matched
    the theoretical ion masses within tolerance.
    '''

    # initialize defaults and update with optional arguments
    _settings = dict(show_theor_mz=True, ms2_mz_figs=2, show_mass_error=False,
                     mass_error_figs=2, mass_error_units='ppm')
    _settings.update(settings)

    if _settings['mass_error_units'] == 'ppm':
        calc_error = lambda exp_mz, theor_mz: (abs(theor_mz - exp_mz) / theor_mz) * 1E6
    else:
        calc_error = lambda exp_mz, theor_mz: abs(theor_mz - exp_mz)

    if _settings['show_theor_mz'] and _settings['show_mass_error']:
        label_text = '%%s [%%.%df - %%.%df %s]' % (_settings['ms2_mz_figs'],
                                                   _settings['mass_error_figs'],
                                                   _settings['mass_error_units'])
    elif _settings['show_theor_mz']:
        label_text = '%%s [%%.%df]' % _settings['ms2_mz_figs']
    elif _settings['show_mass_error']:
        label_text = '%%s [%%.%df %s]' % (_settings['mass_error_figs'],
                                          _settings['mass_error_units'])

    if charge and 'MH' in ions:
        ions = list(ions)
        ions.remove('MH')

        for z in range(int(charge), 0, -1):
            ions.append('MH%s' % ('+'*z))

    calc_masses, label_dict = fragment_legacy(peptide, ions, labels=True)

    scan = sorted(scan)

    match_count = defaultdict(lambda: -1)

    matched_calc = []
    matched_exp = []
    matched_int = []

    # for each mass/int peak in the scan
    for j,(mass,inte) in enumerate(scan):
        # go through each calculated mass and try to match
        for i,cmass in enumerate(calc_masses):
            if abs(mass - cmass) <= tolerance:
                if match_count[i] > -1:
                    if matched_int[match_count[i]] < inte:
                        matched_exp[match_count[i]] = mass
                        matched_int[match_count[i]] = inte
                    continue

                matched_calc.append(cmass)
                matched_exp.append(mass)
                matched_int.append(inte)
                match_count[i] = len(matched_int) - 1

    matched_exp.sort()
    matched_calc.sort()

    if _settings['show_theor_mz'] and _settings['show_mass_error']:
        return tuple((e, (label_text % (label_dict[c], c, calc_error(e,c))))
                     for e,c in zip(matched_exp,matched_calc))
    elif _settings['show_theor_mz']:
        return tuple((e, (label_text % (label_dict[c], c)))
                     for e,c in zip(matched_exp,matched_calc))
    elif _settings['show_mass_error']:
        return tuple((e, (label_text % (label_dict[c], calc_error(e,c))))
                     for e,c in zip(matched_exp,matched_calc))
    else:
        return tuple((e, label_dict[c])
                     for e,c in zip(matched_exp,matched_calc))


def legacy_mz(peptide, charge=1):
    """Returns the mz for an amino acid sequence.

    Takes a peptide sequence (including possible mods) and returns the mz

    Example:
        >>> peptide = 'PEPTIDE'
        >>> totalMass = mzFunctions.mz(peptide, charge=1)
        >>> print totalMass
        799.359964319

    """

    sequence, mods, vm_masses, neutral_loss = mz_pep_decode(peptide)

    masses = dict([(k, v[0]) for k, v in AminoAcidMasses])
    masses['-'] = 0.0
    masses['C-term'] = 17.00274
    masses['N-term'] = 1.007825
    
    # randomly choose AAs for 'Z' and 'B'
    # Is this a good idea? Makes the mz function nondeterministic.
    masses['Z'] = masses[['E','Q'][randint(0,1)]]
    masses['B'] = masses[['N','D'][randint(0,1)]]
    # choose Leu for 'X'
    masses['X'] = masses['L']
    # choose Cys for Selenocysteine
    masses['U'] = masses['C']    

    total_mass = sum(masses[AA] for AA in sequence) + masses['N-term'] + masses['C-term']
    total_mass += sum(vm_masses[m] for mod in mods for m in mod)

    if charge:
        total_mass += (AW["H"] - AW["e"]) * charge
        total_mass /= charge

    return total_mass


def legacy_mw(peptide):
    """Returns the monoisotopic mass for an amino acid sequence.

    Calls mz(), with charge = 0

    """

    return mz(peptide, charge=0)



#mascotVarModPattern = re.compile(r'[A-Z][0-9]{,3}: .*')
#mascotFixModPattern = re.compile(r'.* \(.*\)')

#def mw(peptide, mods = []):        
    #pepmass = sum([AminoAcidMasses[x][0] for x in peptide])
    
    
    #if isinstance(mods, basestring):
        #mods = mods.split('; ')
        
    #modmass = 0
    #for mod in mods:
        #if not mod:
            #continue
        #elif isinstance(mod, float):
            ## Assumed to be present.
            #modmass += mod
        #elif mascotVarModPattern.match(mod):
            ## Given var mods are assumed to be present.
            #modstr = mod.split(' ')[1]
            ##modmass += mod_masses.get(modstr, float(modstr))
            #try:
                #modmass += mod_masses[modstr]
            #except KeyError:
                #modmass += float(modstr)
        #elif mascotFixModPattern.match(mod):
            #site = mod.split('(')[1].split(')')[0]
            #modstr = mod.split(' ')[0]
            #if site in peptide:
                #try:
                    #modmass += mod_masses[modstr] * peptide.count(site)
                #except KeyError:
                    #modmass += float(modstr) * peptide.count(site)
            #elif site in ['N-term', 'C-term']:
                #try:
                    #modmass += mod_masses[modstr]
                #except KeyError:
                    #modmass += float(modstr)
        #else:
            #raise ValueError, "Modification string '%s' is not in a valid format." % mod
        
    ##endMass = 0
    ###if not any(['N-term' in x for x in mods]):
    ##endMass += ElementalMasses['H']
    ###if not any(['C-term' in x for x in mods]):
    ##endMass += ElementalMasses['O'] + ElementalMasses['H']
    ## There's a ~0.004 Da difference in bound H + HO 
    ## to the raw elemental masses, and it matters.  Chemistry!
    
    #return pepmass + modmass + H2Omass





def randAA(length):
    '''
    Returns a random amino acid sequence of given length. Note:
    all amino acids will appear with equal probability.

    Example:
        >>> randomSequence = mzFunctions.randAA(10)
        >>> print randomSequence
        QGSQHQYRGD

    '''

    AA_list = ['A','C','D','E','F',
               'G','H','I','K','L',
               'M','N','P','Q','R',
               'S','T','V','W','Y']

    return ''.join(sample(AA_list, length))


def mz_pep_format(peptide, modifications='', iTRAQ=False, carb=False):
    '''Converts a given sequence and a modifications string to multiplierz format

    The modifications string is in the following format:
    "Ri: Modification; " where R = residue, i = position, Modification = modification name
    Ex: "M3: Oxidation; "

    A boolean option for iTRAQ is provided, where True indicates that all Lysines (K)
    and the N-Terminus is labeled with an 'i' modification.

    A boolean option for carbamidomethyl is provided, where True indicates that all Cysteines (C)
    are labeled with a 'c' modification.

    >>> mz_pep_format('ATLPRTLK', 'T2: Phospho; ', iTRAQ = True)
    'i-ApTLPRTLiK'

    '''

    sorted_var_mods_list = []
    n_term_mod = None
    c_term_mod = None

    if peptide.startswith('i-'):
        peptide = peptide[2:]
    peptide = re.sub('iK', 'K', peptide)
    peptide = re.sub('cC', 'C', peptide)

    if not modifications:
        modifications = ''

    vm_re = re.compile(r'\s*([NC]-term|\w(\d+))\s*:\s*(.+)')
    for item in modifications.split(';'):
        m = vm_re.match(item)
        #m = re.search('[A-Z](\d+)\s*\:\s*(.+)', item)
        if m:
            if m.group(2):
                pos = int(m.group(2)) - 1
                #mod_type = m.group(3)[0].lower()
                #if mod_type in mod_dictionary:
                    #sorted_var_mods_list.append((pos, mod_type))
                #else:
                sorted_var_mods_list.append((pos, m.group(3))) # Pure-UNIMOD style is also acceptable.
            elif m.group(1).startswith('N'):
                mod_type = m.group(3)[0].lower()
                if mod_type in mod_dictionary:
                    n_term_mod = mod_type
            elif m.group(1).startswith('C'):
                mod_type = m.group(3)[0].lower()
                if mod_type in mod_dictionary:
                    c_term_mod = mod_type

    #sorted_var_mods_list.sort()
    #for (index,(pos, mod_type)) in enumerate(sorted_var_mods_list):
        #peptide = peptide[:pos + index] + "[" + mod_type + "]" + peptide[pos + index:]
        # Adding the index is clever, but it won't work for multi-character tags!
    sorted_var_mods_list.sort(reverse = True)
    for pos, mod_type in sorted_var_mods_list:
        peptide = peptide[:pos] + "[" + mod_type + "]" + peptide[pos:]

    if n_term_mod:
        peptide = '%s-%s' % (n_term_mod, peptide)
    if c_term_mod:
        peptide = '%s-%s' % (peptide, c_term_mod)

    if iTRAQ:
        peptide = 'i-' + 'iK'.join(peptide.split('K'))

    if carb:
        peptide = 'cC'.join(peptide.split('C'))

    return peptide


def mz_pep_decode(peptide):
    '''Close to the opposite of mz_pep_format, this takes a multiplierz-format
    peptide and extracts out the basic sequence and information about modifications.

    The primary use of this function is for the mz and fragment functions above. It
    is *not* the inverse of mz_pep_format, as it does not return modification names.
    This is because it is possible for the name to be unknown, e.g. if there was a
    bare delta value in the peptide.

    The output is a 4-tuple:
    - The peptide sequence, with all mods removed
    - An array of mods, one for each amino acid as well as the N and C terminii
      (the first and last elements, respectively). Each element of the mods array
      is a set of integers, each integer being a unique mod of this peptide.
    - A dictionary of mod deltas, with keys being the integers used above.
    - A dictionary of neutral losses, keys as above and values being tuples
      of neutral loss values, (0.0,) being the minimum possibility.

    '''

    # empty mod brackets are not allowed
    peptide = re.sub(r'\[\]', '', peptide)

    pep_cleaner = re.compile(r'-?(\[.+?\]|[a-z]|\d)-?')

    # peptide without any modifications
    sequence = pep_cleaner.sub('', peptide)

    # variable modifications we can handle at the moment
    vm_masses = {}
    neutral_loss = {}

    mods = [set() for i in range(len(sequence) + 2)]

    # first regex finds mods per residue (either a bracketed group or an abbreviation)
    mod_finder_A = re.compile(r'\[(.+?)\]|([cdiop]|s\d{1,2})')
    # second regex parses a bracketed group (a '; '-delimited list of names/masses)
    mod_finder_B = re.compile(r'([-+]?\d+\.?(?:\d+)?|\.\d+)|(.+)')

    mod_set = dict()

    for mA in mod_finder_A.finditer(peptide):
        if mA.end() < len(peptide):
            AA = peptide[mA.end()]
            if AA == '-':
                AA = 'N-term'
        else:
            AA = 'C-term'

        if mA.group(1): # bracketed group
            m = []
            for g in mA.group(1).split('; '):
                mB = mod_finder_B.match(g) # this regex can't fail
                if mB.group(1): # number, no neutral losses known
                    m.append((mB.group(1),
                              float(mB.group(1)),
                              (0.0,)))
                else: # mod name (this will raise an exception if invalid)
                    m.append((mB.group(2),
                              unimod.get_mod_delta(mB.group(2)),
                              (0.0,) + (unimod.get_mod_neutral_loss(mB.group(2), AA),) ))
        else: # abbreviation
            m = [(mA.group(2),
                  unimod.get_mod_delta(mod_dictionary[mA.group(2)]),
                  (0.0,) + (unimod.get_mod_neutral_loss(mod_dictionary[mA.group(2)], AA),) )]

        # note: the first neutral loss is always 0.0, so that the mod delta is reflected in
        # the default set of ions. additional neutral losses (with corresponding labels) are
        # available if you set labels=True
        for mod,mod_delta,mod_neutral_losses in m:
            if mod not in mod_set:
                vm_masses[len(vm_masses) + 1] = mod_delta
                neutral_loss[len(vm_masses)] = mod_neutral_losses
                mod_set[mod] = len(vm_masses)

            if mA.end() < len(peptide):
                mods[len(pep_cleaner.sub('', peptide[:mA.end() + 1]))].add(mod_set[mod])
            else:
                mods[len(sequence) + 1].add(mod_set[mod])

    mods = [frozenset(m) for m in mods]

    return sequence, mods, vm_masses, neutral_loss



def mz_pep_decode_new(peptide):
    brackets = [(x.group(), x.start(), x.end()) for x in
                re.finditer('(\[\w+\])|([a-z]+)', peptide)]
    
    pure_peptide = ''
    prevstop = 0
    for _, start, stop in brackets:
        pure_peptide += peptide[prevstop:start]
        prevstop = stop
    pure_peptide += peptide[prevstop:]
    
    brackets = [(x[0], x[2]) for x in brackets]
    
    # Adjust location indexes to account for shifts from preceeding (and own)
    # mod strings, so that they'll correspond to actual points on the peptide.
    for i in range(len(brackets)-1, -1, -1):
        modstrsumlen = sum([len(x[0]) for x in brackets[:i+1]])
        brackets[i] = brackets[i][0], brackets[i][1] - modstrsumlen
    
    modsets = {}    
    modmasses = {}
    for modstrs, loc in brackets:
        assert pure_peptide[loc] in AminoAcidMasses, "Invalid peptide string: %s" % peptide
        modstrs = modstrs.strip('[]')
        
        modset = set()
        for modstr in modstrs.split(','):
            modmass = 0
            if modstr.islower() and all([x in mod_dictionary for x in modstr]):
                mod_name = mod_dictionary[x]
                modset.add(mod_name)
                modmass += sum([mod_masses[mod_name] for x in modstr])
            elif modstr in mod_masses:
                modset.add(modstr)
                modmass += mod_masses[modstr]
            elif re.match(r'[0-9]+\.?[0-9]*$', modstr):
                modset.add('Specified Mass %.2f' % modstr)
                modmass += float(modstr)
            else:
                raise NotImplementedError("Unimod lookup is disabled in this version.")
        modsets[loc+1] = modset
        modmasses[loc+1] = modmass
    for i in range(1, len(peptide)+1):
        if i not in modsets:
            modsets[i] = None
            modmasses[i] = None
    
    return pure_peptide, modsets, modmasses
        
        


def formulaForMass(mass, tolerance, components = ('C', 'H', 'N', 'O', 'P', 'S')):
    assert all(c in AW for c in components), "Not all component masses available."
    
    return list(map(dict, iter_formulaForMass(mass, tolerance, list(components))))

def iter_formulaForMass(mass, tolerance, components):
    if mass < tolerance:
        return [defaultdict(int)]
    elif not components:
        return []
    
    results = []
    atom, other_comp = components[0], components[1:]
    if mass + tolerance >= AW[atom]:
        submass = mass - AW[atom]
        subresults = iter_formulaForMass(submass, tolerance, components)
        for subr in subresults:
            subr[atom] += 1
        results += subresults
    
    results += iter_formulaForMass(mass, tolerance, other_comp)
    
    return results
    
    
    

def peptideForMass(mass, length, tolerance, pieces = None,
                   add_H2O = False, unique_sets = True): 
    """
    Takes the observed mass, the length of peptide to be searched for, and
    the mass tolerance, and returns a list of peptides within the tolerance
    of the given mass.
    
    'pieces' may be set to a list of token-mass pairs, which will then be used
    as the available peptide elements.  As a default this is set to the 20
    common amino acids and their monoisotopic masses, but it can be changed to,
    e.g., include modifications or exclude certain amino acids.
    
    If add_H2O is set to True, the mass of H2O is added to the bare
    mass of the amino acids.
    
    If unique_sets is set to True, only one example of each set of amino acids will
    be reported (so, instead of reporting "DVV,", "VDV", and "VVD", only "DVV" will
    be reported.)
    """    
    
    if not pieces:
        pieces = [(k, m) for k, (m, _) in list(AminoAcidMasses.items())]
    
    if add_H2O:
        mass -= H2Omass
    
    results = iter_peptideForMass(pieces, mass, length, tolerance)
    
    if add_H2O:
        results = [(x, m + H2Omass) for x, m in results]
        
    strResults = []
    for count, mass in results:
        resstr = ''.join(sum([[x * n] for x, n in list(count.items())], []))
        if unique_sets:
            strResults.append((resstr, mass))
        else:
            strResults += [(''.join(x), mass) for x in set(permutations(resstr))]
            
    
    return sorted(strResults)


            
def iter_peptideForMass(pieces, mass, length, tol):
    """
    Iteration function used by peptideForMass.  Don't use this directly!
    """    
    if not (pieces and length):
        if abs(mass) < tol:
            return [[Counter(), 0]]
        else:
            return []
     
    results = []
    (amino, aminomass), subpieces = pieces[0], pieces[1:]
    if mass + tol >= aminomass:
        submass = mass - aminomass
        subres = iter_peptideForMass(pieces, submass, length - 1, tol)
        for i in range(len(subres)):
            subres[i][0].update([amino]) # List-wrapped to preserve multiletter tokens.
            subres[i][1] += aminomass
        results += subres
    
    results += iter_peptideForMass(subpieces, mass, length, tol)
    
    return results
            
            
def findIonsInData(datafile, targetIons, tolerance = 0.01, 
                   threshold = 10000, thresholdMode = 'absolute',
                   includeColumns = False):
    assert thresholdMode in ['absolute', 'signalToNoise']
    from multiplierz.mzAPI import mzFile
    
    def filterscan(scan):
        foundIons = []
        for targetIon in targetIons:
            try:
                peak = max([x for x in scan if abs(x[0] - targetIon) < tolerance],
                           key = lambda x: x[1])
            except ValueError:
                foundIons.append(None)
                continue
            
            
            if thresholdMode == 'absolute' and peak[1] > threshold:
                foundIons.append((targetIon, peak))
            elif thresholdMode == 'signalToNoise':
                assert len(peak) > 2, 'Non-RAW S/N ion-thresholds not implemented here yet.'
                if peak[1] / peak[2] > threshold:
                    foundIons.append((targetIon, peak))
                else:
                    foundIons.append(None)
            else:
                foundIons.append(None)
        return foundIons
                    
                
    
    data = mzFile(datafile)
    output = []
    scaninfo = data.scan_info(0, 999999)
    for _, _, scanNum, _, _ in [x for x in scaninfo if x[3] == 'MS2']:
        if datafile.lower().endswith('raw'):
            scan = data.lscan(scanNum)
        else:
            scan = data.scan(scanNum)
        ions = filterscan(scan)
        if not all(ions):
            continue
        
        row = {'Scan Number':scanNum,
               'Total Int':sum([x[1][1] for x in ions])}
        row.update([("%s Int" % x[0], x[1][1]) for x in ions])
        output.append(row)
    
    if not includeColumns:
        return output
    else:
        return output, ['Scan Number'] + ["%s Int" % x for x in targetIons] + ['Total Int']
    




def add_protons(mass, charge):
    return (mass + (protonMass * charge)) / float(charge)

def remove_protons(mz, charge):
    return (mz * charge) - (charge * protonMass)


def isotopic_peaks_mass(mass, charge, number_peaks):
    return [add_protons(mass + i*protonMass, charge) 
            for i in range(0, number_peaks+1)]    

def isotopic_peaks_mz(mz, charge, number_peaks):
    mass = remove_protons(mz, charge)
    return isotopic_peaks_mass(mass, charge, number_peaks)
    
    
