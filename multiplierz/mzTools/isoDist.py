from collections import defaultdict
import re

__all__ = ['isotopicDistribution', 'chartDistribution']


iso = {"H":[99.99, 0.015],
       "C":[98.90, 1.10],
       "N":[99.63, 0.37],
       "O":[99.76, 0.038, 0.2],
       "P":[100.0],
       "S":[95.02, 0.75, 4.21],
       "F":[100.0],
       "Cl":[75.78, 24.22]}

def isotopicDistribution(recipe, precision = 0.000001):
    """
    Takes a dict of (atomic symbol):(atomic count) entries and
    returns the isotopic distribution of the corresponding molecule.
    
    precision <- Minimum proportional (out of 100) magnitude of isotopic
    species that the algorithm will consider.
    
    An example:
    > isotopes = isotopicDistribution({"C":10, "H":10, "N":10, "O":10, "P":10, "S":10})
    > print isotopes
    [100.0, 23.26009098259562, 48.91644916655518, 10.608690676209262, 10.877544201023955...
    chartDistribution(isotopes)
    
    """
    
    for el in recipe.keys():
        assert el in iso, "No isotopic distribution for %s" % el
    
    pattern = [100.0]
    for atom, count in recipe.items():
        for i in range(0, count):
            newPattern = defaultdict(float)
            
            for j in range(0, len(pattern)):
                for isoNumber, isoAbundance in enumerate(iso[atom], 1):
                    weightSlot = float(j + isoNumber)
                    newPattern[weightSlot] += pattern[j] * isoAbundance
            
            greatest = max(newPattern.values())
            pattern = []
            for slot in sorted(newPattern.keys()):
                normValue = (newPattern[slot]/greatest) * 100
            
                if normValue > precision:
                    pattern.append(normValue)
            
    return pattern


aminoForms = {
    "A" : [3, 7, 1, 2],
    "R" : [6, 14, 4, 2],
    "N" : [4, 8 , 2, 3],
    "D" : [4, 7, 1, 3],
    "C" : [3, 7, 1, 2, 1],
    "Q" : [5, 10, 2, 3],
    "E" : [5, 9, 1, 4],
    "G" : [2, 5, 1, 2],
    "H" : [6, 9, 1, 4],
    "I" : [6, 13, 1, 1],
    "L" : [6, 13, 1, 2],
    "K" : [6, 14, 2, 2],
    "M" : [5, 11, 1, 2, 1],
    "F" : [9, 11, 1, 2],
    "P" : [5, 9, 1, 2],
    "S" : [3, 7, 1, 3],
    "T" : [4, 9, 1, 3],
    "W" : [11, 12, 2, 2],
    "Y" : [9, 11, 1, 3],
    "V" : [5, 11, 1, 2]
    }

def renderFormula(inputString):
    numbers = re.compile(".*[0-9].*")
    segments = re.compile("([A-Z]+[a-z+]*)[ ]*\(*([0-9]+)\)*")

    recipe = {"H" : 0, "C" : 0, "N" : 0, "O" : 0, "P" : 0, "S" : 0}

    if not numbers.match(inputString): # Then this is a peptide of some sort.
        for letter in inputString:
            aminoRec = aminoForms[letter]
            recipe["C"] += aminoRec[0]
            recipe["H"] += aminoRec[1]
            recipe["N"] += aminoRec[2]
            recipe["O"] += aminoRec[3]
            try:
                recipe["S"] += aminoRec[4]
            except IndexError:
                pass
    else: # This is a chemical formula.
        terms = [x for x in segments.split(inputString) if x != '']
        for i in range(0, len(terms)-1, 2):
            atom = terms[i]
            try:
                count = int(terms[i+1])
            except ValueError:
                raise Exception("Could not parse formula.")
            assert len(atom) == 1
            assert atom in recipe.keys(), "%s currently not supported in this function." % atom

            recipe[atom] += count

    return recipe



def chartDistribution(distribution, outputfile = None):
    """
    Convenience function for displaying a matplotlib chart of an isotopic distribution
    generated by isotopicDistribution().
    """
    
    import matplotlib.pyplot as pyt
    from numpy import arange
    
    while distribution[-1] < 0.001:
        distribution = distribution[:-1]
    
    pyt.bar(range(0, len(distribution) + 1), [0] + distribution, width = 0.3,
            align = 'center')
    pyt.xticks(arange(1, len(distribution)))
    pyt.yticks(arange(0, 110, 10))
    pyt.ylim((0, 110))
    pyt.show()
    
    if outputfile:
        pyt.savefig(outputfile)


def forPeptide(peptide_string):
    formula = renderFormula(peptide_string)
    return isotopicDistribution(formula)

def forFormula(formula_string):
    formdict = {}
    i = 0
    while i < len(formula_string):
        x = formula_string[i]
        if x.isalpha():
            try:
                y = formula_string[i+1]
                if y.isdigit():
                    formdict[x] = int(y)
                    i += 1
                else:
                    formdict[x] = 1                    
            except IndexError:
                formdict[x] = 1
            i += 1
        else:
            raise NotImplementedError, "Unable to parse chemical formula %s" % formula_string      
    return isotopicDistribution(formdict)
