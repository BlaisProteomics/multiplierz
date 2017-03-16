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

# this module is for accessing the unimod database
# (which we'll include with Multiplierz)

# customizations are a bit tricky--they require modifying the unimod.xml file.
# the easiest way to do that is to use a Mascot server's configuration editor
# to add the mods to its own unimod file, and then copy that file into the
# multiplierz directory.

#try:
    #import lxml.etree as ET
#except ImportError:
    #print "Failed to import lxml!  Using xml.etree.cElementTree, which may be slower."
    #import xml.etree.cElementTree as ET
import lxml.etree as ET
from collections import defaultdict

import xml.etree.cElementTree as ETT

# namespace dictionary for the unimod schema. 'f' is just to distinguish
# the 'lower' function, using the Blais URL is not really necessary
NS = {'umod': 'http://www.unimod.org/xmlns/schema/unimod_2',
      'f': 'http://blais.dfci.harvard.edu'}

def lower(context, a):
    return a[0].lower()


# For pyComet lookup table.  Remove when obsolete!
siteTypeLookup = {'G' : 'add_G_glycine',
                  'A' : 'add_A_alanine',
                  'S' : 'add_S_serine',
                  'P' : 'add_P_proline',
                  'V' : 'add_V_valine',
                  'T' : 'add_T_threonine',
                  'C' : 'add_C_cysteine',
                  'L' : 'add_L_leucine',
                  'I' : 'add_I_isoleucine',
                  'N' : 'add_N_asparagine',
                  'D' : 'add_D_aspartic_acid',
                  'Q' : 'add_Q_glutamine',
                  'K' : 'add_K_lysine',
                  'E' : 'add_E_glutamic_acid',
                  'M' : 'add_M_methionine',
                  'O' : 'add_O_ornithine',
                  'H' : 'add_H_histidine',
                  'F' : 'add_F_phenylalanine',
                  'R' : 'add_R_arginine',
                  'Y' : 'add_Y_tyrosine',
                  'W' : 'add_W_tryptophan',
                  'N-term' : 'add_Nterm_peptide',
                  'C-term' : 'add_Cterm_peptide'}

class LookupError(Exception):
    '''Exception thrown when a modification does not exist. Also
    thrown if a mod does not have a given specificity.'''

    def __init__(self, mod, site=None):
        self.mod = mod
        self.site = site

    def __str__(self):
        if self.site:
            return 'Modification %s at site %s was not found' % (self.mod, self.site)
        else:
            return 'Modification %s was not found' % self.mod

class UnimodDatabase(object):
    '''Class to represent the Unimod XML file, with some methods
    for looking up modification deltas, specificities, etc, as well
    as things like element masses and amino acid masses. All masses
    are monoisotopic by default.'''

    def __init__(self, file_name):
        self.file_name = file_name

        # this is a special class for evaluating XPath expressions quickly.
        # we have to give it the extra function 'lower' because lxml doesn't
        # have xpath 2.0 functions (yet)
        #self.utree = ET.XPathEvaluator(ET.parse(self.file_name),
                                       #namespaces=NS,
                                       #extensions={(NS['f'], 'lower'): lower})
        self.tree = ET.parse(self.file_name)
        evaluator = ET.XPathEvaluator(self.tree,
                                       namespaces=NS,
                                       extensions={(NS['f'], 'lower'): lower})
        
        #etreeeval = ETT.parse(self.file_name).findall        
        #def doublecheck(expression):
            #lxmlbased = evaluator(expression)
            #etreebased = etreeeval(expression, NS)
            #assert map(lambda x: x.attrib, lxmlbased) == map(lambda x: x.attrib, etreebased)
            #return lxmlbased
        
        #self.utree = doublecheck
        self.utree = evaluator
        
        
    


    def get_amino_acid_composition(self, aa):
        '''Gets the amino acid composition. 'aa' should be the single-letter abbreviation
        for an amino acid.'''

        aa_list = self.utree('.//umod:amino_acids/umod:aa[@title="%s"]' % aa)
        new_aa_list = self.tree.findall('umod:amino_acids', NS)
        

        if aa_list:
            return dict((e.get('symbol'), int(e.get('number')))
                        for e in aa_list[0].xpath('.//umod:element[@symbol][@number]',
                                                  namespaces=NS))
        else:
            return None






    def mod_exists(self, mod_name, site=None):
        '''Searches for a mod by name--which name applies depends on the mod
        (see unimod.org help for details). Returns True is the mod is found.

        'site' is an optional argument--if given, the function only returns true
        if the modification exists at the given site (AA code or 'N|C-term').
        '''

        modifications = self.utree('.//umod:modifications')[0]
        matches = [x for x in modifications if
                   mod_name.lower() in [x.get('title').lower(), x.get('full_name').lower()]]
        
        if not site:
            return bool(matches)
        else:
            specificities = [x for x in matches[0] if 'specificity' in x.tag]
            return any([site in [x.get('site'), x.get('position')] for x in specificities])
        

    def get_mod_delta(self, mod_name, mass='monoisotopic'):
        '''Searches for a mod by name--which name applies depends on the mod
        (see unimod.org help for details). Returns monoisotopic or average mass.

        'mass' should be 'monoisotopic' or 'average'. Defaults to monoisotopic mass.'''

        modifications = self.utree('.//umod:modifications')[0]
        matches = [x for x in modifications if
                   mod_name.lower() in [x.get('title').lower(), x.get('full_name').lower()]]        
        
        #assert len(matches) == 1, str([x.attrib for x in matches])
        if len(matches) > 1:
            print "Ambiguous entries for %s" % mod_name
        
        delta = [x for x in matches[0] if 'delta' in x.tag][0]
        
        if mass == 'monoisotopic':
            return float(delta.get('mono_mass'))
        elif mass == 'average':
            return float(delta.get('avge_mass'))
        else:
            raise ValueError("%s is not a valid argument for mass type (should be 'monoisotopic' or 'average'.)")
        

    def get_mod_specificities(self, mod_name):
        '''Searches for a mod by name and returns a dictionary of specificities, organized by
        group number. Keys are integers and values are lists of AAs'''

        modifications = self.utree('.//umod:modifications')[0]
        matches = [x for x in modifications if
                   mod_name.lower() in [x.get('title').lower(), x.get('full_name').lower()]]        
        
        #assert len(matches) <= 1
        if len(matches) > 1:
            print "Ambiguous entries for %s" % mod_name        
        
        if not matches:
            raise ValueError
        
        specificities = defaultdict(list)
        for match in [x for x in matches[0] if 'specificity' in x.tag]:
            specificities[match.get('spec_group')].append((match.get('site'), match.get('position')))
        return specificities


    def get_mod_neutral_loss(self, mod_name, site, mass='monoisotopic'):
        '''For a given modification and a site specificity, returns the neutral
        loss list if possible. This returns None due to the mod not existing,
        or the site not existing, and it returns an empty tuple if the site has
        no neutral losses.

        'mass' should be 'monoisotopic' or 'average'. Defaults to monoisotopic mass.'''

        specificity_list = self.utree(('.//umod:modifications/'
                                       'umod:mod[f:lower(@title) = "%s"]/'
                                       'umod:specificity[@site="%s"]') % (mod_name.lower(), site))
        specificity_list += self.utree(('.//umod:modifications/'
                                       'umod:mod[f:lower(@full_name) = "%s"]/'
                                       'umod:specificity[@site="%s"]') % (mod_name.lower(), site))

        if not specificity_list and site in ('N-term', 'C-term'):
            specificity_list = self.utree(('.//umod:modifications/'
                                           'umod:mod[f:lower(@title) = "%s"]/'
                                           'umod:specificity[fn:ends-with(f:lower(@position), "%s")]') % (mod_name.lower(), site))
            specificity_list += self.utree(('.//umod:modifications/'
                                           'umod:mod[f:lower(@full_name) = "%s"]/'
                                           'umod:specificity[fn:ends-with(f:lower(@position), "%s")]') % (mod_name.lower(), site))            

        if specificity_list:
            if mass == 'monoisotopic':
                foo = tuple(float(nl.get('mono_mass')) for s in specificity_list
                             for nl in s.xpath('.//umod:NeutralLoss[@mono_mass]',
                                               namespaces=NS))
            elif mass == 'average':
                foo = tuple(float(nl.get('avge_mass')) for s in specificity_list
                             for nl in s.xpath('.//umod:NeutralLoss[@avge_mass]',
                                               namespaces=NS))
            else:
                raise ValueError("%s is not a valid argument for mass type, "
                                 "should be 'monoisotopic' or 'average'." % mass)
            return foo
        else:
            raise LookupError(mod_name, site)

    def get_mod_composition(self, mod_name):
        '''Searches for a mod, and returns the atomic composition of the delta
        as a dictionary of element symbol -> number pairs'''

        mod_list = self.utree(('.//umod:mod[f:lower(@title)="%s"]/'
                               'umod:delta[@mono_mass]') % mod_name.lower())
        if not mod_list:
            mod_list = self.utree(('.//umod:mod[f:lower(@full_name)="%s"]/'
                                   'umod:delta[@mono_mass]') % mod_name.lower())            

        if mod_list:
            return dict((e.get('symbol'), int(e.get('number')))
                        for e in mod_list[0].xpath('.//umod:element[@symbol][@number]',
                                                   namespaces=NS))
        else:
            raise LookupError(mod_name)
        
    
    def get_all_mod_names(self, both = True):
        modifications = self.utree('.//umod:modifications')[0]
        if both:
            return [(x.attrib.get('full_name', None), x.attrib['title'])
                    for x in modifications]
        else:
            return [x.attrib['full_name'] if x.attrib['full_name'] else x.attrib['title']
                    for x in modifications]
    
    def get_mod_weight_lookup(self):
        modifications = self.utree('.//umod:modifications')[0]
        lookup = {}
        for mod in modifications:
            modName = mod.attrib['full_name']
            delta = mod.iterfind('{http://www.unimod.org/xmlns/schema/unimod_2}delta').next()
            modMass = delta.attrib['avge_mass']
            lookup[modName] = modMass
            
        return lookup
    

    def get_pycomet_lookup(self):
        modifications = self.utree('.//umod:modifications')[0]
        lookup = {}
        for mod in modifications:
            modName = mod.attrib['full_name']
            delta = mod.iterfind('{http://www.unimod.org/xmlns/schema/unimod_2}delta').next()
            modMass = delta.attrib['avge_mass']

            massLookup = {}
            specificities = mod.iterfind('{http://www.unimod.org/xmlns/schema/unimod_2}specificity')
            for specificity in specificities:
                site = specificity.attrib['site']
                longSite = siteTypeLookup[site]
                massLookup[longSite] = modMass
            
            lookup[modName] = massLookup
        
        return lookup
                
    
