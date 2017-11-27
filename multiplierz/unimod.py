import sqlite3
import base64
import cPickle as pickle


def marshal(data):
    return base64.b64encode(pickle.dumps(data, protocol = 2))
def demarshal(data):
    return pickle.loads(base64.b64decode(data))


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
                  'C-term' : 'add_Cterm_peptide',
                  'U' : 'add_U_user_amino_acid' # Not actually used.
                  }

class UnimodDatabase(object):
    def __init__(self, file_name):
        self.file_name = file_name
        
        self.con = sqlite3.connect(file_name)
        self.cur = self.con.cursor()
    
    def mod_exists(self, mod_name, site = None):
        mod_name = mod_name.lower()
        
        if site:
            command = """SELECT EXISTS(SELECT 1 FROM modsites WHERE
            (name ='%s' OR full_name='%s') AND site = '%s')""" % (mod_name, mod_name, site)
        else:
            command = """SELECT EXISTS(SELECT 1 FROM mods WHERE
            name ='%s' OR full_name='%s')""" % (mod_name, mod_name)
        self.cur.execute(command)
        return bool(self.cur.fetchall()[0]) #?
    
    def get_mod_delta(self, mod_name, mass = 'monoisotopic'):
        mod_name = mod_name.lower()
        
        if mass == 'monoisotopic':
            command = """SELECT mon_delta FROM mods WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        elif mass == 'average':
            command = """SELECT avg_delta FROM mods WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        else:
            raise ValueError, "Mass type must be 'monoisotopic' or 'average'."
        self.cur.execute(command)
        try:
            return self.cur.fetchall()[0][0]
        except IndexError:
            raise ValueError, "%s not found in Unimod database." % mod_name
    
    def get_mod_specificities(self, mod_name):
        command = """SELECT spec FROM mods WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        self.cur.execute(command)
        try:
            spec = self.cur.fetchall()[0][0]
        except IndexError:
            raise ValueError, "%s not found in Unimod database." % mod_name
        
        return demarshal(spec)
    
    def get_mod_neutral_loss(self, mod_name, site, mass = 'monoisotopic'):
        mod_name = mod_name.lower()
        
        if mass == 'monoisotopic':
            command = """SELECT mon_n_loss FROM modsites WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        elif mass == 'average':
            command = """SELECT avg_n_loss FROM modsites WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        else:
            raise ValueError, "Mass type must be 'monoisotopic' or 'average'."            
        self.cur.execute(command)
        return demarshal(self.cur.fetchall()[0][0])
    
    def get_mod_composition(self, mod_name):
        mod_name = mod_name.lower()
        
        command = """SELECT comp FROM mods WHERE name = '%s' OR full_name = '%s'""" % (mod_name, mod_name)
        self.cur.execute(command)
        return demarshal(self.cur.fetchall()[0][0])
    
    def get_all_mod_names(self):
        command = """SELECT name FROM mods"""
        self.cur.execute(command)
        return [x[0] for x in self.cur.fetchall()]
    
    def get_all_mod_titles(self):
        command = """SELECT full_name FROM mods"""
        self.cur.execute(command)
        return [x[0] for x in self.cur.fetchall()]    
    
    def site_form_mod_names(self):
        from multiplierz.internalAlgorithms import collectByCriterion
        command = """SELECT name, site, specset FROM modsites"""
        self.cur.execute(command)
        specmoddata = self.cur.fetchall()
        specmodsets = collectByCriterion(specmoddata, lambda x: (x[0], x[2]))
        modnames = []
        for specmodset in specmodsets.values():
            modname = specmodset[0][0]
            sites = ''.join(sorted(set(zip(*specmodset)[1])))
            modnames.append('%s (%s)' % (modname, sites))
        return sorted(modnames)
            
    def get_pycomet_lookup(self):
        from multiplierz.internalAlgorithms import collectByCriterion
        command = """SELECT name, site, specset FROM modsites"""
        self.cur.execute(command)
        specmoddata = self.cur.fetchall()
        specmodsets = collectByCriterion(specmoddata, lambda x: (x[0], x[2]))
        lookup = {}
        for specmodset in specmodsets.values():
            modname = specmodset[0][0]
            delta = self.get_mod_delta(modname)
            sites = sorted(set(zip(*specmodset)[1]))
            massLookup = {}
            for site in sites:
                massLookup[siteTypeLookup[site]] = delta
            lookup[modname] = massLookup
        return lookup