from multiplierz import __version__, SettingsFile, logger_message
from collections import defaultdict
import os








#if not os.path.exists(SettingsFile):
    #newfile = open(SettingsFile, 'w')
    #for line in defaultSettingsFile.split('\n'):
        #newfile.write(line)
    #newfile.close()


defaultKeys = [('logger verbosity', int),
               ('format default', str),
               ('image width', float),
               ('image height', float),
               ('XIC_view time', float),
               ('XIC_gen time', float),
               ('XIC_gen mz', float),
               ('MS1_view mz', float),
               ('sig_figs ms1_mz', int),
               ('sig_figs ms1_int', int),
               ('sig_figs ms2_mz', int),
               ('sig_figs ms2_int', int),
               ('sig_figs xic_time', int),
               ('sig_figs xic_int', int),
               ('ion_labels theor', bool),
               ('ion_labels error', bool),
               ('ion_labels figs', int),
               ('ion_labels units', str),
               ('mzServer use', str),
               ('comet directory', str),
               ('xtandem directory', str),
               ('mascot server', str),
               ('user email', str)]

mascotKeys = [('mascot var_mods', bool),
              ('mascot max hits', str),
              ('mascot ion cutoff', str),
              ('mascot bold red', bool),
              ('mascot show input query', bool),
              ('mascot rank one only', bool),
              ('mascot pep quant', bool),
              ('mascot retain DAT', bool),
              ('mascot sub-set hits', bool),
              ('mascot protein summary', bool),
              ('mascot p2g', bool),
              ('mascot mzid', bool),
              ('mascot same-set hits', bool)]

class MultiplierzSettings(object):
    '''Class to represent the settings file, with properties for getting and
    setting the various values, and default values in case the settings file gets
    corrupted somehow.'''
    
    def __init__(self, file_name, logger):
        self.file_name = file_name
        self._logger = logger
        
        self.load()

    
    def load(self):
        "Loads settings from the file."
        
        self.fields = {}
        with open(self.file_name, 'r') as file:
            for line in file:
                if not line.strip():
                    continue                
                if line.strip()[0] == '#':
                    continue
                if '=' not in line:
                    print(("Invalid settings line: %s" % line))
                    continue
                
                field, value = [x.strip() for x in line.split('=')]
                self.fields[field] = value
                
                
    def save(self):
        output = []
        seenFields = set()
        with open(self.file_name, 'r') as file:
            for line in file:
                if not line.strip():
                    output.append(line)
                    continue
                if line.strip()[0] == '#':
                    output.append(line)
                    continue
                if '=' not in line:
                    output.append(line)
                    continue       
                
                field, value = [x.strip() for x in line.split('=')]
                newval = self.fields[field]
                output.append('%s=%s\n' % (field, newval))
                seenFields.add(field)
        output.append('\n')
        
        with open(self.file_name, 'w') as file:
            for line in output:
                file.write(line)
            for field, value in list(self.fields.items()):
                if field in seenFields:
                    continue
                file.write("%s=%s\n" % (field, value))
        
    def get(self, field, typespec):
        if typespec == bool:
            try:
                return not (self.fields[field] == '0' or self.fields[field].lower() == 'false')
            except KeyError:
                assert field in [x[0] for x in mascotKeys]
                return False
        else:
            return typespec(self.fields[field])

            
            
    def set(self, field, value):
        self.fields[field] = str(value)
    
    
        
    def get_logger_level(self):
        return self.get('logger verbosity', int)
    def set_logger_level(self, verbosity):
        self.set('logger verbosity', verbosity)
    logger_level = property(get_logger_level, set_logger_level)
    
    def get_default_format(self):
        return self.get('format default', str) 
    def set_default_format(self, value):
        self.set('format default', value)
    default_format = property(get_default_format, set_default_format)
                           
    def get_image_size(self):
        return (self.get('image width', float), self.get('image height', float))
    def set_image_size(self, value):
        self.set('image width', value[0])
        self.set('image height', value[1])
    image_size = property(get_image_size, set_image_size)   
    
    def get_XIC_view_time_window(self):
        return self.get('XIC_view time', float)
    def set_XIC_view_time_window(self, value):
        self.set('XIC_view time', value)
    XIC_view_time_window = property(get_XIC_view_time_window, set_XIC_view_time_window)
    
    def get_XIC_gen_time_window(self):
        return self.get('XIC_gen time', float)
    def set_XIC_gen_time_window(self, value):
        self.set('XIC_gen time', value)
    XIC_gen_time_window = property(get_XIC_gen_time_window, set_XIC_gen_time_window)
    
    def get_XIC_gen_mz_window(self):
        return self.get('XIC_gen mz', float)
    def set_XIC_gen_mz_window(self, value):
        self.set('XIC_gen mz', value)
    XIC_gen_mz_window = property(get_XIC_gen_mz_window, set_XIC_gen_mz_window) 
    
    def get_MS1_view_mz_window(self):
        return self.get('MS1_view mz', float)
    def set_MS1_view_mz_window(self, value):
        self.set('MS1_view mz', value)
    MS1_view_mz_window = property(get_MS1_view_mz_window, set_MS1_view_mz_window)
    
    def get_ms1_mz_figs(self):
        return self.get('sig_figs ms1_mz', int)
    def set_ms1_mz_figs(self, value):
        self.set('sig_figs ms1_mz', value)
    ms1_mz_figs = property(get_ms1_mz_figs, set_ms1_mz_figs)
    
    def get_ms1_int_figs(self):
        return self.get('sig_figs ms1_int', int)
    def set_ms1_int_figs(self, value):
        self.set('sig_figs ms1_int', value)
    ms1_int_figs = property(get_ms1_int_figs, set_ms1_int_figs)
    
    def get_ms2_mz_figs(self):
        return self.get('sig_figs ms2_mz', int)
    def set_ms2_mz_figs(self, value):
        self.set('sig_figs ms2_mz', value)
    ms2_mz_figs = property(get_ms2_mz_figs, set_ms2_mz_figs)
    
    def get_ms2_int_figs(self):
        return self.get('sig_figs ms2_int', int)
    def set_ms2_int_figs(self, value):
        self.set('sig_figs ms2_int', value)
    ms2_int_figs = property(get_ms2_int_figs, set_ms2_int_figs)
    
    def get_xic_time_figs(self):
        return self.get('sig_figs xic_time', int)
    def set_xic_time_figs(self, value):
        self.set('sig_figs xic_time', value)
    xic_time_figs = property(get_xic_time_figs, set_xic_time_figs)
    
    def get_xic_int_figs(self):
        return self.get('sig_figs xic_int', int)
    def set_xic_int_figs(self, value):
        self.set('sig_figs xic_int', value)
    xic_int_figs = property(get_xic_int_figs, set_xic_int_figs)
    
    def get_show_theor_mz(self):
        return self.get('ion_labels theor', bool)
    def set_show_theor_mz(self, value):
        self.set('ion_labels theor', value)
    show_theor_mz = property(get_show_theor_mz, set_show_theor_mz)
    
    def get_show_mass_error(self):
        return self.get('ion_labels error', bool)
    def set_show_mass_error(self, value):
        self.set('ion_labels error', value)
    show_mass_error = property(get_show_mass_error, set_show_mass_error)
    
    def get_mass_error_figs(self):
        return self.get('ion_labels figs', int)
    def set_mass_error_figs(self, value):
        self.set('ion_labels figs', value)
    mass_error_figs = property(get_mass_error_figs, set_mass_error_figs)
    
    def get_mass_error_units(self):
        return self.get('ion_labels units', str)
    def set_mass_error_units(self, value):
        self.set('ion_labels units', value)
    mass_error_units = property(get_mass_error_units, set_mass_error_units)
    
    
    def get_mzServer(self):
        return self.get('mzServer use', str)
    def set_mzServer(self, value):
        self.set('mzServer use', value)
    mzServer = property(get_mzServer, set_mzServer)    
    
    def get_mascot_server(self):
        return self.get('mascot server', str)
    def set_mascot_server(self, value):
        self.set('mascot server', value)
    mascot_server = property(get_mascot_server, set_mascot_server)
    
    def get_mascot_version(self):
        return self.get('mascot version', str)
    def set_mascot_version(self, value):
        self.set('mascot version', value)
    mascot_version = property(get_mascot_version, set_mascot_version)
    
    def get_mascot_ms2(self):
        return self.get('mascot ms2', bool)
    def set_mascot_ms2(self, value):
        self.set('mascot ms2', value)
    mascot_ms2 = property(get_mascot_ms2, set_mascot_ms2)
    
    def get_mascot_security(self):
        return self.get('mascot security', bool)
    def set_mascot_security(self, value):
        self.set('mascot security', value)
    mascot_security = property(get_mascot_security, set_mascot_security)
    
    def get_mascot_var_mods(self):
        return self.get('mascot var_mods', bool)
    def set_mascot_var_mods(self, value):
        self.set('mascot var_mods', value)
    mascot_var_mods = property(get_mascot_var_mods, set_mascot_var_mods)
    
    def retrieve_mascot_report_controls(self):
        defaultCtrl = {}
        #defaultCtrl['max hits'] = self.get('mascot max hits', str)
        #defaultCtrl['ion cutoff'] = self.get('mascot ion cutoff', str)
        #defaultCtrl['bold red'] = self.get('mascot bold red', bool)
        #defaultCtrl['show input query'] = self.get('mascot show input query', bool)
        for field, type in mascotKeys:
            defaultCtrl[field[7:]] = self.get(field, type)
        return defaultCtrl
        
    def save_mascot_report_controls(self, controls):
        #self.set('mascot max hits', controls['max hits'])
        #self.set('mascot ion cutoff', controls['ion cutoff'])
        #self.set('mascot bold red', controls['bold red'])
        #self.set('mascot show input query', controls['show input query'])
        for key, value in list(controls.items()):
            self.set('mascot ' + key, value)
        self.save()
        
    def get_comet(self):
        return os.path.normpath(self.fields['comet directory'])
    def set_comet(self, value):
        self.fields['comet directory'] = value
    comet = property(get_comet, set_comet)
    
    def get_xtandem(self):
        return os.path.normpath(self.fields['xtandem directory'])
    def set_xtandem(self, value):
        self.fields['xtandem directory'] = value
    xtandem = property(get_xtandem, set_xtandem)
 
    def get_user_email(self):
        return self.fields['user email']
    def set_user_email(self, value):
        self.fields['user email'] = value
    user_email = property(get_user_email, set_user_email)
    

settings = MultiplierzSettings(SettingsFile, logger_message)
#settings.save()
#print "FOO"