from mascot import interface, report
from multiplierz.settings import settings
import os

# This is a "user-friendly" interface to the multiplierz Mascot-search
# capabilities; the internals of the process are located in the
# mascot submodule.



def bestType(value):
    try:
        try:
            return int(value)
        except ValueError:
            return float(value)
    except ValueError:
        return str(value)

class MascotSearch(object):
    def __init__(self, parameters):
        self.fields = []
        self.values = {}
        with open(parameters, 'r') as par:
            for line in par:
                if '=' in line:
                    words = [x.strip() for x in line.split('=')]
                    if len(words) != 2:
                        print "ERROR on line: %s" % line
                    else:
                        self.fields.append(words[0])
                        self.values[words[0]] = bestType(words[1])
    
    
    #def __getattr__(self, field):
        #return self.values[field]
    
    def write(self, filename):
        with open(filename, 'w') as output:
            for field in self.fields:
                output.write('%s\%s' % (field, self.values[field]))
            
            for field in [x for x in self.values.keys() if x not in self.fields]:
                output.write('%s\%s' % (field, self.values[field]))
    
    
    def run_search(self, data_file = None, user = None, password = None):
        if data_file:
            self.values['FILE'] = data_file
        assert 'FILE' in self.values, "Target data must be specified!"
        
        if not self.values['FILE'].lower().endswith('mgf'):
            # This will perform extraction, but currently without much control.
            # Add this to parameters files?
            from multiplierz.mgf import extract
            self.values['FILE'] = extract(self.values['FILE'],
                                          default_charge = self.values['CHARGE'])
            
        if 'mascot_server' in self.values:
            mascot_server = self.values['mascot_server']
            mascot_version = self.values['mascot_version']
        else:
            mascot_server = settings.mascot_server
            mascot_version = settings.mascot_version
            
        search = interface.MascotSearcher(mascot_server, version = mascot_version)
        if user and password:
            search.login(user, password)
        dat_id, error = search.search(self.values)
        if error:
            raise Exception, error
        assert dat_id, "No dat_id; no error raised. %s" % error
        
        output_dir = os.path.dirname(data_file)
        resultfile = report.retrieveMascotReport(mascot_ids = [dat_id],
                                                 mascot_server = mascot_server,
                                                 mascot_version = mascot_version,
                                                 login_name = user,
                                                 password = password,
                                                 ext = '.xlsx',
                                                 chosen_folder = output_dir,
                                                 ion_cutoff = bestType(self.values['ion_cutoff']),
                                                 show_query_data = bestType(self.values['show_query_data']) != 0,
                                                 show_same_set = bestType(self.values['show_same_set']) != 0,
                                                 show_sub_set = bestType(self.values['show_sub_set']) != 0,
                                                 rank_one = bestType(self.values['rank_one']) != 0,
                                                 bold_red = bestType(self.values['bold_red']) != 0,
                                                 ) 
        #                      Other parameters should be externally suppliable!
        
        assert len(resultfile) == 1, "Invalid result count: %s" % len(resultfile) # retrieveMascotReport can do multiple at a time; this doesn't.
        resultfile = resultfile[0]
        assert os.path.exists(resultfile), "Result file not present: %s" % resultfile
        return resultfile
        