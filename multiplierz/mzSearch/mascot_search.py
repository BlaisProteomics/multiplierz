from mascot import interface, report
from multiplierz.settings import settings
from multiplierz import myData
from datetime import datetime
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
    
    
    def __getitem__(self, field):
        return self.values[field]
    def __setitem__(self, field, value):
        self.values[field] = value
    
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
        assert os.path.exists(self.values['FILE']), (("Search target %s not "
                                                      "found!") % self.values['FILE'])
        
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
        
        if ':' in dat_id:
            dat_id, datestr = dat_id.split(':')
        else:
            datestr = datetime.now().strftime("%Y%m%d")
        
        output_dir = os.path.dirname(data_file)
        #resultfile = retrieveMascotReport(mascot_ids = [dat_id],
                                          #dates = [datestr],
                                          #mascot_server = mascot_server,
                                          #mascot_version = mascot_version,
                                          #login_name = user,
                                          #password = password,
                                          #ext = '.xlsx',
                                          #chosen_folder = output_dir,
                                          #ion_cutoff = bestType(self.values['ion_cutoff']),
                                          #show_query_data = bestType(self.values['show_query_data']) != 0,
                                          #show_same_set = bestType(self.values['show_same_set']) != 0,
                                          #show_sub_set = bestType(self.values['show_sub_set']) != 0,
                                          #rank_one = bestType(self.values['rank_one']) != 0,
                                          #bold_red = bestType(self.values['bold_red']) != 0,
                                          #) 
        #                      Other parameters should be externally suppliable!
        reportManager = report.MascotReport(mascot_server, mascot_version,
                                            user, password,
                                            cleanup = True)
        resultfile = reportManager.get_reports(mascot_ids = [dat_id],
                                               dates = [datestr],
                                               outputfile = self.values['FILE'] + '.xlsx',
                                               ion_cutoff = bestType(self.values['ion_cutoff']),
                                               show_query_data = bestType(self.values['show_query_data']) != 0,
                                               show_same_set = bestType(self.values['show_same_set']) != 0,
                                               show_sub_set = bestType(self.values['show_sub_set']) != 0,
                                               rank_one = bestType(self.values['rank_one']) != 0,
                                               bold_red = bestType(self.values['bold_red']) != 0)      
        if not isinstance(resultfile, basestring):
            assert len(resultfile) == 1, "Invalid result count: %s" % len(resultfile) # retrieveMascotReport can do multiple at a time; this doesn't.
            resultfile = resultfile[0]
        assert os.path.exists(resultfile), "Result file not present: %s" % resultfile
        return resultfile
        
        
        
        
        
        
        
        
        
        
        
def retrieveMascotReport(mascot_ids = None,
                         dat_file_list = None,
                         dates = None,
                         chosen_folder = myData,
                         mascot_server = settings.mascot_server,
                         mascot_version = settings.mascot_version,
                         combined_file = False,
                         rank_one = False,
                         max_hits = 1000,
                         ion_cutoff = 20,
                         bold_red = True,
                         show_query_data = True,
                         show_same_set = False,
                         show_sub_set = False,
                         protein_report = False,
                         quant = False,
                         ext = '.xlsx',
                         login_name = None,
                         password = None,
                         keep_dat = False,
                         pep2gene = None):
    """    
    Retrieves a search results from Mascot and performs several optional
    post-processing and annotation steps; use this to download reports 
    separately from the main MascotSearcher routine (e.g., to access old
    searches that hadn't been saved off-server, or to download reports
    with different output settings than the original run.)

    Note that if the requested output format (ext) is '.mzid' (mzIdentML), then 
    most of the optional arguments are nonfunctional (and the Mascot server 
    version must be >= 2.4.)

    Parameters:

    mascot_ids - A list of the Mascot search IDs that will be retreived.  (Should be
    left blank if dat_file_list is given.)

    dat_file_list - A list of paths to locally-located .DAT files. (Should be left blank
    if mascot_ids is given.)

    chosen_folder - Local directory path that the result files will be placed into.

    mascot_server - Address of the Mascot server containing the search results
    to be retrieved.

    mascot_version - Version number of the Mascot server; necessary due to small
    differences in the CGI protocol between versions.

    combined_file - If true, returns all search results in one file.  (CONFIRM?)

    rank_one - If true, search results will only contain first-rank peptides.

    max_hits - Maximum number of search hits per result file.

    ion_cutoff - (ASK?)

    bold_red - (ASK?)

    show_query_data - (ASK?)

    show_same_set - (ASK?)

    protein_report - (ASK?)

    quant - (ASK?)

    ext - Result file type; can be one of: .csv, .xls, .xlsx, .mzd, mzid . 
    (.xls and .xlsx require Microsoft Excel to be installed.)

    login_name - Username on the target Mascot server; required if security is
    enabled on the server.

    password - Password corresponding to given username.

    keep_dat - If true, the Mascot .DAT file is downloaded and not removed after
    processing.

    pep2gene - The path of a Pep2Gene database created by the Pep2Gene
    Utility, or None. If a database is provided, the Pep2Gene algorithm is
    used to annotate search result data with the genes corresponding to each
    peptide match.
    """
    import sys
    if max_hits > sys.maxint:
        print "Warning: max_hits must be of type int (for msparser.)"
        print "Reducing max_hits from %s to maximum int size (%s) ." % (max_hits, sys.maxint)
        max_hits = sys.maxint

    if mascot_ids and dat_file_list:
        raise NotImplementedError, ("Reports from server-located and local "
                                    "searches must be processed separately.")


    if mascot_ids:
        mascot_reporter = report.MascotReport(mascot_server,
                                              mascot_version,
                                              login_name,
                                              password,
                                              cleanup = not keep_dat)

        if any([':' in mid for mid in mascot_ids]):
            if not dates or not any(dates):
                dates = [mid.split(':')[1] if ':' in mid else None
                         for mid in mascot_ids]
            mascot_ids = [mid.split(':')[0] for mid in mascot_ids]
    else:
        mascot_reporter = report.MascotReport("No server",
                                              mascot_version)

    #Check Login
    if not login_name and password:
        mascot_reporter.mascot.logged_in = True
    else:
        if mascot_ids:
            check_login = mascot_reporter.login_mascot()
            if check_login == 'error':
                raise IOError, ("Could not log in to Mascot; check your Mascot settings "
                                "and username/password, if applicable.")
    
    # Currently always processing requests sequentially/separately;
    # get_reports for more than one report at once generates a file
    # combined by sheets.
    ret_vals = []
    for mascot_id in mascot_ids:
        ret_vals.append(mascot_reporter.get_reports(mascot_ids = [mascot_id],
                                                    dates = dates,
                                                    local_dat_files = dat_file_list,
                                                    chosen_folder = chosen_folder,
                                                    combined_file = combined_file,
                                                    rank_one = rank_one,
                                                    max_hits = max_hits,
                                                    ion_cutoff = ion_cutoff,
                                                    bold_red = bold_red,
                                                    show_query_data = show_query_data,
                                                    show_same_set = show_same_set,
                                                    show_sub_set = show_sub_set,
                                                    protein_report = protein_report,
                                                    quant = quant,
                                                    ext = ext,
                                                    mascotIDInResultName = True
                                                    )
                        )

    if pep2gene:
        p2gDB = pep2gene
        import multiplierz.mzTools.pep2gene as p2g
        for filename in ret_vals:
            p2g.add_gene_ids(filename, p2gDB, inPlace = True)

    return ret_vals