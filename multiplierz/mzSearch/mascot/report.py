

import os
import shutil

from tempfile import mkstemp
from collections import defaultdict

import multiplierz.mzTools.mz_image as mz_image

#import multiplierz.mzTools.precursor_peaks
# Plots in Excel files are REALLY REALLY not supported anymore.



#import multiplierz.mzSearch.mascot.interface as interface # This should be made more proper!
import interface
import multiplierz.mzReport as mzReport

from multiplierz import myData, myTemp, logger_message
from multiplierz.settings import settings





class MascotReport:
    def __init__(self, server=None, version=None, login=None, password=None, cleanup=True):
        if server:
            self.mascot = interface.mascot(server, version)
            self.login = login
            self.password = password
            #set to False to leave dat file in place for further queries
            self.cleanup = cleanup
            self.mascot.login(login, password)
        else:
            self.mascot = interface.mascot()
            self.cleanup = cleanup
            print "Note: No server specified, Mascot-independent mode."

    def login_mascot(self, *args, **kwargs):
        print "Deprecated function login_mascot called."
        pass

    def mascot_header(self, report_file, mascot_header):
        '''Adds a Mascot Header page to a report (XLS or MZD)'''

        logger_message(30, 'Adding Mascot Header...')

        if report_file.lower().endswith('.mzd'):
            mascot_rep = mzReport.writer(report_file, columns=mascot_header[0], table_name='MascotHeader')
        else:
            mascot_rep = mzReport.writer(report_file, columns=mascot_header[0], sheet_name='Mascot_Header')

        for line in mascot_header[1:]:
            mascot_rep.write(line)

        mascot_rep.close()

        logger_message(20, 'Mascot Header Complete')

    def mascot_headers(self, report_file, mascot_headers):
        '''Combines a list of Mascot header pages into a single page (XLS or MZD)'''

        logger_message(30, 'Adding Mascot Headers...')

        first_columns = set()
        columns = defaultdict(dict)

        report_files = [f for f,mh in mascot_headers]

        for f,mascot_header in mascot_headers:
            first_col = []

            for line in mascot_header[1:]:
                if line[0] == ' ' or (isinstance(line[1], (str,unicode))
                                      and line[1].startswith('-----')):
                    continue
                columns[f][line[0]] = line[1]
                first_col.append(line[0])

            first_columns.add(tuple(first_col))

        if len(first_columns) > 1:
            logger_message(20, 'Headers seems different, will try to merge them')

            main_h = list(max(first_columns, key=len))
            main_h += sorted(reduce(set.union, first_columns, set()).difference(main_h))
        else:
            main_h = list(first_columns.pop())

        cols = ['Header'] + [(os.path.basename(f) if f else ('-' * 50))
                             for f in report_files]

        if report_file.lower().endswith('.mzd'):
            mascot_rep = mzReport.writer(report_file, columns=cols,
                                         table_name='MascotHeader')
        else:
            mascot_rep = mzReport.writer(report_file, columns=cols,
                                         sheet_name='Mascot_Header')

        for col in main_h:
            row = [col]

            for f in report_files:
                if col in columns[f]:
                    row.append(columns[f][col])
                else:
                    row.append(None)

            mascot_rep.write(row)

        mascot_rep.close()

    def mascot_prot_coverage(self, mascot_id, ion_cutoff, date=None):
        logger_message(30,'Downloading Mascot Protein Coverage...')

        prots = {}

        row = (yield None)


        #for row in dataFile:
        while row:
            acc = row['Accession Number']
            start_pos = row['Start Position']
            end_pos = row['End Position']
            db_idx = row.get('protein database', None)
            if db_idx:
                db_idx = db_idx.split('::')[0]

            if acc in prots:
                (coverage, start_end, cov_per) = prots[acc]
            else:
                logger_message(10, "Getting coverage for: %s" % acc)
                (coverage, start_end, cov_per) = self.mascot.get_coverage(acc, mascot_id=mascot_id, date=date, cutoff=ion_cutoff, db_idx=db_idx)
                prots[acc] = (coverage, start_end, cov_per)

            #check if connection error
            if coverage == 'error':
                logger_message(30,'Mascot Login Error. Please check your username and password and try again')
                fh = open('multiplierz_error.txt','a')
                fh.write("Mascot Login Error. Please check your username and password and try again\n")
                fh.close()
                break

            # convert protein-index into text-index
            conv_x = lambda i: i + 6 + ((i - 1) / 50) * 11 + ((i - 1) % 50) / 10

            nse = self.mascot.get_range([conv_x(i) for i in range(start_pos, end_pos + 1)])
            metadatatuple = (coverage, start_end, nse)

            row = (yield (cov_per, ('Accession Number', 'prot coverage', metadatatuple)))

    def mascot_web(self, mascot_id, ms2_img, mascot_var_mods=True,
                   instrument='ESI-TRAP', date=None,
                   isMZD=False, dat_file=None, im_size=(8.0,6.0)):

        logger_message(30, 'Downloading Mascot Spectra...')

        if mascot_var_mods and not ms2_img:
            mods_only = True
        if ms2_img:
            mods_only = False
        if mascot_var_mods and self.mascot.version >= '2.2':
            logger_message(40, 'Warning: Variable modifications already present, overwriting')

        if dat_file:
            info_gen = interface.mascot_ms2(dat_file)
            info_gen.next()

        row = (yield None)


        while row:
            query = row['Query']
            logger_message(10, 'Downloading Mascot Spectra for query: %s' % query)

            peptide = row['Peptide Sequence']
            charge = row['Charge']
            score = row['Peptide Score']
            pep_rank = row['Peptide Rank']

            if dat_file:
                (peptide, var_text, mass_int, mass_labels, ion_list) = info_gen.send((query, pep_rank))

                mass_int = tuple(sorted(mass_int))

                if ms2_img:
                    if isMZD:
                        image_tup = ('Peptide Sequence', 'ms2', (mass_int, 'c', peptide, mass_labels,
                                                                 ion_list, charge, score))
                    else:
                        (h, img_file) = mkstemp(suffix='.png', prefix='mascot_', dir=myTemp)
                        os.close(h)

                        logger_message(20, 'Drawing MS MS Mass Plot...')

                        mz_image.make_ms2_im(img_file[:-4], mass_int, 'c', peptide,
                                             mass_labels, ion_list, charge, score, im_size=im_size)

                        image_tup = ('Peptide Sequence', 'image', img_file)
                else:
                    image_tup = None
            else:
                (h, img_file) = mkstemp(suffix='.gif', prefix='mascot_', dir=myTemp)
                os.close(h)

                (var_text,ionInfo,phosphoPos,oxidPos) = self.mascot.get_msms(query, img_file, mascot_id=mascot_id,
                                                                             score=score, pep_rank=pep_rank,
                                                                             date=date, sequence=peptide,
                                                                             mod_only=mods_only, dat_file=dat_file)

                if ms2_img:
                    (h, gifpath) = mkstemp(suffix='.png', prefix='combo_', dir=myTemp)
                    os.close(h)

                    mz_image.make_mascot_ms2_im(gifpath, img_file, peptide, charge, score,
                                                ionInfo, phosphoPos, oxidPos, instrument)

                    image_tup = ('Peptide Sequence', 'image', gifpath)
                else:
                    image_tup = None

            if not mascot_var_mods:
                var_text = None

            row = (yield (var_text, image_tup))


    def prot_report(self, report_file, prot_report):
        '''Adds a protein page to a report (XLS or MZD)'''

        logger_message(30, 'Adding Protein Info...')

        if report_file.lower().endswith('.mzd'):
            prot_rep = mzReport.mzDB.sqlite3.connect(report_file)
            prot_rep.execute('create view ProteinData as select '
                             '"Protein Rank","Accession Number","Protein Description",'
                             '"Protein Mass","Protein Matches","Protein Score",'
                             'count(distinct "Peptide Sequence") as "Unique Peptides"'
                             ' from PeptideData group by "Protein Rank","Accession Number"')

            prot_rep.close()
        else:
            cols = ['Protein Rank', 'Accession Number', 'Protein Description',
                    'Protein Mass', 'Protein Matches', 'Protein Score', 'Unique Peptides']

            prot_rep = mzReport.writer(report_file, columns=cols, sheet_name='Protein')

            for line in prot_report:
                prot_rep.write(line)

            prot_rep.close()

    #def get_report(self, mascot_id, date=None, chosen_folder=None, 
                   #local_dat_file = None, mascotIDInResultName = False, **kwargs):
        #return self.get_reports(mascot_ids=[mascot_id], dates=[date] if date else None,
                                #chosen_folder=chosen_folder,
                                #combined_file=False, 
                                #local_dat_files = [local_dat_file] if local_dat_file else None,
                                #mascotIDInResultName = mascotIDInResultName,
                                #**kwargs)[0]

    def legacy_get_reports(self, mascot_ids, dates=None,
                    chosen_folder=None, combined_file=False, rank_one=False,
                    protein_report=False,
                    mascot_options=None,
                    peaks=False, peaks_options=None,
                    mascot_web=False, mascot_web_options=None,
                    mascot_prot_cov=False, ext='.xlsx',
                    local_dat_files = None,
                    mascotIDInResultName = False,
                    percolatorDirectory = None,
                    **kwargs):

        # mascot_ids should be a list/tuple of IDs. dates should be a matching list/tuple of dates,
        # or False. combined_file should be None for individual files or an output file name

        # mascot options: (max_hits, ion_cutoff, bold_red, unassigned_queries,
        #                  show_query_data, show_same_set, show_sub_set, quant) + mascot_id, date
        # mascot_web options: (ms2_img, mascot_ms2, mascot_var_mods,
        #                      draw_pep, instrument, im_size) + mascot_id, date
        # mascot_prot_cov options: ion_cutoff, mascot_id, date

        # defaults and overrides. The priority is:  keyword > option_dict > default

        # Using local .DATs means you don't have access to certain fancy
        # Mascot features.
        if local_dat_files:
            mascot_web = False
            mascot_prot_cov = False



        assert not peaks, ("precursor_peaks and images in result files are no longer supported; "
                           "peaks argument to get_reports must be False.")

        # defaults
        _mascot_options = dict(max_hits=1000, ion_cutoff=20, bold_red=True,
                               unassigned_queries=False, show_query_data=True,
                               show_same_set=False, show_sub_set=False, quant=False)
        # option_dict
        if mascot_options:
            _mascot_options.update(mascot_options)
        # keywords
        _mascot_options.update((k,kwargs[k]) for k in kwargs if k in _mascot_options)
        for k in _mascot_options:
            if k in kwargs:
                _mascot_options[k] = kwargs[k]

        if peaks:
            # defaults
            _peaks_options = dict(time_window=(0.5,0.5), mz_window = (0.1,0.1),
                                  plot_ms1=False, plot_xic=False, plot_ms2=False,
                                  peak_area=False, reporter_ions=False, peakfilter=None,
                                  ion_list=['b','y'], instrument='ESI-TRAP', im_size=(8.0,6.0))
            # option_dict
            if peaks_options:
                _peaks_options.update(peaks_options)
            # keywords
            _peaks_options.update((k,kwargs[k]) for k in kwargs if (k in _peaks_options
                                                                    or k == 'peak_data_path'))

            # need a path (file or directory) to actually do this,
            # so we raise an exception if it's not present
            if 'peak_data_path' not in _peaks_options:
                raise ValueError('peak_data_path value is required for peak extraction')

        if mascot_web:
            # defaults
            _mascot_web_options = dict(ms2_img=True, mascot_ms2=True,
                                       mascot_var_mods=True, instrument='ESI-TRAP', im_size=(8.0,6.0))
            # option_dict
            if mascot_web_options:
                _mascot_web_options.update(mascot_web_options)
            # keywords
            _mascot_web_options.update((k,kwargs[k]) for k in kwargs if k in _mascot_web_options)

        # if version is 2.2+, mod positions are extracted automatically
        if mascot_web and self.mascot.version >= '2.2':
            _mascot_web_options['mascot_var_mods'] = False
            if not _mascot_web_options['ms2_img']:
                mascot_web = False

        # require agreement 'instrument' between two dictionaries
        if peaks and mascot_web:
            if _peaks_options['instrument'] != _mascot_web_options['instrument']:
                raise ValueError('instrument value must be consistent; input dictionaries disagree')

        # Getting both of these would be redundant, so force at most one
        if mascot_web and peaks and _mascot_web_options['ms2_img']:
            _peaks_options['plot_ms2'] = False

        if chosen_folder is None:
            chosen_folder = myData


        # if creating a single file, we'll create the writer now
        if combined_file:
            # figuring out the report columns. start with defaults...
            repcols = mzReport.default_columns[:]

            # mascot 2.3 can have multiple databases so add a column for that
            if self.mascot.version >= '2.3':
                repcols.insert(1, 'Protein Database')

            # these are the columns coming out of the dat file, need them separate
            res_cols = repcols[:]

            # add columns for peak extraction
            if peaks:
                repcols.extend(c for c in ['MS2 Time', 'Peak Time', 'Peak Intensity',
                                           'Peak Width (sec)', 'Peak Comment']
                               if c not in repcols)
                if _peaks_options['peak_area'] and 'Peak Area' not in repcols:
                    repcols.append('Peak Area')
                if _peaks_options['reporter_ions']:
                    repcols.extend(c for c in ['Rep114', 'Rep115', 'Rep116', 'Rep117']
                                   if c not in repcols)           

            repcols.insert(0, 'File')

            report_file = os.path.join(chosen_folder, combined_file)

            if os.path.exists(report_file):
                os.remove(report_file)

            report = mzReport.writer(report_file, columns=repcols)
            isMZD = isinstance(report, mzReport.mzDB.SQLiteWriter)

            mascot_headers = []
        else:
            report_files = []

        if dates:
            mid_d = zip(mascot_ids, dates, [None]*len(mascot_ids))
        elif local_dat_files:
            mid_d = zip(["Local File"]*len(local_dat_files), [None]*len(local_dat_files), local_dat_files)
        else:
            mid_d = [(mid,None, None) for mid in mascot_ids]

        for mascot_id,date,local in mid_d:
            mascot_id = str(mascot_id)

            if ':' in mascot_id:
                (mascot_id, date) = mascot_id.split(':', 1)

            mascot_id = str(mascot_id).zfill(6)

            if not (date or local):
                date = self.mascot.get_date(mascot_id)

            logger_message(30, 'Generating Multiplierz-Mascot Report for JobID %s...' % mascot_id)

            if ext == '.mzid':
                logger_message(30, 'Downloading mzIdentML File...')
                destination = chosen_folder if chosen_folder else myData

                reportfilename = "F%s.mzid" % mascot_id
                outputfile = os.path.join(destination, reportfilename)
                report_file = self.mascot.download_mzid(mascot_id,
                                                        save_file = outputfile,
                                                        date = date)
                assert report_file == outputfile

                report_files.append(report_file)
                continue
                # mzIdentML files don't use the rest of this function;
                # what they contain is essentially fixed, to multiplierz.        

            if not local:
                logger_message(30, 'Downloading Mascot DAT File...')
                dat_file = self.mascot.download_dat(chosen_folder, mascot_id, date)
            else:
                dat_file = os.path.abspath(local)
                mascot_id = os.path.basename(local).split('.')[0]

            if dat_file:
                logger_message(20, 'Mascot DAT File Downloaded!')
                mascot_dat_file = interface.MascotDatFile(dat_file, **_mascot_options)

                if percolatorDirectory and mascot_dat_file.hasDecoyHits():
                    print "Running Mascot Percolator..."
                    mascot_dat_file.close()
                    percolatedDatFile = runPercolator(dat_file, percolatorDirectory)
                    mascot_dat_file = interface.MascotDatFile(dat_file, **mascot_options)

                    if self.cleanup:
                        os.remove(dat_file)
                    dat_file = percolatedDatFile


            else:
                logger_message(40, 'Failed to download DAT file for %s' % mascot_id)
                continue




            if self.mascot.version != mascot_dat_file.res_file.getMascotVer()[:len(self.mascot.version)]:
                print ("Mascot version mismatch detected; changing version from %s to %s" 
                       % (self.mascot.version, mascot_dat_file.res_file.getMascotVer()[:len(self.mascot.version)]))
                self.mascot.version = mascot_dat_file.res_file.getMascotVer()[:len(self.mascot.version)]

            if not combined_file:
                # Report column stuff moved from above, in order to handle version dependency.  (Heavy sigh.)
                # figuring out the report columns. start with defaults...
                repcols = mzReport.default_columns[:]

                # mascot 2.3 can have multiple databases so add a column for that
                if self.mascot.version >= '2.3':
                    repcols.insert(1, 'Protein Database')

                # these are the columns coming out of the dat file, need them separate
                res_cols = repcols[:]

                # add columns for peak extraction
                if peaks:
                    repcols.extend(c for c in ['MS2 Time', 'Peak Time', 'Peak Intensity',
                                               'Peak Width (sec)', 'Peak Comment']
                                   if c not in repcols)
                    if _peaks_options['peak_area'] and 'Peak Area' not in repcols:
                        repcols.append('Peak Area')
                    if _peaks_options['reporter_ions']:
                        repcols.extend(c for c in ['Rep114', 'Rep115', 'Rep116', 'Rep117']
                                       if c not in repcols)

                if mascot_prot_cov:
                    repcols.append('Protein Coverage')            


            #Get MS File Name
            mascot_header = mascot_dat_file.mascot_header()

            ms_file_name = mascot_header[7][1] or ('F%s' % mascot_id)


            if not combined_file:
                filename = os.path.basename(ms_file_name)
                if mascotIDInResultName and filename.endswith('.mgf'):
                    filename = filename[:-4] + "." + mascot_id

                report_file = os.path.join(chosen_folder, filename + ext)

                if os.path.exists(report_file):
                    os.remove(report_file)

                report = mzReport.writer(report_file, columns=repcols)
                isMZD = isinstance(report, mzReport.mzDB.SQLiteWriter)

                     
            if mascot_web and (_mascot_web_options['ms2_img'] or _mascot_web_options['mascot_var_mods']):
                gen_options = {}
                try:
                    gen_options['ms2_img'] = _mascot_web_options['ms2_img']
                except KeyError:
                    pass
                try:
                    gen_options['mascot_var_mods'] = _mascot_web_options['mascot_var_mods']
                except KeyError:
                    pass

                mascot_web_gen = self.mascot_web(mascot_id, date=date,
                                                 dat_file=(dat_file if _mascot_web_options['mascot_ms2'] else None),
                                                 isMZD=isMZD,
                                                 **gen_options)
                mascot_web_gen.next()

            if mascot_prot_cov:
                prot_cov_gen = self.mascot_prot_coverage(mascot_id,
                                                         _mascot_options['ion_cutoff'],
                                                         date)
                prot_cov_gen.next()

            prot_desc_dict = {}

            if self.mascot.version != mascot_dat_file.res_file.getMascotVer()[:len(self.mascot.version)]:
                raise TypeError, "Incorrect version of Mascot selected. %s %s" % (self.mascot.version, mascot_dat_file.res_file.getMascotVer()[:len(self.mascot.version)])

            missing_desc_count = 0
            for row in mascot_dat_file.peptide_report():
                row = mzReport.ReportEntry(res_cols, row)

                if rank_one and row['Peptide Rank'] != 1:
                    continue

                if (not local) and not (row['Protein Description'] or row['Protein Mass']):
                    if row['Accession Number'] not in prot_desc_dict:
                        missing_desc_count += 1
                        # Very slow!
                        #prot_desc_dict[row['Accession Number']] = self.mascot.get_description(row['Accession Number'],
                                                                                              #row.get('protein database', '1').split('::')[0],
                                                                                              #mascot_id,
                                                                                              #date)
                    row['Protein Description'],row['Protein Mass'] = prot_desc_dict.get(row['Accession Number'], ('-', '-'))

                md = []
                #if peaks:
                    #(new_row, img_tuples) = peak_gen.send(row)
                    #row.update(new_row)
                    #md.extend(img_tuples)

                if mascot_web and (_mascot_web_options['ms2_img']
                                   or _mascot_web_options['mascot_var_mods']):
                    (vartext, img_tup) = mascot_web_gen.send(row)
                    if _mascot_web_options['mascot_var_mods']:
                        row['Variable Modifications'] = vartext
                    if _mascot_web_options['ms2_img']:
                        md.append(img_tup)

                if mascot_prot_cov:
                    (prot_cov, md_tup) = prot_cov_gen.send(row)
                    row['Protein Coverage'] = prot_cov
                    md.append(md_tup)

                if combined_file:
                    row['File'] = ms_file_name

                report.write(row, metadata=md)
            if missing_desc_count:
                print "Missing protein info for %d PSMs." % missing_desc_count
                
            if peaks:
                peak_gen.close()
            if mascot_web and (_mascot_web_options['ms2_img']
                               or _mascot_web_options['mascot_var_mods']):
                mascot_web_gen.close()
            if mascot_prot_cov:
                prot_cov_gen.close()

            # Mascot-decoy-data finder!
            if mascot_dat_file.hasDecoyHits():
                decoy_dat_file = interface.MascotDatFile(dat_file, decoyMode = True, **_mascot_options)
                for row in decoy_dat_file.peptide_report():
                    report.write(row)
                decoy_dat_file.close()

            if not combined_file:
                if os.path.splitext(report_file)[1].lower() in ('.xls', '.xlsx', '.mzd'):
                    report.close()

                    self.mascot_headers(report_file, [(None, mascot_header)])
                    if protein_report:
                        self.prot_report(report_file, mascot_dat_file.protein_report())
                else:
                    report.close()

                report_files.append(report_file)
            else:
                if os.path.splitext(report_file)[1].lower() in ('.xls', '.xlsx', '.mzd'):
                    mascot_headers.append((ms_file_name, mascot_header))

            mascot_dat_file.close()
            if self.cleanup and not local_dat_files:
                os.remove(dat_file)

            logger_message(30, 'Multiplierz-Mascot Report for JobID %s Generated!' % mascot_id)

        if combined_file:
            if os.path.splitext(report_file)[1].lower() in ('.xls', '.xlsx', '.mzd'):
                report.close()

                self.mascot_headers(report_file, mascot_headers)
                # not supported right now: protein reports for XLS.
                if isMZD and protein_report:
                    self.prot_report(report_file, None)
            else:
                report.close()
           
        return [report_file] if combined_file else report_files
    
    def retrieve_report_data(self, mascot_id, report_columns,
                             date = None,
                             rank_one = False,
                             ext = '.xlsx',
                             local_dat_file = None,
                             retain_dat_file = False,
                             **mascot_options):
        """
        Retrieves a .DAT file (or opens one local) and reads
        into tabular format, handling several issues such as
        decoy mode data, header data, and version-dependent
        rows.
        """
        
        if local_dat_file:
            datfile = os.path.abspath(local_dat_file)
        else:
            datfile = self.mascot.download_dat(myData, mascot_id, date)
            
        dat_interface = interface.MascotDatFile(datfile, **mascot_options)
        
        header_sheet = dat_interface.mascot_header()
        data_sheet = []
        for values in dat_interface.peptide_report():
            row = mzReport.ReportEntry(mzReport.default_columns,
                                       values)
            
            if rank_one and row['Peptide Rank'] != 1:
                continue
            
            data_sheet.append(row)
        dat_interface.close()        
        
        # With 2.5 Mascot, the decoy mode test crashes Python outright.
        #if dat_interface.hasDecoyHits():
            #decoy_dat_interface = interface.MascotDatFile(dat_file, 
                                                          #decoyMode = True,
                                                          #**_mascot_options)
            #for row in decoy_dat_interface.peptide_report():
                #row = mzReport.ReportEntry(report_columns, values)
                #data_sheet.append(row)
            #decoy_dat_interface.close()        
        
        if not retain_dat_file:
            os.remove(datfile)
        
        return header_sheet, data_sheet
        
    
    def get_reports(self, mascot_ids,
                    dates = None,
                    outputfile = None,
                    ext = None,
                    chosen_folder = '',
                    **report_kwargs):
        

        
        report_columns = mzReport.default_columns
        if float(self.mascot.version) >= 2.3 and 'Protein Database' not in report_columns:
            report_columns.insert(1, 'Protein Database')        
            
        if dates:
            assert len(dates) == len(mascot_ids), "Mismatched date list provided."
            mascot_searches = zip(mascot_ids, dates)
        else:
            mascot_searches = [(x, None) for x in mascot_ids]
            
        reports = []
        for mascot_id, date in mascot_searches:
            header, psms = self.retrieve_report_data(mascot_id, report_columns,
                                                     date,
                                                     **report_kwargs)
            datafilename = header[7][1] or mascot_id
            reports.append((mascot_id, datafilename, header, psms))
            
        if outputfile and ext:
            if not outputfile.lower().endswith(ext):
                outputfile += '.' + ext        
        elif not outputfile:
            if not ext:
                ext = 'xlsx'
            outputfile = '_'.join(mascot_ids) + '.' + ext.strip('.')            
            #outputfile = '.'.join([reports[0][1], ext.strip('.')])
        elif outputfile and not ext:
            ext = outputfile.split('.')[-1]
            assert ext in ['csv', 'xlsx', 'xls', 'mzd', 'mzid']

          
        assert outputfile 
        if chosen_folder and not os.path.isabs(outputfile):
            outputfile = os.path.join(chosen_folder, outputfile)
        
        if ext and 'mzid' in ext:
            assert len(mascot_ids) == 1, ("Combined result file not supported "
                                          "for mzIdentML files.")
            self.mascot.download_mzid(mascot_id,
                                      save_file = outputfile,
                                      date = date)
        elif len(mascot_ids) == 1:
            mascot_id, datafilename, header, psms = reports[0]
            if not outputfile:
                outputfile = datafilename + '.' + ext
            output = mzReport.writer(outputfile,
                                     columns = header[0],
                                     sheet_name = 'Mascot Header')
            for line in header[1:]:
                output.write(line)
            output.close()
            
            output = mzReport.writer(outputfile, columns = report_columns,
                                     sheet_name = 'Data')
            for psm in psms:
                output.write(psm)
            output.close()
        
        else:
            #report_columns.insert(0, 'File')
            if not outputfile:
                raise IOError("Combined report file name must be specified.")
            
            if (outputfile.lower().endswith('xls') or
                outputfile.lower().endswith('xlsx') or
                outputfile.lower().endswith('mzd')):            
                for m_id, datafilepath, header, _ in reports:
                    datafilename = os.path.basename(datafilepath)
                    output = mzReport.writer(outputfile, columns = header[0],
                                    sheet_name = '%s Mascot Header' % datafilename)
                    for line in header[1:]:
                        output.write(line)
                    output.close()
            else:
                extension = outputfile.split('.')[-1]
                print "Omitting header tables due to %s format." % extension
            
            output = mzReport.writer(outputfile, columns = ['File'] + report_columns,
                                     sheet_name = 'Data')
            for _, datafilepath, _, psms in reports:
                datafilename = os.path.basename(datafilepath)
                for psm in psms:
                    psm['File'] = datafilename
                    output.write(psm)
            
            output.close()


        return outputfile




def runPercolator(datfile, mascotPercolatorPath):
    """
    Going by the known-successful-run:

    """

    assert os.path.exists(mascotPercolatorPath), "Cannot find Mascot Percolator at %s" % mascotPercolatorPath

    if not mascotPercolatorPath.lower().endswith('MascotPercolator.jar'):
        mascotPercolatorPath = os.path.join(mascotPercolatorPath, 'MascotPercolator.jar')

    resultTemplate = datfile + '.percolation'
    expectedDatFile = resultTemplate + '.datp'
    result = subprocess.call(['java', '-cp', mascotPercolatorPath, 'cli.MascotPercolator',
                              '-target', datfile, '-decoy', datfile, '-newDat', '-out',
                              resultTemplate])

    assert not result
    assert os.path.exists(expectedDatFile)

    return expectedDatFile



