import os, sys
from subprocess import call

def convert_report(reportfile, outputfile = None):
    # Reads .mzid file to generate a Mascot-like report which has most
    # important fields; doesn't report protein-level statistics, missed cleavages,
    # or Mascot-specific fields (e.g., query.)
    
    from multiplierz.mzTools.mzIdentMLAPI import mzIdentML
    from multiplierz.mass_biochem import remove_protons
    from multiplierz.mzReport import reader, writer
    
    if not outputfile:
        outputfile = reportfile + '.csv'
    
    report = mzIdentML(reportfile)
    out = writer(outputfile,
                 columns = ['File', 'Rank', 'Accession Number', 'Protein Description',
                            'Peptide Sequence', 'Variable Modifications',
                            'Experimental mz', 'Predicted mr', 'Charge', 'Delta',
                            'Peptide Score', 'Expectation Value',
                            'Start Position', 'End Position',
                            'Preceding Residue', 'Following Residue', 
                            'Spectrum Description'])
    for row in report.peptideSummary():
        row['Predicted mr'] = remove_protons(row['Calculated mz'], int(row['Charge']))
        del row['Calculated mz']
        del row['Passed Threshold']
        del row['Peptide ID']
        del row['Spectrum ID']
        out.write(row)
    out.close()
        


class MSGFPlusSearch(dict):
    # Really basic for now, just takes a parameter file and applies it.
    def __init__(self, parameterfile = None, database = None, program_path = None):
        if not all([parameterfile, database, msgf_path]):
            raise NotImplementedError
        
        self.parfile = parameterfile
        self.fasta = database
        self.msgf_path = msgf_path
    
    def run_search(self, datafile, output = None,
                   extraction_parameters = {},
                   convert_output = False):
        if not datafile.lower().endswith('.mgf'):
            datafile = extract(datafile, **extraction_parameters)
        
        if not output:
            output = datafile + '.mzid'
        elif output and not output.lower().endswith('mzid'):
            output += '.mzid'
        
        
        cmd = ["java",
               "-jar", self.msgf_path,
               "-conf", self.parfile,
               "-o", output,
               "-s", datafile]
        if self.fasta:
            cmd += ["-d", self.fasta]
        
        retcode = call(cmd)
        assert not retcode, "MSGF+ return value was nonzero: %s" % retcode
        assert os.path.exists(output), "Output file was not found."
        
        if convert_output:
            conv_output = output + '.tsv'
            tsv_file = call(["java", "-cp", self.msgf_path,
                             "edu.ucsd.msjava.ui.MzIDToTsv",
                             "-i", output,
                             "-o", conv_output])
            assert os.path.exists(conv_output)
            return conv_output
        else:
            return output




if __name__ == '__main__':
    convert_report(r'\\rc-data1\blaise\ms_data_share\Max\MSGF_test\2019-04-07-Xcal1.raw.hcd.mzid',
                         r'\\rc-data1\blaise\ms_data_share\Max\MSGF_test\2019-04-07-Xcal1.raw.hcd.mzid.TEST.csv')    