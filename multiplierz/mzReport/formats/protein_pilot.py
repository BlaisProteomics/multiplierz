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


"""Formats Protein Pilot spreadsheet into multiplierz format

"""

import re
import os

from collections import defaultdict

import multiplierz.mzReport

from multiplierz import logger_message


class ProteinPilot():
    reps = {'N': 'Protein Rank',
            'Accessions': 'Accession Number',
            'Names': 'Protein Description',
            'Total': 'Protein Score',
            'Sequence': 'Peptide Sequence',
            'Modifications': 'Variable Modifications',
            #'Prec m/z': 'mz',
            'Prec m/z': 'Experimental mz',
            'Theor z': 'Charge',
            'Theor MW': 'Predicted mr',
            'dMass': 'Delta',
            'Conf': 'Peptide Score',
            'Cleavages': 'Missed Cleavages',
            'Spectrum': 'Spectrum Description',
            'Time': 'MS2 Time',
            '%Cov': 'Protein Coverage'}

    def __init__(self, file_name):
        self.orig_file = file_name

    def format(self, new_file_name=None):
        fh = open(self.orig_file)

        #Change header names
        header_line = fh.readline()
        headers = header_line.strip().split('\t')

        new_headers = multiplierz.mzReport.default_columns[:]
        new_headers.extend(('MS2 Time', 'Protein Coverage'))
        new_headers.remove('Peptide Rank')
        new_headers.remove('Query')
        new_headers.extend(h for h in headers if h not in self.reps)
        if 'Unused' in new_headers:
            new_headers.remove('Unused')
        if 'Contrib' in new_headers:
            new_headers.remove('Contrib')
        if 'Sc' in new_headers:
            new_headers.remove('Sc')

        rows = []

        protein_matches = defaultdict(int)

        for line in fh:
            new_data = dict((h,None) for h in new_headers)
            data = dict(list(zip(headers,line[:-1].split('\t'))))

            protein_matches[data['Accessions']] += 1

            for h in data:
                if h in self.reps and self.reps[h] in new_data:
                    new_data[self.reps[h]] = data[h]
                elif h in new_data:
                    new_data[h] = data[h]
                elif h not in ('Unused','Contrib','Sc'):
                    logger_message(10, 'Missing key: %s' % h)

            if 'Modifications' in data:
                new_data['Variable Modifications'] = self.convert_var_mod(data['Modifications'])

            if 'Cleavages' in data:
                new_data['Missed Cleavages'] = data['Cleavages'].count('missed')

            if 'Time' in data and 'Spectrum' in data:
                new_data['Spectrum Description'] =  self.convert_spectrum(data['Spectrum'], data['Time'])

            rows.append(new_data)

        fh.close()

        for row in rows:
            row['Protein Matches'] = protein_matches[row['Accession Number']]

        dir_split = os.path.split(self.orig_file)
        if not new_file_name:
            new_file_name = os.path.join(dir_split[0], "mz_" + dir_split[1])

        report = multiplierz.mzReport.writer(new_file_name, columns=new_headers)

        for row in rows:
            report.write(row)

        report.close()

    def convert_var_mod(self, mods):
        mods = mods.split(";")
        new_var_mod = []
        mod_re = re.compile("(.+)\(([A-Z])\)\@(\d+)")
        term_re = re.compile("(.+)\@([NC]\-term)")
        for mod in mods:
            mod = re.sub("\s+", "", mod)
            m1 = mod_re.match(mod)
            m2 = term_re.match(mod)
            if m1:
                new_var_mod.append('%s%s: %s' % (m1.group(2), m1.group(3), m1.group(1)))
            elif m2:
                new_var_mod.append('%s: %s' % (m2.group(2), m2.group(1)))

        new_var_mod = '; '.join(new_var_mod)

        return new_var_mod

    def convert_spectrum(self, spectrum, time):
        data = spectrum.split('.')
        file_name = data[0]
        cycle = int(data[-2])
        experiment = int(data[-1])
        sample = int(data[1])
        time = float(time)
        spec = ("File: %s.wiff, Sample: (sample number %d),"
                " Elution: %0.3f min, Period: 1, Cycle(s):"
                " %d (Experiment %d)") % (file_name, sample, time, cycle, experiment)

        return spec

if __name__ == "__main__":
    import mzGUI
    pilot = ProteinPilot(mzGUI.file_chooser())
    pilot.format()



