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

"""Formats OMSSA CSV into multiplierz format

"""

import csv
import os

import multiplierz.mzReport
import multiplierz.mass_biochem

class OMSSA_CSV():
    def __init__(self, file):
        self.orig_file = file

    def format(self, new_file_name=None):
        fin = open(self.orig_file)
        csvReader = csv.reader(fin)

        headers = [x.strip() for x in next(csvReader)]
        new_headers = self.convert_headers(headers)
        spec_desc = new_headers.index('Spectrum Description')
        mz        = new_headers.index('Experimental mz')
        charge    = new_headers.index('Charge')
        seq       = new_headers.index('Peptide Sequence')
        var_mods  = new_headers.index('Variable Modifications')
        new_data = []
        for prerow in csvReader:
            row = [x.strip() for x in prerow]
            row[spec_desc] = self.convert_spectrum(row[spec_desc])
            row[mz] = self.convert_mass(row[mz],row[charge])
            (the_seq,the_mods) = self.convert_seq(row[seq])
            row[seq] = the_seq
            row[var_mods] = the_mods
            new_data.append( row )
        fin.close()


        dir_split = os.path.split(self.orig_file)
        if not new_file_name:
            new_file_name = os.path.join(dir_split[0], "mz_" + dir_split[1])

        report = multiplierz.mzReport.writer(new_file_name, columns = new_headers)
        for data in new_data:
            report.write(data)

        report.close()

    def convert_headers(self, headers):
        replacements = [
            ('Mass', 'Experimental mz'),
            ('Theo Mass', 'Predicted mr'),
            ('Accessions', 'Accession Number'),
            ('Peptide', 'Peptide Sequence'),
            ('E-value', 'Peptide Score'),
            ('Charge', 'Charge'),
            ('Start', 'Start Position'),
            ('Stop', 'End Position'),
            ('Mods', 'Variable Modifications'),
            ('Filename/id', 'Spectrum Description'),
            ('Spectrum number', 'Query'),
        ]
        replacements = dict(replacements)

        new_headers = []
        for header in headers:
            if header in replacements:
                new_headers.append(replacements[header])
            else:
                new_headers.append(header)

        return new_headers

    def convert_mass(self,mw,charge):
        mz = (float(mw) + multiplierz.mass_biochem.AW['H']*float(charge))/float(charge)
        return repr(mz)

    def convert_seq(self,seq):
        new_seq = seq[:]
        the_mods = ""
        while new_seq.find("s")> -1 :
            offset = new_seq.find("s")
            new_seq = new_seq[:offset] + "S" + new_seq[(offset+1) :]
            the_mods += "S%d: Phospho;" % (offset+1)
        while new_seq.find("t")> -1 :
            offset = new_seq.find("t")
            new_seq = new_seq[:offset] + "T" + new_seq[(offset+1) :]
            the_mods += "T%d: Phospho;" % (offset+1)
        while new_seq.find("y")> -1 :
            offset = new_seq.find("y")
            new_seq = new_seq[:offset] + "Y" + new_seq[(offset+1) :]
            the_mods += "Y%d: Phospho;" % (offset+1)
        while new_seq.find("m")> -1 :
            offset = new_seq.find("m")
            new_seq = new_seq[:offset] + "M" + new_seq[(offset+1) :]
            the_mods += "M%d: Oxid;" % (offset+1)
        return (new_seq,the_mods)

    def convert_spectrum(self, spectrum):
        spec = spectrum[ (spectrum.rfind("\\")+1) : ]
        return spec
