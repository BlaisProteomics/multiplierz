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

# this module is for parsing an mzIdentML file. Mascot-only for the moment.

# This used to prefer lxml-based ElementTree; there *shouldn't* be differences?
import xml.etree.cElementTree as ET
from collections import defaultdict

# namespace dictionary for the unimod schema. 'f' is just to distinguish
# the 'lower' function, using the Blais URL is not really necessary.
NS = { 'mzid': 'http://psidev.info/psi/pi/mzIdentML/1.0',
       'f': 'http://blais.dfci.harvard.edu' }

def lower(context, a):
    return a[0].lower()


class mzIdentML(object):
    '''Class to represent an mzIdentML (.mzid) file. There are methods
    for iterating through the file in various ways, and for extracting
    other information (such as search methods).'''

    def __init__(self, file_name):
        self.file_name = file_name

        # this is a special class for evaluating XPath expression quickly
        # we have to give it the extra function 'lower' because lxml doesn't
        # have xpath 2.0 functions
        self.mzid = ET.XPathEvaluator(ET.parse(self.file_name),
                                      namespaces=NS,
                                      extensions={(NS['f'], 'lower'): lower})

        self._proteins = None

    def __iter__(self):
        #protein_report = []

        prot_hits = self.mzid('./mzid:DataCollection/mzid:AnalysisData/mzid:ProteinDetectionList'
                              '/mzid:ProteinAmbiguityGroup/mzid:ProteinDetectionHypothesis')

        prot_hits.sort(key=lambda p: float(p.xpath(r'./mzid:cvParam[@name="mascot:score"]',
                                                   namespaces=NS)[0].get('value')),
                       reverse=True)

        for prot_rank,prot in enumerate(prot_hits):
            db_seq = prot.get('DBSequence_ref')

            prot_acc,prot_desc = self.proteins[db_seq]

            prot_score = prot.xpath(r'./mzid:cvParam[@name="mascot:score"]',
                                    namespaces=NS)

            prot_score = float(prot_score[0].get('value')) if prot_score else None

            peptides = prot.xpath(r'./mzid:cvParam[@name="distinct peptide sequences"]',
                                  namespaces=NS)

            peptides = int(peptides[0].get('value')) if peptides else None

            pep_evs = prot.xpath(r'./mzid:PeptideHypothesis/@PeptideEvidence_Ref',
                                 namespaces=NS)

            #prot_matches = len(pep_evs)

            prot_matches = sum([1 for si in [self.mzid('./mzid:DataCollection/mzid:AnalysisData'
                                                       '/mzid:SpectrumIdentificationList'
                                                       '/mzid:SpectrumIdentificationResult'
                                                       '/mzid:SpectrumIdentificationItem'
                                                       '[child::mzid:PeptideEvidence[@id="%s"]]' % pe)
                                             for pe in pep_evs]
                                if si and si[0].get('passThreshold') == 'true'])


            for pe in pep_evs:
                pep_evidence = self.mzid('./mzid:DataCollection/mzid:AnalysisData'
                                         '/mzid:SpectrumIdentificationList'
                                         '/mzid:SpectrumIdentificationResult'
                                         '/mzid:SpectrumIdentificationItem'
                                         '/mzid:PeptideEvidence[@id="%s"]' % pe)

                if not pep_evidence:
                    raise ValueError('This mzIdentML file is missing a PeptideEvidence node with id: %s' % pe)
                else:
                    pep_evidence = pep_evidence[0]

                spec = pep_evidence.getparent()

                spec_title = spec.xpath('parent::*/mzid:cvParam[@name="spectrum title"]',
                                        namespaces=NS)
                spec_title = spec_title[0].get('value') if spec_title else None

                pep_id = spec.get('Peptide_ref')
                pep_node = self.mzid('./mzid:SequenceCollection/mzid:Peptide[@id="%s"]' % pep_id)[0]

                pep_seq = self._peptide(pep_node)
                pep_var_mod = self._modifications(pep_node)

                #prot_acc,prot_desc = self.proteins[pep_evidence[0].get('DBSequence_Ref')]
                pep_start = pep_evidence.get('start', None)
                pep_end = pep_evidence.get('end', None)
                pep_res_before = pep_evidence.get('pre', None)
                pep_res_after = pep_evidence.get('post', None)
                pep_miss = pep_evidence.get('missedCleavages', None)

                pep_exp_mz = float(spec.get('experimentalMassToCharge'))
                pep_exp_z = int(spec.get('chargeState'))

                pep_score = spec.xpath('./mzid:cvParam[@name="mascot:score"]',
                                       namespaces=NS)
                pep_score = float(pep_score[0].get('value')) if pep_score else None

                #print pep_score, spec.get('passThreshold')

                pep_rank = int(spec.get('rank'))

                if spec.get('passThreshold') == 'true':
                    yield [prot_rank + 1, #prot_hit_num,
                           prot_acc,
                           prot_desc,
                           None, # prot_mass,
                           prot_matches, # this is too high...counting bad peptides
                           prot_score,
                           pep_seq,
                           pep_var_mod,
                           pep_exp_mz,
                           pep_exp_z,
                           None, #pep_calc_mr,
                           None, #pep_delta,
                           pep_score,
                           pep_rank,
                           pep_start,
                           pep_end,
                           pep_res_before,
                           pep_res_after,
                           pep_miss,
                           spec_title,
                           None]

    @property
    def proteins(self):
        if self._proteins is None:
            self._proteins = {}

            prot_list = self.mzid(r'./mzid:SequenceCollection/mzid:DBSequence[@id][@accession]')

            for prot in prot_list:
                prot_id = prot.get('id')
                prot_acc = prot.get('accession')

                prot_desc = prot.xpath(r'./mzid:cvParam[@name="protein description"]',
                                       namespaces=NS)
                prot_desc = prot_desc[0].get('value') if prot_desc else None

                self._proteins[prot_id] = (prot_acc, prot_desc)

        return self._proteins

    def _peptide(self, pep_node, mods=False):
        '''Returns a peptide (with or without modifications) from a peptide node'''

        peptide = [''] + list(pep_node.xpath('./mzid:peptideSequence/text()', namespaces=NS)[0]) + ['']

        ## Should I replace unknown residues with the matched amino acid?
        #rep_aa = pep_node.xpath('./mzid:SubstitutionModification[@location]', namespaces=NS)
        #for rep in rep_aa:
            #loc = int(rep.get('location'))
            #if loc == 0 or loc == len(peptide) - 1:
                #raise ValueError('Substitution modification on a terminus?')

            #peptide[loc] = rep.get('replacementResidue')

        if mods:
            mod_list = pep_node.xpath('./mzid:Modification', namespaces=NS)
            for mod in mod_list:
                mod_pos = int(mod.get('location'))
                u = mod.xpath('./mzid:cvParam[@cvRef="UNIMOD"]', namespaces=NS)
                if u:
                    mod_name = u[0].get('name')
                else:
                    mod_name = mod.get('monoisotopicMassDelta')

                if mod_pos == 0:
                    peptide[0] = '[%s]-' % mod_name
                elif mod_pos == len(peptide) - 1:
                    peptide[-1] = '-[%s]' % mod_name
                else:
                    peptide[mod_pos] = '[%s]%s' % (mod_name, peptide[mod_pos])

        return ''.join(peptide)

    def _modifications(self, pep_node):
        '''Returns a modification string from a peptide node'''

        peptide = pep_node.xpath('./mzid:peptideSequence/text()', namespaces=NS)[0]

        mods = []

        mod_list = pep_node.xpath('./mzid:Modification', namespaces=NS)
        for mod in mod_list:
            mod_pos = int(mod.get('location'))
            u = mod.xpath('./mzid:cvParam[@cvRef="UNIMOD"]', namespaces=NS)
            if u:
                mod_name = u[0].get('name')
            else:
                mod_name = 'Unknown Mod (%s)' % mod.get('monoisotopicMassDelta')

            if mod_pos == 0:
                mods.append('N-term: %s' % mod_name)
            elif mod_pos == len(peptide) + 1:
                mods.append('C-term: %s' % mod_name)
            else:
                mods.append('%s%d: %s' % (peptide[mod_pos-1], mod_pos, mod_name))

        return '; '.join(mods)
