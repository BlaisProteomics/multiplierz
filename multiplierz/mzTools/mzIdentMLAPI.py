import xml.etree.cElementTree as xml
import os
from itertools import chain as ch

from multiplierz.mzAPI import mzFile
from multiplierz import __version__



__all__ = ['mzIdentML']

def deRef(elem, tag):
    childs = elem.findall(tag)
    assert len(childs) == 1
    return childs[0]

def unzip(thing): return [list(t) for t in zip(*thing)]




def transplantSpectraData(mzidFile, newTargets):
    mzid = xml.parse(mzidFile)
    root = mzid.getroot()

    def title(filename):
        return os.path.basename('.'.join(filename.split(".")[0]))

    targetLookup = dict([(title(target), target) for target in newTargets])

    prefix = root.tag[:-9]
    for fileElement in root.iter(prefix + "SpectraData"):
        oldfilename = title(fileElement.get("location"))
        try:
            newfilepath = targetLookup[oldfilename]
            fileElement.set("location", os.path.normpath(newfilepath))
        except KeyError:
            print(targetLookup)
            print(("No match for {0}".format(oldfilename)))

    mzid.write(mzidFile)



def parseSpectrumTitle(spectrumTitle):
    dotsplit = spectrumTitle.split(".")
    if len(dotsplit) == 4 and dotsplit[1] == dotsplit[2]:
        # Then this is a .RAW file title.
        return dotsplit[0] + ".raw", int(dotsplit[1])
    else:
        raise NotImplementedError("Unrecognized spectrum title type!")

def renderModificationString(modList):
    modStrs = []
    for pos, name, res in modList:
        if res == None:
            if (pos == '0') or not pos:
                posStr = "N-Term: "
            else:
                posStr = "C-Term: "
        else:
            posStr = res + pos + ": "

        modStrs.append(posStr + name)
    return "; ".join(modStrs)



class mzIdentML(object):
    def __init__(self, filename):
        self.filename = filename
        self.mzid = open(filename, "r")
        self.tree = xml.parse(self.mzid)
        self.root = self.tree.getroot()

        assert "mzIdentML" in self.root.tag, "%s does not appear to be a valid mzIdentML file."
        if self.root.get("version") != "1.1.0":
            print(("Multiplierz mzIdentML reader was written against the specification "
                   "for mzIdentML version 1.1.0; this file is version %s, so there may "
                   "be discrepancies in the results." % self.root.get("version")))

        self.getSourceMode()

        self.pfx = self.root.tag[:-9] # Take out "mzIdentML", when present.
        # (There's probably a more proper way to do this!)

        # The below be done as-needed
        # or at least in a single iterator (but it seems fast enough.)
        self.parentMap = dict(sum([[(c.get("id"),p.get("id")) for c in p if c.get("id") and p.get("id")]
                                   for p in ch(self.tree.getiterator(self.pfx + "ProteinAmbiguityGroup"),
                                               self.tree.getiterator(self.pfx + "SpectrumIdentificationResult"))], []))

        proteinElements = self.root.getiterator(self.pfx + "ProteinDetectionHypothesis")
        self.proteinLookup = {prot.get("id"):prot for prot in proteinElements}

        spectrumElements = self.root.getiterator(self.pfx + "SpectrumIdentificationResult")
        self.spectrumLookup = {spec.get("id"):spec for spec in spectrumElements}        

        pepIdElements = self.root.getiterator(self.pfx + "SpectrumIdentificationItem")
        self.pepIdLookup = {pepid.get("id"):pepid for pepid in pepIdElements}

        evidenceElements = self.root.getiterator(self.pfx + "PeptideEvidence")
        self.evidenceLookup = {evi.get("id"):evi for evi in evidenceElements}

        dbElements = self.root.getiterator(self.pfx + "DBSequence")
        self.dbLookup = {db.get("id"):db for db in dbElements}

        peptideElements = self.root.getiterator(self.pfx + "Peptide")
        self.peptideLookup = {pep.get("id"):pep for pep in peptideElements}

        fileElements = self.root.getiterator(self.pfx + "SpectraData")
        self.fileLookup = {data.get("id"):data for data in fileElements}

        self.dataFileScans = {}
        self.filePointers = {}

    def getSourceMode(self):
        analysissoft_list = next((x for x in self.root if 'AnalysisSoftwareList' in x.tag), [])
        self.mode = None
        for software in analysissoft_list:
            soft_name = software.attrib['name']
            if 'Mascot' in soft_name:
                self.mode = "mascot"
                break
            elif 'MS-GF+' in soft_name:
                self.mode = 'msgf+'
                break

    def numberOfProteins(self):
        return len(self.proteinLookup)
    def numberOfAmbiguityGroups(self):
        return len(self.root.findall(self.pfx + "ProteinAmbiguityGroups"))
    def numberOfSpectra(self):
        return len(self.spectrumLookup)

    def giveCVs(self, element):
        cvs = element.findall(self.pfx + "cvParam")
        cvData = {}
        for cv in cvs:
            vals = cv.attrib
            try:
                value = vals["value"]
            except KeyError:
                value = None
            cvData[vals["accession"]] = value
            cvData[vals["name"]] = value

        return cvData    

    def proteinInfo(self, proteinID):
        pEl = self.proteinLookup[proteinID]
        pepHypEls = pEl.findall(self.pfx + "PeptideHypothesis")
        elements = []
        for pepHyp in pepHypEls:
            spectra = []
            for spectrumItemRef in pepHyp.findall(self.pfx + "SpectrumIdentificationItemRef"):
                specItemRefref = spectrumItemRef.get("spectrumIdentificationItem_ref")
                specItem = self.pepIdLookup[specItemRefref]
                spectra.append(specItem)

            # Only for mascot files!  What are equivalents from other search engines?
            if mode == 'mascot':
                spectrum = max(spectra, key = lambda x: float(self.giveCVs(x)["Mascot:score"]))
            elif mode == 'msgf+':
                spectrum = min(spectra, key = lambda x: float(self.giveCVs(x)["MS-GF:SpecEValue"]))
            else:
                raise NotImplementedError("Only implemented for Mascot and MSGF+ reports.")

            pepEvidence = self.evidenceLookup[pepHyp.get("peptideEvidence_ref")]
            peptide = self.peptideLookup[pepEvidence.get("peptide_ref")]
            dbSeq = self.dbLookup[pepEvidence.get("dBSequence_ref")]

            elements.append((spectrum, pepEvidence, peptide, dbSeq))

        dbEntry = unzip(elements)[3][0]    
        assert all([dbEntry == x for x in unzip(elements)[3]])

        ambiguity = self.parentMap[proteinID]


        protData = {"name" : pEl.get("name"),
                    "id" : pEl.get("id"),
                    "Spectrum IDs" : [pep.get("id") for pep in unzip(elements)[0]],
                    "Peptides" : [pep.get("id") for pep in unzip(elements)[2]],
                    "Sequence" : (dbEntry.find(self.pfx + "Seq").text if 
                                  dbEntry.find(self.pfx + "Seq") else None),
                    "Accession" : dbEntry.get("accession"),
                    "database" : dbEntry.get("searchDatabase_ref"),
                    "Ambiguity Group" : ambiguity
                    }

        return protData


    def spectrumInfo(self, spectrumId):
        """
        Returns the available information for the specified scan.
        """
        spectEl = self.pepIdLookup[spectrumId]
        assert spectEl.tag == self.pfx + "SpectrumIdentificationItem"
        resultEl = self.spectrumLookup[self.parentMap[spectrumId]]
        dataEl = self.fileLookup[resultEl.get("spectraData_ref")]

        spectCVs = self.giveCVs(spectEl)
        # Annoying that these fields are Mascot specific.
        if self.mode == 'mascot':
            try:
                spectrum_score = spectCVs['Mascot:score']
            except KeyError:
                # Mascot spectral library search result.
                spectrum_score = spectCVs['MSPepSearch:score']
            spectrum_expect = spectCVs['Mascot:expectation value']
        elif self.mode == 'msgf+':
            spectrum_expect = spectCVs['MS-GF:SpecEValue']
            spectrum_score = spectCVs['MS-GF:RawScore']
        else:
            raise NotImplementedError("Not an MSGF or Mascot report.")

        resultCVs = self.giveCVs(resultEl)
        try:
            spectDesc = resultCVs["spectrum title"]
        except KeyError:
            spectDesc = resultEl.get("spectrumID")


        spectrum_info = {"Calculated mz" : float(spectEl.get("calculatedMassToCharge")),
                         "Charge" : int(spectEl.get("chargeState")),
                         "Experimental mz" : float(spectEl.get("experimentalMassToCharge")),
                         "Passed Threshold" : spectEl.get("passThreshold") == 'true',
                         "Rank" : int(spectEl.get("rank")),
                         "Spectrum ID" : spectEl.get("id"),
                         "Spectrum Description" : spectDesc,
                         "File" : dataEl.get("location").replace(r'file:///', ''),
                         "Peptide Score":spectrum_score,
                         "Expectation Value":spectrum_expect
                        }
        # How is this not a function of the peptide assignment?
        spectrum_info['Delta'] = (float(spectrum_info['Experimental mz']) -
                                  float(spectrum_info['Calculated mz']))

        psms = []
        for evidenceRef in spectEl.findall(self.pfx + "PeptideEvidenceRef"):
            evidence = self.evidenceLookup[evidenceRef.get("peptideEvidence_ref")]
            dbSeq = self.dbLookup[evidence.get("dBSequence_ref")]
            accessions = dbSeq.get("accession")
            try:
                desc = self.giveCVs(dbSeq)['protein description']
            except KeyError:
                desc = 'No Description'
            pepInfo = self.peptideInfo(evidence.get("peptide_ref"))
            
            start, stop = evidence.get('start'), evidence.get('stop')
            pre, post = evidence.get('pre'), evidence.get('post')
            
            psm = spectrum_info.copy()
            psm.update(pepInfo)
            psm['Accession Number'] = accessions
            psm['Protein Description'] = desc
            psm['Preceding Residue'] = pre
            psm['Following Residue'] = post
            psm['Start Position'] = start
            psm['End Position'] = stop
            psms.append(psm)
            
        return psms


    def evidenceInfo(self, evidenceId):
        eviEl = evidenceLookup[evidenceId]

        dbSeq = dbLookup[eviEl.get("dBSequence_ref")]

        evidenceInfo = {"Protein Accession" : dbSeq.get("accession"),
                        "Database" : dbSeq.get("searchDatabase_ref")
                        } # Leaving out data about peptide location within protein.

        evidenceInfo.update(petideInfo(eviEl.get("peptide_ref")))
        return evidenceInfo


    def peptideInfo(self, peptideId):
        pepEl = self.peptideLookup[peptideId]
        pep_seq = pepEl.find(self.pfx + "PeptideSequence").text
        
        modificationData = []
        for modification in pepEl.findall(self.pfx + "Modification"):
            cvs = modification.findall(self.pfx + "cvParam")
            kind = next((x.get("name") for x in cvs if x.get("cvRef") == 'UNIMOD'), "NOT-UNIMOD")
            modificationData.append((modification.get("location"), kind, 
                                     modification.get("residues")))
        for modification in pepEl.findall(self.pfx + "SubstitutionModification"):
            modificationData.append((modification.get("location"),
                                     modification.get("replacementResidue") + "-replacement",
                                     modification.get("originalResidue")))

        # MS-GF+ (and possibly others?) don't always supply residue location in
        # mo
        for i in range(len(modificationData)):
            if modificationData[i][2] is None and modificationData[i][0] is not None:
                pos = int(modificationData[i][0])
                residue = pep_seq[pos-1] # 1-indexed position.
                modificationData[i] = str(pos), modificationData[i][1], residue
            
        peptideInfo = {"Peptide Sequence" : pep_seq,
                       "Variable Modifications" : renderModificationString(modificationData),
                       "Peptide ID" : pepEl.get("id")}

        return peptideInfo


    def componentInfo(pepID, protID):
        """
        Info on the relation between a protein and one of its component found peptides.
        Returns None if peptide is not part of specified protein.
        """

        protEl = self.proteinLookup(protID)
        pepEv = None
        for evidence in protEl.findall(self.prefix + "PeptideHypothesis"):
            eviID = evidence.get("peptideEvidence_ref")
            evidence = self.evidenceLookup[eviID]
            if evidence.get("peptide_ref") == pepID:
                pepEv = evidence        
        if pepEv == None:
            return None

        componentInfo = {"Start Position" : pepEv.get("start"),
                         "End Position" : pepEv.get("end"),
                         "Preceding Residue" : pepEv.get("pre"),
                         "Following Residue" : pepEv.get("post"),
                         }

        return componentInfo


    def peptideSummary(self, top_rank_only = True):
        """
        Info on the top-ranking peptide match for each spectrum in the file."

        Written into a report file, this produces an output roughly equivalent
        to the standard Mascot XLS result file.
        """

        spectraData = []
        for spectrumResult in list(self.spectrumLookup.values()):
            items = spectrumResult.findall(self.pfx + "SpectrumIdentificationItem")
            if top_rank_only:
                items = [s for s in items if int(s.get("rank")) == 1]
            for spectrumItem in items:
                spectraData += self.spectrumInfo(spectrumItem.get("id"))

        return spectraData

    def proteinSummary(self, reportAll = False):
        """
        Info on proteins found in the file; if "reportAll" is true, all proteins
        indicated are retrieved, otherwise some effort is made to take the best
        protein of each ambiguity group.
        """

        assert self.mode == 'mascot', "Protein summary only implemented for Mascot reports."

        proteinData = []
        for ambiguity in self.root.getiterator(self.pfx + "ProteinAmbiguityGroup"):
            bestProtein = None
            bestRating = None
            ambiguityLength = len(ambiguity.findall(self.pfx + "ProteinDetectionHypothesis"))
            ambiguityId = ambiguity.get("id")
            for protein in ambiguity.findall(self.pfx + "ProteinDetectionHypothesis"):
                if reportAll:
                    protInfo = self.proteinInfo(protein.get("id"))
                    protInfo.update({"Protein Redundancy":ambiguityLength,
                                     "Protein Group ID":ambiguityId})
                    proteinData.append(protInfo)
                else:
                    rating = self.giveCVs(protein)["Mascot:score"]
                    if rating > bestRating:
                        bestRating = rating
                        bestProtein = protein

            if not reportAll:
                protInfo = self.proteinInfo(bestProtein.get("id"))
                protInfo.update({"Protein Redundancy":ambiguityLength,
                                 "Protein Group ID":ambiguityId})
                proteinData.append(protInfo)

        return proteinData


    def writeIonAnnotations(self, datafile = None, in_place = False):
        for spectrumList in self.root.getiterator(self.pfx + "SpectrumIdentificationList"):
            try:
                fragtab = [x for x in spectrumList if x.tag == self.pfx + "FragmentationTable"][0]
            except IndexError:
                fragtab = xml.SubElement(spectrumList, "FragmentationTable")

            intMeasure = xml.SubElement(fragtab, "Measure")
            intMeasure.set("id", "m_intensity")
            intMeasureKind = xml.SubElement(intMeasure, "cvParam")
            intMeasureKind.set("cvRef", "PSI-MS")
            intMeasureKind.set("accession", "MS:1001226")
            intMeasureKind.set("name", "product ion intensity")

            mzMeasure = xml.SubElement(fragtab, "Measure")
            mzMeasure.set("id", "m_mz")
            mzMeasureKind = xml.SubElement(mzMeasure, "cvParam")
            mzMeasureKind.set("cvRef", "PSI-MS")
            mzMeasureKind.set("accession", "MS:1001225")
            mzMeasureKind.set("name", "product ion m/z")

        for spectrumResult in self.root.getiterator(self.pfx + "SpectrumIdentificationResult"):
            dataEl = self.fileLookup[spectrumResult.get("spectraData_ref")]


            spectrumTitle = self.giveCVs(spectrumResult)['spectrum title']
            derivedData, scanNum = parseSpectrumTitle(spectrumTitle)            
            if not datafile:
                datafile = derivedData
            try:
                data = self.filePointers[datafile]
            except KeyError:
                data = mzFile(datafile)
                self.filePointers[datafile] = data

            for spectrumItem in spectrumResult.getiterator(self.pfx + 'SpectrumIdentificationItem'):
                scan = data.scan(scanNum, centroid = True)

                if len(scan) > 500:
                    scan = sorted(scan, key = lambda x: x[1], reverse = True)[:500]                


                try:
                    fragmentation = next(spectrumItem.getiterator(self.pfx + 'Fragmentation'))
                except StopIteration:
                    fragmentation = xml.SubElement(spectrumItem, self.pfx + "Fragmentation")
                iontype = xml.SubElement(fragmentation, self.pfx + "IonType")
                iontype.set("index", "0 " * len(scan))
                iontype.set("charge", "0")

                ionKind = xml.SubElement(iontype, self.pfx + "cvParam")
                ionKind.set("cvRef", "PSI-MS")
                ionKind.set("accession", "MS:1001240")
                ionKind.set("name", "non-identified ion")

                def listStr(thing):
                    out = ""
                    for x in thing:
                        out += (str(x) + " ")
                    return out

                mzArray = xml.SubElement(iontype, self.pfx + "FragmentArray")
                mzArray.set("values", listStr(unzip(scan)[0]))
                mzArray.set("measure_ref", "m_mz")

                intArray = xml.SubElement(iontype, self.pfx + "FragmentArray")
                intArray.set("values", listStr(unzip(scan)[1]))
                intArray.set("measure_ref", "m_intensity")

        if not in_place:
            outputFile = self.filename[:-5] + "_annotated.mzid"
        else:
            outputFile = self.filename


        softwareUsed = next(self.root.getiterator(self.pfx + "AnalysisSoftwareList"))
        mzDesktopEl = xml.SubElement(softwareUsed, self.pfx + "AnalysisSoftware")
        mzDesktopEl.set("id", "DFCI Multiplierz v1.1.0")
        mzDesktopEl.set("name", "Multiplierz")
        mzDesktopEl.set("uri", "http://sourceforge.net/projects/multiplierz/")
        mzDesktopEl.set("version", __version__)
        softwareName = xml.SubElement(mzDesktopEl, "SoftwareName")
        nameParam = xml.SubElement(softwareName, "userParam")
        nameParam.set("name", "Multiplierz")

        self.mzid.close()

        output = open(outputFile, "w")
        self.tree.write(output)
        output.close()

        self.mzid = open(outputFile, "r")

    def isAnnotated(self):
        # This will have to be updated if multiplierz gains the ability
        # to write mzIdentML directly aside from annotation.
        software = self.root.getiterator(self.pfx + "AnalysisSoftware")
        return any([el.get("name") == "Multiplierz" for el in software])

    def getAnnotation(self, spectrumItemId):
        spectrum = self.pepIdLookup[spectrumItemId]
        fragmentation = next(spectrum.getiterator(self.pfx + "Fragmentation"))

        annotationArray = []
        for iontype in fragmentation:
            try:
                self.giveCVs(iontype)['non-identified ion']
            except KeyError:
                pass
            if not all([x == 0 for x in iontype.get("index").split()]):
                pass

            mzArray = None
            intArray = None
            for fragArray in iontype.getiterator(self.pfx + "FragmentArray"):
                if fragArray.get("measure_ref") == "m_mz":
                    mzArray = fragArray.get("values").split()
                elif fragArray.get("measure_ref") == "m_intensity":
                    intArray = fragArray.get("values").split()
            assert mzArray and intArray, "Incomplete annotation!"

            mzArray = [float(x) for x in mzArray]
            intArray = [float(x) for x in intArray]
            annotationArray = list(zip(mzArray, intArray))
            break

        return annotationArray

    def close(self):
        self.mzid.close()
