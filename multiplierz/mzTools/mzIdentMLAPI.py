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
            print targetLookup
            print "No match for {0}".format(oldfilename)

    mzid.write(mzidFile)



def parseSpectrumTitle(spectrumTitle):
    dotsplit = spectrumTitle.split(".")
    if len(dotsplit) == 4 and dotsplit[1] == dotsplit[2]:
        # Then this is a .RAW file title.
        return dotsplit[0] + ".raw", int(dotsplit[1])
    else:
        raise NotImplementedError, "Unrecognized spectrum title type!"

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
            print ("Multiplierz mzIdentML reader was written against the specification "
                   "for mzIdentML version 1.1.0; this file is version %s, so there may "
                   "be discrepancies in the results." % self.root.get("version"))

        self.pfx = self.root.tag[:-9] # Take out "mzIdentML"        

        #print "Indexing file..." # Really this could all be done on an as-needed basis.
        self.parentMap = dict(sum([[(c.get("id"),p.get("id")) for c in p if c.get("id") and p.get("id")]
                                   for p in ch(self.tree.getiterator(self.pfx + "ProteinAmbiguityGroup"),
                                               self.tree.getiterator(self.pfx + "SpectrumIdentificationResult"))], []))

        proteinElements = self.root.getiterator(self.pfx + "ProteinDetectionHypothesis")
        #self.proteinLookup = {prot.get("id"):prot for prot in proteinElements}
        self.proteinLookup = dict([(prot.get("id"), prot) for prot in proteinElements])

        spectrumElements = self.root.getiterator(self.pfx + "SpectrumIdentificationResult")
        #self.spectrumLookup = {spec.get("id"):spec for spec in spectrumElements}
        self.spectrumLookup = dict([(spec.get("id"), spec) for spec in spectrumElements])

        pepIdElements = self.root.getiterator(self.pfx + "SpectrumIdentificationItem")
        #self.pepIdLookup = {pepid.get("id"):pepid for pepid in pepIdElements}
        self.pepIdLookup = dict([(pepid.get("id"), pepid) for pepid in pepIdElements])

        evidenceElements = self.root.getiterator(self.pfx + "PeptideEvidence")
        #self.evidenceLookup = {evi.get("id"):evi for evi in evidenceElements}
        self.evidenceLookup = dict([(evi.get("id"), evi) for evi in evidenceElements])

        dbElements = self.root.getiterator(self.pfx + "DBSequence")
        #self.dbLookup = {db.get("id"):db for db in dbElements}
        self.dbLookup = dict([(db.get("id"), db) for db in dbElements])

        peptideElements = self.root.getiterator(self.pfx + "Peptide")
        #self.peptideLookup = {pep.get("id"):pep for pep in peptideElements}
        self.peptideLookup = dict([(pep.get("id"), pep) for pep in peptideElements])

        fileElements = self.root.getiterator(self.pfx + "SpectraData")
        #self.fileLookup = {data.get("id"):data for data in fileElements}
        self.fileLookup = dict([(data.get("id"), data) for data in fileElements])
        #print "Indexing complete."

        self.dataFileScans = {}
        self.filePointers = {}

    def numberOfProteins():
        return len(self.proteinLookup)
    def numberOfAmbiguityGroups():
        return len(self.root.findall(self.pfx + "ProteinAmbiguityGroups"))
    def numberOfSpectra():
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
            spectrum = max(spectra, key = lambda x: self.giveCVs(x)["Mascot:score"])

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
        spectrum_score = spectCVs['Mascot:score']
        spectrum_expect = spectCVs['Mascot:expectation value']

        resultCVs = self.giveCVs(resultEl)
        try:
            spectDesc = resultCVs["spectrum title"]
        except KeyError:
            spectDesc = resultEl.get("spectrumID")


        spectrum_info = {"Calculated mz" : float(spectEl.get("calculatedMassToCharge")),
                         #"Calculated PI" : spectEl.get("calculatedPI"),
                         "Charge" : int(spectEl.get("chargeState")),
                         "Experimental mz" : float(spectEl.get("experimentalMassToCharge")),
                         "Passed Threshold" : spectEl.get("passThreshold") == 'true',
                         "Rank" : int(spectEl.get("rank")),
                         #"Spectrum Name" : spectEl.get("name"),
                         "Spectrum ID" : spectEl.get("id"),
                         "Spectrum Description" : spectDesc,
                         #"MS2 Time" : retentionTime,
                         "File" : dataEl.get("location").replace(r'file:///', ''),
                        }
        # How is this not a function of the peptide assignment?
        spectrum_info['Delta'] = (float(spectrum_info['Experimental mz']) -
                                  float(spectrum_info['Calculated mz']))
        #spectrumInfo.update(pepInfo)
        ## All peptide data will be the same in all evidence elements?
        ## No, this is designed to accomodate multiple peptide matches!
        #evidenceRef = spectEl.find(self.pfx + "PeptideEvidenceRef")
        #pepEvidence = self.evidenceLookup[evidenceRef.get("peptideEvidence_ref")]

        #pepInfo = self.peptideInfo(pepEvidence.get("peptide_ref"))

        #accessions = []
        #scores = []
        #for evidenceRef in spectEl.findall(self.pfx + "PeptideEvidenceRef"):
            #evidence = self.evidenceLookup[evidenceRef.get("peptideEvidence_ref")]
            #dbSeq = self.dbLookup[evidence.get("dBSequence_ref")]
            #accessions.append(dbSeq.get("accession"))

            #try:
                #cvparams = self.giveCVs(evidence)
                #score = [x['value'] for x in cvparams if x['name'] == 'Mascot:score'][0]
                #scores.append(score)
            #except IndexError:
                #pass
        #accessions = "; ".join(accessions)
        #scores = "; ".join(scores)

        psms = []
        for evidenceRef in spectEl.findall(self.pfx + "PeptideEvidenceRef"):
            evidence = self.evidenceLookup[evidenceRef.get("peptideEvidence_ref")]
            dbSeq = self.dbLookup[evidence.get("dBSequence_ref")]
            accessions = dbSeq.get("accession")
            pepInfo = self.peptideInfo(evidence.get("peptide_ref"))
            
            try:
                cvparams = self.giveCVs(evidence)
                score = float([x['value'] for x in cvparams if x['name'] == 'Mascot:score'][0])
                expect = float([x['value'] for x in cvparams if x['name'] == 'Mascot:expectation value'][0])
            except IndexError:
                score = spectrum_score
                expect = spectrum_expect
            
            psm = spectrum_info.copy()
            psm.update(pepInfo)
            psm['Accession Number'] = accessions
            psm['Peptide Score'] = float(score)
            psm['Expectation Value'] = float(expect)
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

        modificationData = []
        for modification in pepEl.findall(self.pfx + "Modification"):
            cvs = modification.findall(self.pfx + "cvParam")
            try:
                kind = [x.get("name") for x in cvs if x.get("cvRef") == 'UNIMOD'][0]
            except IndexError:
                kind = "NOT-UNIMOD"
            modificationData.append((modification.get("location"), kind, 
                                     modification.get("residues")))
        for modification in pepEl.findall(self.pfx + "SubstitutionModification"):
            modificationData.append((modification.get("location"),
                                     modification.get("replacementResidue") + "-replacement",
                                     modification.get("originalResidue")))

        peptideInfo = {"Peptide Sequence" : pepEl.find(self.pfx + "PeptideSequence").text,
                       "Variable Modifications" : renderModificationString(modificationData),
                       "Peptide ID" : pepEl.get("id"),
                       #"Peptide Name" : pepEl.get("name")
                       }

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





    def peptideSummary(self):
        """
        Info on the top-ranking peptide match for each spectrum in the file."

        Written into a report file, this produces an output roughly equivalent
        to the standard Mascot XLS result file.
        """

        spectraData = []
        for spectrumResult in self.spectrumLookup.values():
            items = spectrumResult.findall(self.pfx + "SpectrumIdentificationItem")
            topRank = [s for s in items if int(s.get("rank")) == 1]
            for spectrumItem in topRank:
                #info = self.spectrumInfo(spectrumItem.get("id"))
                #info['Delta'] = str(float(info['Experimental mz']) - float(info['Calculated mz']))
                #spectraData.append(info)
                spectraData += self.spectrumInfo(spectrumItem.get("id"))

        return spectraData

    def proteinSummary(self, reportAll = False):
        """
        Info on proteins found in the file; if "reportAll" is true, all proteins
        indicated are retrieved, otherwise some effort is made to take the best
        protein of each ambiguity group.
        """

        if not reportAll:
            assert self.giveCVs(self.proteinLookup.values()[0])["Mascot:score"], \
                   "Protein ranking currently only supported via mascot scoring."

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
                #protInfo = self.proteinInfo(protein.get("id"))
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
                #datafile = dataEl.get("location")
                datafile = derivedData
            try:
                data = self.filePointers[datafile]
            except KeyError:
                data = mzFile(datafile)
                self.filePointers[datafile] = data
            #rT = float(self.giveCVs(spectrumResult)["MS:1001114"]) / 60.0

            #scanName = spectrumResult.get("spectrumID") # Perhaps?  Not entirely clear.
            #try:
                #scanNum = int(scanName)
            #except ValueError:
                #scanNum = int(scanName.split("=")[1])



            #scan = data.cscan(data.scan_time_from_scan_name(scanNum))

            #for spectrumItem in [x for x in spectrumResult 
                                    #if x.tag == (self.pfx + 'SpectrumIdentificationItem')]:
            for spectrumItem in spectrumResult.getiterator(self.pfx + 'SpectrumIdentificationItem'):
                #mz = float(spectrumItem.get("experimentalMassToCharge"))
                #scanHeader = min([x for x in data.scan_info(rT - 0.1, rT + 0.1, mz - 1, mz + 1)
                                    #if x[3] == 'MS2'],
                                    #key = lambda x: abs(x[1] - mz))
                #scan = data.scan(scanHeader[0], centroid = True)
                scan = data.scan(scanNum, centroid = True)

                if len(scan) > 500:
                    scan = sorted(scan, key = lambda x: x[1], reverse = True)[:500]                

                #try:
                    #scans = self.dataFileScans[datafile]
                #except KeyError:
                    #scans = data.scan_info()
                    #self.dataFileScans[datafile] = scans
                try:
                    #fragmentation = [x for x in spectrumItem
                                        #if x.tag == (self.pfx + "Fragmentation")][0]
                    fragmentation = spectrumItem.getiterator(self.pfx + 'Fragmentation').next()
                except StopIteration:
                    fragmentation = xml.SubElement(spectrumItem, self.pfx + "Fragmentation")
                    #fragmentation = [x for x in spectrumItem
                                        #if x.tag == (self.pfx + "Fragmentation")][0]
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


        softwareUsed = self.root.getiterator(self.pfx + "AnalysisSoftwareList").next()
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
        fragmentation = spectrum.getiterator(self.pfx + "Fragmentation").next()

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
            #for fragArray in [x for x in iontype if x.tag == "FragmentArray"]:
            for fragArray in iontype.getiterator(self.pfx + "FragmentArray"):
                if fragArray.get("measure_ref") == "m_mz":
                    mzArray = fragArray.get("values").split()
                elif fragArray.get("measure_ref") == "m_intensity":
                    intArray = fragArray.get("values").split()
            assert mzArray and intArray, "Incomplete annotation!"

            mzArray = [float(x) for x in mzArray]
            intArray = [float(x) for x in intArray]
            annotationArray = zip(mzArray, intArray)
            break

        return annotationArray


    #def fullReport(self):
        #"""An attempt to replicate the content of basic Mascot output."""

        #protsum = self.proteinSummary()
        #pepsum = self.peptideSummary()

        #protLookup = {} # Using this to not lose ambiguity group info?
        #for protein in protsum:
            #protLookup[protein['id']] = protein            

        #results = []
        #for peptide in pepsum:
            #protIds = peptide['Accession Number']

            #proteins = []
            #for protId in protIds.split('; '):
                #proteins.append(protLookup['PDH_' + protId + '_0'])
            #protein = max(proteins, key = lambda x: float(x['Rank']))

            #try:
                #del protein['id']
            #except KeyError:
                #pass # How does this happen???
            #peptide.update(protein)
            #results.append(peptide)

        #return results


    def close(self):
        self.mzid.close()




















