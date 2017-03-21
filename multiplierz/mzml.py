from multiplierz.mzAPI import mzFile
from multiplierz.internalAlgorithms import gzOptOpen
import xml.etree.ElementTree as xml
import xml.parsers.expat as expat
import tempfile, sqlite3
import struct
import zlib
import sqlite3 as sqlite
import cPickle as pickle
import multiprocessing
import base64
import sys, os

namespace = '{http://psi.hupo.org/ms/mzml}'
def ns(tag):
    assert namespace
    return namespace + tag
def unns(tag):
    return tag[len(namespace):]

# Compression level is a tradeoff between speed and file size;
# 2 seems to get the benefits of both.
def marshal(data):
    return base64.b64encode(zlib.compress(pickle.dumps(data), 2))
def demarshal(data):
    return pickle.loads(zlib.decompress(base64.b64decode(data)))

def child(element, tag):
    return [x for x in element if x.tag == ns(tag)][0]

def cvp(parent, cvRef = 'MS', unitCvRef = 'MS',
        **parameters):
    return xml.SubElement(parent, ns('cvParam'), cvRef = cvRef, unitCvRef = unitCvRef,
                          **parameters)

def uncvp(el):
    assert el.tag == ns('cvParam')
    return {'name':el.get('name'),
            'units':el.get('unitName'),
            'value':el.get('value')}


    
def encodeBinaryFloats(floats, compression = False):
    output = ''
    for flt in floats:
        output += struct.pack('d', flt)
        
    if compression:
        output = zlib.compress(output)
    hexput = base64.b64encode(output)
    return hexput

def decodeBinaryFloats(binary, bit64 = True, compression = True):
    hexary = base64.b64decode(binary)
    if compression:
        hexary = zlib.decompress(hexary)    
        
    iterinc = 8 if bit64 else 4
    structtype = 'd' if bit64 else 'f'
    output = []
    for i in range(0, len(hexary), iterinc):
        try:
            output.append(struct.unpack(structtype, hexary[i:i+iterinc])[0])
        except struct.error as err:
            if len(hexary) - i < iterinc:
                break
            else:
                raise err
    return output

    

    
def makeSpectrumXML(specData):
    specEl = xml.Element(ns('spectrum'))
    specEl.set('defaultArrayLength', len(specData['spectrum']))
    specEl.set('id', specData['Spectrum Description'])
    specEl.set('index', specData['index'])
    specEl.set('sourceFileRef', specData['Source'])
    
    #cvp(specEl, accession="MS:1000130", name = 'positive scan')
    cvp(specEl, accession="MS:1000789", name = 'enhanced multiply charged spectrum') # I guess?
    if specData['centroid']:
        cvp(specEl, accession="MS:000127", name = 'centroid spectrum')
    else:
        cvp(specEl, accession="MS:1000128", name = 'profile spectrum')
    
    if specData['MS Level'] == 'MS1':
        cvp(specEl, accession='MS:1000579', name = 'MS1 spectrum')
    elif specData['MS Level'] == 'MS2':
        cvp(specEl, accession='MS:1000580', name = 'MSn spectrum')
    else:
        raise NotImplementedError, specData['MS Level']
    

    # I'm assuming there's only one scan.
    # No idea what multiple scans per spectrum would even mean.
    scanlistdataEl = xml.SubElement(specEl, ns('scanList'),
                                    count = 1)
    scanEl = xml.SubElement(scanlistdataEl, ns('scan'),
                            sourceFileRef = specData['Source'],
                            spectrumRef = specData['Spectrum Description'])

    windowList = xml.SubElement(scanEl, ns('scanWindowList'), count = 2)
    cvp(windowList, accession="MS:1000501", name="scan window lower limit",
        unitAccession="MS:1000040", unitName="m/z", value=specData['mz range'][0])
    cvp(windowList, accession="MS:1000500", name="scan window upper limit",
        unitAccession="MS:1000040", unitName="m/z",  value=specData['mz range'][1])
    
    if specData['MS Level'] != 'MS1':
    # Not include this for MS1 scans?
        precList = xml.SubElement(scanEl, ns('precursorList'), count = 1)
        precEl = xml.SubElement(precList, ns('precursor'), 
                                sourceFileRef = specData['Source'],
                                spectrumRef = specData['previous MS1'])
        isoWinEl = xml.SubElement(precEl, ns('isolationWindow'))
        cvp(isoWinEl, accession="MS:1000827", name="isolation window target m/z",
            unitAccession="MS:1000040", unitName="m/z", value = specData['precursor'])
        # Dunno how to get isolation width from most formats?
    
        selIonList = xml.SubElement(precEl, ns('selectedIonList'), count = 1)
        selIonEl = xml.SubElement(selIonList, ns('selectedIon'), count = 1)
        cvp(selIonEl, accession="MS:1000744", name="selected ion m/z",
            unitAccession="MS:1000040", unitName="m/z", value = specData['precursor'])
        # This 'precursor' value is the ion itself, whereas the previous is the
        # isolation window; different?  Where to derive the difference?
        if specData['charge']:
            cvp(selIonEl, accession='MS:1000041', name = 'charge state',
                value = specData['charge'])
    
        # How to deal with MS1s here?
        activationEl = xml.SubElement(precEl, ns('precursor'))
        if specData['dissociation mode'] == 'cid':
            cvp(activationEl, accession="MS:1000133", name="collision-induced dissociation")
        elif specData['dissociation mode'] == 'hcd':
            cvp(activationEl, accession="MS:1000422", name="high-energy collision-induced dissociation")
        elif specData['dissociation mode'] == 'etd':
            cvp(activationEl, accession="MS:1000250", name="electron capture dissociation")
        else:
            raise NotImplementedError, "Can't encode dissociation mode %s" % specData['dissociation mode']
        
    
    encodedMZs = encodeBinaryFloats(x[0] for x in specData['spectrum'])
    encodedInts = encodeBinaryFloats(x[1] for x in specData['spectrum'])
    
    binaryDataArrayList = xml.SubElement(specEl, ns('binaryDataArrayList'), count = 2)
    mzDataArray = xml.SubElement(binaryDataArrayList, ns('binaryDataArray'), 
                                 arrayLength = len(specData['spectrum']),
                                 encodedLength = len(encodedMZs))
    cvp(mzDataArray, accession='MS:1000574', name='zlib compression')
    cvp(mzDataArray, accession='MS:1000514', name='m/z array')
    mzBinary = xml.SubElement(mzDataArray, ns('binary'))
    mzBinary.text = encodedMZs
    intDataArray = xml.SubElement(binaryDataArrayList, ns('binaryDataArray'),
                                  arrayLength = len(specData['spectrum']),
                                  encodedLength = len(encodedInts))
    cvp(intDataArray, accession='MS:1000574', name='zlib compression')
    cvp(intDataArray, accession='MS:1000515', name='intensity array')
    intBinary = xml.SubElement(intDataArray, ns('binary'))
    intBinary.text = encodedInts
    
    return specEl


def readSpectrumXML(spectrumEl):
    cvparams = [uncvp(x) for x in list(spectrumEl) if x.tag == ns('cvParam')]
    cvparams = dict([(x['name'], x) for x in cvparams])

    specdata = {}
    specdata['index'] = spectrumEl.attrib['index']
    specdata['Spectrum Description'] = spectrumEl.get('id')
    specdata['Source'] = spectrumEl.get('sourceFileRef')
    specdata['centroid'] = 'profile spectrum' not in cvparams
    try:
        specdata['mz range'] = cvparams['lowest observed m/z']['value'], cvparams['highest observed m/z']['value']
    except KeyError:
        specdata['mz range'] = None, None

    scan = child(child(spectrumEl, 'scanList'), 'scan')
    scancvps = dict([(uncvp(x)['name'], uncvp(x)) for x in list(scan) if x.tag == ns('cvParam')])
    specdata['time'] = scancvps['scan start time']['value']
    if 'filter string' in scancvps:
        specdata['filter'] = scancvps['filter string']['value']    

    if int(cvparams['ms level']['value']) == 2:
        specdata['precursor'] = cvparams['base peak m/z']['value']
        specdata['charge'] = 0 # How does one get the precursor charge??
    else:
        specdata['precursor'] = 0

    arrays = child(spectrumEl, 'binaryDataArrayList')
    for array in arrays:
        assert array.tag == ns('binaryDataArray')
        arraycvps = [uncvp(x)['name'] for x in list(array) if x.tag == ns('cvParam')]
        decompress = 'no compression' not in arraycvps
        doublefloats = '32-bit float' not in arraycvps

        binary = child(array, 'binary')
        if binary.text:
            spectrum = decodeBinaryFloats(binary.text, doublefloats, decompress)
        else:
            spectrum = []
        
        if 'm/z array' in arraycvps:
            mzspectrum = spectrum
        elif 'intensity array' in arraycvps:
            intspectrum = spectrum
        else:
            raise Exception, 'Unidentified array type: %s' % arraycvps

    specdata['spectrum'] = zip(mzspectrum, intspectrum)


    return specdata

def readChromatoXML(chromatoEl):
    cvparams = [uncvp(x) for x in list(chromatoEl) if x.tag == ns('cvParam')]
    cvparams = dict([(x['name'], x) for x in cvparams])

    if 'total ion current chromatogram' in cvparams:
        span = 'total'
    else:
        raise NotImplementedError, "Non-total XIC element."

    for array in list(child(chromatoEl, 'binaryDataArrayList')):
        assert array.tag == ns('binaryDataArray')
        arraycvps = [uncvp(x)['name'] for x in list(array) if x.tag == ns('cvParam')]
        decompress = 'no compression' not in arraycvps
        doublefloats = '32-bit float' not in arraycvps

        binary = child(array, 'binary')
        if 'time array' in arraycvps:
            timepts = decodeBinaryFloats(binary.text, doublefloats, compression = decompress)
        elif 'intensity array' in arraycvps:
            intpts = decodeBinaryFloats(binary.text, doublefloats, compression = decompress)
        else:
            raise Exception, 'Unidentified array type: %s' % arraycvps

    # This may want to be a more inclusive return type.
    return span, zip(timepts, intpts)


            
        
def readSpectrumData(data):
    datafile = data.data_file
    scaninfo = data.scan_info(0, 999999)
    filters = data.filters()
    
    spectrumDataObjs = []
    prev_MS1 = None
    for (rt, mz, scanNum, level, mode), (_, filter)  in zip(scaninfo, filters):
        filterWords = filter.split()
        lowerBound, upperBound = [float(x) for x in filterWords[-1][1:-1].split('-')]
        
        specData = {}
        specData['index'] = scanNum
        specData['Spectrum Description'] = filter
        specData['Source'] = datafile
        specData['centroid'] = mode == 'c'
        specData['precursor'] = mz
        specData['mz range'] = lowerBound, upperBound
        
        try:
            mz2, z = data.scanPrecursor(scanNum)
            specData['charge'] = int(z)
        except AttributeError:
            specData['charge'] = None
        
        specData['MS Level'] = level
        if level == 'MS1':
            prev_MS1 = filter 
        else:
            specData['previous MS1'] = prev_MS1
    
            if '@cid' in filter:
                specData['dissociation mode'] = 'cid'
            elif '@hcd' in filter:
                specData['dissociation mode'] = 'hcd'
            elif '@etd' in filter:
                specData['dissociation mode'] = 'etd'
            else:
                raise NotImplementedError, "Can't extract dissociation mode from '%s'." % filter
        
        specData['spectrum'] = data.scan(scanNum)
        
        
        spectrumDataObjs.append(specData)
        
    
    return spectrumDataObjs
    
    



    
def readXML(xmlfile):
    tree = xml.parse(open(xmlfile))
    root = tree.getroot()
    if root.tag == ns('indexedmzML'):
        root = child(root, 'mzML')
    run = child(root, 'run')

    spectrumList = child(run, 'spectrumList')
    chromatoList = child(run, 'chromatogramList')

    spectrumData = []
    chromatoData = []
    for spectrumEl in spectrumList:
        if not spectrumEl.tag == ns('spectrum'):
            continue
        spectrumData.append(readSpectrumXML(spectrumEl))
    for chromatoEl in chromatoList:
        if not chromatoEl.tag == ns('chromatogram'):
            continue
        chromatoData.append(readChromatoXML(chromatoEl))

    additionalData = {} # There's probably more to be found in the top-level.

    return additionalData, spectrumData, chromatoData



    
    


        
def mzmlToSqlite_writer(sqlitefile, inputs):            
    connection = sqlite.connect(sqlitefile)
    cursor = connection.cursor()
    
    createSpectrumTable = "CREATE TABLE spectra(ind int, mz real, rt int, data blob)"
    cursor.execute(createSpectrumTable)
    createChromatoTable = "CREATE TABLE chromato(ind int, startmz real, stopmz real, startrt real, stoprt real, data blob)"
    cursor.execute(createChromatoTable)
    connection.commit()
    
    i = 0
    while True:
        task, data = inputs.get()
        if task == 'spectrum':
            pdata = marshal(data)
            index, prec, time = data['index'], data['precursor'], data['time']
            cursor.execute('INSERT INTO spectra VALUES (%s,%s,%s,"%s")' % (index, prec, time, pdata))
        elif task == 'chromatogram':
            pdata = marshal(data)
            if data[0] == 'total':
                cursor.execute('INSERT INTO chromato VALUES (0,0,0,0,0,"%s")' % pdata)
            else:
                cursor.execute('INSERT INTO chromato VALUES (%s,%s,%s,%s,%s,"%s")'
                               % span + (pdata,))
        elif task == 'stop':
            break
        else:
            raise NotImplementedError, task
    
        if i > 100:
            connection.commit()
            i = 0
        else:
            i += 1
            
    connection.commit()
    connection.close()
        
    
    
def mzmlToSqlite(xmlfile, sqlitefile):
    parser = xml.iterparse(xmlfile)
    
    writeQueue = multiprocessing.Queue()
    writerProc = multiprocessing.Process(target = mzmlToSqlite_writer,
                                         args = (sqlitefile, writeQueue))
    writerProc.start()
    
    
    for evt, obj in parser:
        if obj.tag == ns('spectrum'):
            writeQueue.put(('spectrum', readSpectrumXML(obj)))
            obj.clear()
        elif obj.tag == ns('chromatogram'):
            writeQueue.put(('chromatogram', readChromatoXML(obj)))
            obj.clear()
        
    writeQueue.put(('stop', None))
            
    writerProc.join()
    return sqlitefile
    
    
    
    
    

class mzmlsql_reader(object):
    """Random-access reader for sql DBs written by mzmlToSqlite."""
    def __init__(self, sqlitefile):
        self.connection = sqlite.connect(sqlitefile)
        self.cursor = self.connection.cursor()
    
    def scan(self, index):
        self.cursor.execute('SELECT data FROM spectra WHERE ind=%s' % str(index))
        return demarshal(self.cursor.fetchone()[0])
    
    def xic(self, index):
        self.cursor.execute('SELECT data FROM chromato WHERE ind=%s' % str(index))
        return demarshal(self.cursor.fetchone()[0])
    
    def scans_by_mz(self, mzstart, mzstop):
        self.cursor.execute('SELECT ind, mz, data FROM spectra WHERE mz >= %s AND mz <= %s' % (mzstart, mzstop))
        return [(x[0], x[1], demarshal(x[2])) for x in self.cursor.fetchall()]
    
    def scans_by_rt(self, rtstart, rtstop):
        self.cursor.execute('SELECT ind, rt, data FROM spectra WHERE rt >= %s AND rt <= %s' % (rtstart, rtstop))
        return [(x[0], x[1], demarshal(x[2])) for x in self.cursor.fetchall()]
    
    def scan_manifest(self):
        self.cursor.execute('SELECT ind, mz, rt FROM spectra')
        return list(self.cursor.fetchall())
    
    def xic_manifest(self):
        self.cursor.execute('SELECT ind, mzstart, mzstop, rtstart, rtstop FROM chromato')
        return list(self.cursor.fetchall())
    
    def total_xic(self):
        self.cursor.execute('SELECT data FROM chromato WHERE ind=0 AND startmz=0 AND stopmz=0 AND startrt=0 AND stopmz=0')
        return demarshal(self.cursor.fetchone()[0])
    
    def close():
        self.connection.close()

    

def iterate_spectra_simple(xmlfile):
    parser = xml.iterparse(gzOptOpen(xmlfile))
    for evt, obj in parser:
        if obj.tag == ns('spectrum'):
            yield readSpectrumXML(obj)
            obj.clear()

def iterate_chromatograms(xmlfile):
    parser = xml.iterparse(xmlfile)
    for evt, obj in parser:
        if obj.tag == ns('spectrum'):
            yield readChromatoXML(obj)
            obj.clear()


#def reduceObject(descent, obj):
    #lineage = [obj.tag]
        

#class iterate_spectra(object):
    #def __init__(self, xmlfile, get_attrbiutes = True):
        #self.parser = xml.iterparse(xmlfile)
        #self.get_attributes = get_attributes
        #self.attributes = {}
        #self.descent = {}
        
    #def __iter__(self, xmlfile):
        #for evt, obj in self.parser:
            #self.descent.update({obj:x for x in obj})
            #if obj.tag == ns('spectrum'):
                #yield readSpectrumXML(obj)
            #elif self.get_attributes:
                #objatts = reduceObject(self.descent, obj)
                #self.attributes.update(objatts)
    
        
