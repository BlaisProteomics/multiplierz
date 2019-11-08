from multiplierz.spectral_process import deisotope_reduce_scan, peak_pick
from multiplierz.spectral_process import centroid as centroid_func # Distinct from 'centroid' argument.
from multiplierz.mass_biochem import remove_protons
from multiplierz.internalAlgorithms import floatrange, aggregate_points
from multiplierz.mgf import standard_title_write
from collections import defaultdict
from numpy.linalg import solve     
import os
import bisect
from itertools import chain
from numpy import average
import re

RAW_CAL_MASS = 445.120025
calibrant_tolerance = 0.003
MS1_LENGTH = 0.5

def compile_correction_matrix(channel_corrections, labels):
    assert labels, "No isobaric tag specified, but correction matrix given."
    from numpy import zeros
    cormat = zeros(shape = (len(labels), len(labels)))
    for i, froml in enumerate(labels):
        for j, tol in enumerate(labels):
            cormat[j, i] += channel_corrections.get(froml, {}).get(tol, 0)
        if froml not in channel_corrections[froml]:
            cormat[i, i] += 100.0 - sum(channel_corrections.get(froml, {}).values())
    return cormat.transpose() / 100



def parse_prec_info(precfile):
    precs = {}
    for line in open(precfile, 'r'):
        accessid, mz, chg = line.strip().strip('__').split()[:3]
        precs[float(accessid)] = float(mz), int(chg)
    return precs

class _extractor_(object):
    """
    Internal class for extracting MGF files; users should only be concerned with
    the higher-level extract() function.
    """
    def __init__(self, data, filename, default_charge, centroid,
                 scan_type, deisotope_and_reduce_charge, derive_precursor_via,
                 maximum_precursor_mass, long_ms1,
                 deisotope_and_reduce_MS1_args, deisotope_and_reduce_MS2_args,
                 min_mz, precursor_tolerance, isobaric_labels, label_tolerance,
                 channel_corrections,
                 prec_info_file = None,
                 region_based_labels = False):
        self.data = data
        self.filename = filename
        self.default_charge = default_charge
        self.centroid = centroid
        self.scan_type = scan_type
        self.deisoreduce = deisotope_and_reduce_charge
        self.derive_precursor_via = derive_precursor_via
        self.deisoreduce_MS1_args = deisotope_and_reduce_MS1_args
        self.deisoreduce_MS2_args = deisotope_and_reduce_MS2_args
        self.min_mz = min_mz
        self.precursor_tolerance = precursor_tolerance
        #self.isobaric_labels = isobaric_labels
        self.label_tolerance = label_tolerance
        self.channel_corrections = channel_corrections
        self.long_ms1 = long_ms1
        self.maximum_mass = maximum_precursor_mass
        self.region_based_labels = region_based_labels
        
        self.prec_info = parse_prec_info(prec_info_file) if prec_info_file else None
        
        self.set_isobaric_labels(isobaric_labels)
        self.initialize_scan_info()
        if channel_corrections is not None:
            if isinstance(channel_corrections, dict):
                self.correction_matrix = compile_correction_matrix(channel_corrections,
                                                                   list(zip(*self.labels))[0])
            else:
                self.correction_matrix = channel_corrections
        else:
            self.correction_matrix = None
        self.calibrant = RAW_CAL_MASS
        self.MS1_scan_batch = {}
        
        if ('tolerance' not in self.deisoreduce_MS2_args
            or not self.deisoreduce_MS2_args['tolerance']):
            self.deisoreduce_MS2_args['tolerance'] = precursor_tolerance        
        
   
    def set_isobaric_labels(self, isobaric_labels):
        if not isobaric_labels:
            labels = []
        elif isobaric_labels == 4 or isobaric_labels == '4plex':
            labels = list(zip(['114', '115', '116', '117'], [114.11123, 115.10826, 116.11162, 117.11497]))
        elif isobaric_labels == 6 or isobaric_labels == '6plex':
            labels = list(zip(['126', '127', '128', '129', '130', '131'],
                         [126.127726, 127.131081, 128.134436,
                          129.137790, 130.141145, 131.138180]))
        elif isobaric_labels == '6plex2':
            labels = list(zip(['126', '127', '128', '129', '130', '131'],
                         [126.127726, 127.124761, 128.128116,
                          129.131471, 130.134825, 131.138180]))            
        elif isobaric_labels == 8 or isobaric_labels == '8plex':
            labels = list(zip(['113', '114', '115', '116', '117', '118', '119', '121'],
                         [113.11, 114.11, 115.11, 116.11, 117.12, 118.12, 119.12, 121.12]))
        elif isobaric_labels == 10 or isobaric_labels == '10plex':
            labels = list(zip(['126', '127N', '127C', '128N', '128C', 
                          '129N', '129C', '130N', '130C', '131'],
                         [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                          129.131471, 129.137790, 130.134825, 130.141145, 131.138180]))
        elif isobaric_labels == 11:
            labels = list(zip(['126', '127N', '127C', '128N', '128C',
                          '129N', '129C', '130N', '130C', '131N', '131C'],
                         [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                          129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144499]))
        else:
            raise NotImplementedError("Labels of type %s not recognized.\n"
                                        "Should be one of [4,6,8,10] or None." % isobaric_labels)
        
        if isobaric_labels in [10, 11, '10plex', '11plex'] and not self.region_based_labels:
            assert self.label_tolerance < 0.005, ("label_tolerance must be lower "
                                                  "than 0.005 for 10-plex experiments! "
                                                  "(Currently %s)"
                                                  % self.label_tolerance)        
        self.labels = labels
    
    def initialize_scan_info(self):
        self.scanInfo = self.data.scan_info()
        self.ms1_list = sorted([x[2] for x in self.scanInfo if x[3] == 'MS1'])
        
        if self.filename.lower().endswith('.raw'): # May also exist for WIFF?
            self.filters = dict(self.data.filters())
            
            # For RAW files only, there's the option to filter by a given
            # scan type.  (It would be more efficient in many cases to
            # actually split files in a single run, though.)
            scan_type = self.scan_type
            if scan_type and isinstance(scan_type, str):
                typestr = "@%s" % scan_type.lower()
                self.scanInfo = [x for x in self.scanInfo if x[3] == 'MS1' or
                                 typestr in self.filters.get(x[0], '')]
        else:
            self.filters = None
            assert not scan_type, "Scan type filtering only enabled with .RAW format files."             
    
    def get_precursor(self, mz, charge):
        if self.possible_precursors == None:
            self.calculate_precursors()        
        try:
            return min([x for x in self.possible_precursors
                        if (charge == None or x[1] == charge)],
                       key = lambda x: abs(x[0] - mz))
        except ValueError:
            return None, None    

    def read_labels(self, scan):
        partscan = [x for x in scan if x[0] < self.labels[-1][1] + 3]
        if not partscan:
            return dict([(str(l), '0') for l in zip(*self.labels)[0]])

        ## This should probably actually sum all points within
        ## the tolerance range.
        #scan_values = {}
        #for label, mz in self.labels:
            #nearpt = min(partscan, key = lambda x: abs(x[0] - mz))
            #if abs(nearpt[0] - mz) < self.label_tolerance:
                #scan_values[str(label)] = nearpt[1]
            #else:
                ## Report noise value?  (Easier to determine bad reads this way though.)
                #scan_values[str(label)] = 0
                
        # This sums all points within the tolerance range; resistant to things
        # like "6plex" mixtures made with mixed isotopalogues, less good for 
        # distinguishing between 10- and 11-plex ions (use region_based_labels = True).
        scan_values = {}
        for label, mz in self.labels:
            zone = [x for x in partscan if abs(x[0] - mz) < self.label_tolerance]
            scan_values[str(label)] = sum([x[1] for x in zone])
        
        if self.channel_corrections is not None:
            reporter_vector = [scan_values[x[0]] for x in self.labels]
            corrected_vector = solve(self.correction_matrix, reporter_vector)
            return dict(list(zip(list(zip(*self.labels))[0],
                            [max(x, 0) for x in corrected_vector.flatten()])))
        else:
            return scan_values    
    
    def read_labels_by_region(self, scan):
        # For wonky 10- and 11-plex scans, a more robust way of getting the desired
        # points.
        
        # Label strings must be in order of ascending MZ.
        if len(self.labels) == 10:
            zones = [(126.127726, ('126', ), (126.127726,)),
                     (127.127921, ('127N', '127C'), (127.124761, 127.131081)),
                     (128.131276, ('128N', '128C'), (128.128116, 128.134436)),
                     (129.1346305, ('129N', '129C'), (129.131471, 129.137790)),
                     (130.137985, ('130N', '130C'), (130.134825, 130.141145))]
        elif len(self.labels) == 11:
            zones = [(126.127726, ('126', ), (126.127726,)),
                     (127.127921, ('127N', '127C'), (127.124761, 127.131081)),
                     (128.131276, ('128N', '128C'), (128.128116, 128.134436)),
                     (129.1346305, ('129N', '129C'), (129.131471, 129.137790)),
                     (130.137985, ('130N', '130C'), (130.134825, 130.141145)),
                     (131.1413395, ('131N', '131C'), (131.138180, 131.144499))]
        else:
            raise Exception
        
        values = {}
        for center, labels, labelpts in zones:
            region = center - self.label_tolerance, center + self.label_tolerance
            pts = [x for x in scan if region[0] <= x[0] <= region[1]]
            tops = sorted(pts, key = lambda x: x[1], reverse = True)[:len(labels)]
            tops.sort(key = lambda x: x[0]) # Sort by MZ to match order of labels.
            tops += [(0, 0)] * len(labels) # In case there's missing values.
            #(No attempt to figure out *which* are missing.)
            values.update(list(zip(labels, [x[1] for x in tops])))
        
        if self.channel_corrections is not None:
            reporter_vector = [values[x[0]] for x in self.labels]
            corrected_vector = solve(self.correction_matrix, reporter_vector)
            return dict(list(zip(list(zip(*self.labels))[0],
                            [max(x, 0) for x in corrected_vector.flatten()])))
        else:
            return values         
        
        
    
    def raw_scan_recalibration(self, scan):
        cal_region = [x for x in scan if abs(self.calibrant - x[0]) < calibrant_tolerance]
        if cal_region:
            self.calibrant = max(cal_region, key = lambda x: x[1])[0]
            cal_factor = RAW_CAL_MASS / self.calibrant
            return [(x[0] * cal_factor, x[1]) for x in scan]
        else:
            return scan
    
    def get_long_MS1_byProfile(self, scanNum):
        # This seems to take a really long time.  Profile-mode scans are heavy!
        ms1_index = bisect.bisect_left(self.ms1_list, scanNum)
        scannumbers = [self.ms1_list[ms1_index + i] for i
                       in [-1, 0, 1]
                       if ms1_index + i >= 0 and ms1_index + i < len(self.ms1_list)]

        ms1s = []
        for scannum in scannumbers:
            if scannum in self.MS1_scan_batch:
                ms1s.append(self.MS1_scan_batch[scannum])
            else:
                scan = self.data.scan(scannum, centroid = False)
                self.MS1_scan_batch[scannum] = scan
                ms1s.append(scan)
        to_del = [k for k in self.MS1_scan_batch if k <= scannumbers[0]]
        for scannum in to_del:
            del self.MS1_scan_batch[scannum]
        
        # MS1 profile scans aren't consistent across the entire MZ range,
        # apparently dependent upon whether there is signal in a particular
        # location; we don't want to cover the entire range evenly, since
        # this would expand the data and slow things down, but we also don't
        # want to miss/mis-assign signal that occurs in only one scan of the
        # batch.
        long_mzs = set()
        for b in ms1s:
            long_mzs.update(round(x[0], 4) for x in b)
            
        long_ms1 = []
        inds = [0]*len(ms1s)
        for lmz in sorted(long_mzs):
            sumint = 0
            for j in range(len(ms1s)):
                bat, ind = ms1s[j], inds[j]
                while ind < len(bat) and bat[ind][0] < lmz:
                    sumint += bat[ind][1]
                    inds[j] += 1
                    ind = inds[j]
            long_ms1.append((lmz, sumint))
            
        return centroid_func(long_ms1)
    
    def get_long_MS1_old(self, scanNum): # By centroid mode subscans.
        if scanNum == 0: # List slicing doesn't work right otherwise.
            scanNum = 1
        ms1_index = bisect.bisect_left(self.ms1_list, scanNum)
        scannumbers = self.ms1_list[ms1_index-1:ms1_index+2]

        ms1s = []
        for scannum in scannumbers:
            if scannum in self.MS1_scan_batch:
                ms1s.append(self.MS1_scan_batch[scannum])
            else:
                scan = self.data.scan(scannum, centroid = True)
                self.MS1_scan_batch[scannum] = scan
                ms1s.append(scan)
                
        if len(self.MS1_scan_batch) > 100:
            to_del = [k for k in self.MS1_scan_batch if k <= scannumbers[0]]
            for scannum in to_del:
                del self.MS1_scan_batch[scannum]
            
        long_ms1 = []
        inds = [0]*len(ms1s)
        agg = []
        while all(inds[j] < len(ms1s[j]) for j in range(len(ms1s))):
            next_i = min(list(range(len(ms1s))),
                         key = lambda j: ms1s[j][inds[j]][0])
            next_pt = ms1s[next_i][inds[next_i]]
            inds[next_i] += 1
            if (not agg) or abs(next_pt[0] - agg[0][0]) < 0.005:
                agg.append(next_pt)
            else:
                # Believe it or not, breaking out 1- or 2-length cases makes the
                # whole extraction process substantially faster.
                if len(agg) == 1:
                    long_ms1.append(agg[0])
                elif len(agg) == 2:
                    long_ms1.append((
                        ((agg[0][0] * agg[0][1] + agg[1][0] * agg[1][1]) / (agg[0][1] + agg[1][1])),
                        agg[0][1] + agg[1][1]))
                else:
                    mzs, ints = list(zip(*agg))
                    long_ms1.append((average(mzs, weights = ints), sum(ints)))
                agg = [next_pt]
                
        return long_ms1
                
    def get_long_MS1(self, scanNum): # By centroid mode subscans.
        if scanNum == 0: # List slicing doesn't work right otherwise.
            scanNum = 1
        ms1_index = bisect.bisect_left(self.ms1_list, scanNum)
        scannumbers = self.ms1_list[ms1_index-1:ms1_index+2]

        ms1s = []
        for scannum in scannumbers:
            if scannum in self.MS1_scan_batch:
                ms1s.append(self.MS1_scan_batch[scannum])
            else:
                scan = self.data.scan(scannum, centroid = True)
                self.MS1_scan_batch[scannum] = scan
                ms1s.append(scan)
                
        if len(self.MS1_scan_batch) > 100:
            to_del = [k for k in self.MS1_scan_batch if k <= scannumbers[0]]
            for scannum in to_del:
                del self.MS1_scan_batch[scannum]                
        
        agg_points = aggregate_points(list(chain(*ms1s)), MAX_WIDTH = 0.005)
        long_ms1 = []
        for agg in agg_points:
            if len(agg) == 1:
                long_ms1.append(agg[0])
            elif len(agg) == 2:
                long_ms1.append((
                    ((agg[0][0] * agg[0][1] + agg[1][0] * agg[1][1]) / (agg[0][1] + agg[1][1])),
                    agg[0][1] + agg[1][1]))
            else:
                mzs, ints = list(zip(*agg))
                long_ms1.append((average(mzs, weights = ints), sum(ints)))
        return long_ms1

   
    def calculate_precursors(self):
        if not self.lastMS1ScanName:
            lastMS1 = []        
        else:
            if self.data.format == 'raw':
                if self.long_ms1:
                    lastMS1 = self.get_long_MS1(self.lastMS1ScanName)
                    # Is calibration valid in this case?
                else:
                    try:
                        lastMS1 = self.data.lscan(self.lastMS1ScanName)
                    except IOError: # No lscan data in file.
                        lastMS1 = self.data.scan(self.lastMS1ScanName, 
                                                 centroid = True)
                        lastMS1 = [(mz, i, 0, 0) for mz, i in lastMS1]
                lastMS1 = self.raw_scan_recalibration(lastMS1)
            else:
                try:
                    lastMS1 = self.data.scan(self.lastMS1ScanName,
                                             centroid = True)
                except NotImplementedError:
                    lastMS1 = centroid_func(self.data.scan(self.lastMS1ScanName))      
                    
            envelopes = peak_pick(lastMS1, **self.deisoreduce_MS1_args)[0]
            self.possible_precursors = sum([[(x[0][0], c) for x in xs]
                                            for c, xs in list(envelopes.items())], [])    
        
    
    def run(self):
        self.inconsistent_precursors = 0
        self.scans_written = 0
        
        self.lastMS1ScanName = None
        self.possible_precursors = None
        for time, mz, scanNum, scanLevel, scanMode in self.scanInfo:
            scanName = scanNum if isinstance(scanNum, int) else time
            
            if scanLevel == 'MS1':
                self.lastMS1ScanName = scanName
                self.possible_precursors = None
                continue
            elif scanLevel == 'MS3':
                continue
            elif self.lastMS1ScanName == None:
                continue                
            
            # Each file type handles centroiding differently (or not at all.)
            if self.data.format == 'raw':
                scan = self.data.scan(scanName, centroid = self.centroid)
                scan = self.raw_scan_recalibration(scan)
            elif self.data.format == 'wiff':
                # explicit_numbering, of course, can't be active here.
                scan = self.data.scan(scanName)
                if self.centroid:
                    scan = centroid_func(scan)
            elif self.data.format == 'd':
                scan = self.data.scan(scanName, centroid = self.centroid)
                if self.centroid and not scan:
                    # mzAPI.D returns empty if centroid data is not present in
                    # the file, but that can be corrected by external centroiding.
                    scan = centroid_func(self.data.scan(scanName, centroid = False))
            else:
                raise NotImplementedError("Extractor does not handle type %s"
                                            % self.data.format)
            
            if self.filters and not mz:
                mz = float(self.filters[time].split('@')[0].split(' ')[-1])            
            
            
            if self.prec_info:
                accessid = self.data.extra_info(scanName)['Access Id']
                assert accessid, scanName
                mzP, chargeP = self.prec_info[accessid]
                assert mzP and chargeP, (scanName, accessid, mzP, chargeP)
                mz = mzP
                charge = chargeP
            else:
                mzP = None
                chargeP = None
                if ("scanPrecursor" in dir(self.data) and 
                    self.derive_precursor_via in ['All', 'Thermo']):
                    assert isinstance(scanName, int)
                    mzP, chargeP = self.data.scanPrecursor(scanName)
                    
                if (not mzP) or self.derive_precursor_via in ['Direct']: # 'and derive_precursor_via not in ['Thermo']', except why would you?
                    # .scanPrecursor sometimes returns charge and not mzP.
                    mzP, chargeP = self.get_precursor(mz, chargeP)
                    if not mzP:
                        # Release presumed charge possibly obtained from scanPrecursor.
                        mzP, chargeP = self.get_precursor(mz, None)
                        if mz and chargeP:
                            self.inconsistent_precursors += 1     
                
                if mzP and (abs(mz - mzP) < 2 or not mz): 
                    mz = mzP
                    charge = chargeP
                else:
                    charge = self.default_charge
                
                if not charge:
                    charge = self.default_charge                
            
            if not mz:
                import warnings
                errmgf = os.path.abspath(self.filename)
                warnings.warn('Unable to recover all precursor masses from %s' % errmgf)
            else:    
                if (self.maximum_mass and
                    remove_protons(mz, charge) > self.maximum_mass):
                    continue
                    
                if self.labels:
                    if self.region_based_labels:
                        scan_labels = self.read_labels_by_region(scan)
                    else:
                        scan_labels = self.read_labels(scan)
                else:
                    scan_labels = {}
                
                if self.filters:
                    try:
                        filt = self.filters[time]
                    except KeyError:
                        # Very rare occurrence; mismatch between list of filters
                        # and list of scans.  Presumed to be due to data corruption.
                        print(("Corrupt scan at %s; ignoring." % time))
                        continue
                    
                    detector = re.search('[FI]TMS', filt)
                    frag_e = re.search(r'@(hcd|etd|cid)([0-9]+.[0-9]+)', filt)
                    if detector:
                        scan_labels['type'] = detector.group(0)
                    if frag_e:
                        try:
                            scan_labels['ce'] = frag_e.group(2)
                        except IndexError:
                            pass
                        try:
                            scan_labels['frgmntr'] = frag_e.group(1) # Better key for that?
                        except IndexError:
                            pass
                    
                
                title = standard_title_write(self.filename, rt = time, mz = mz,
                                             mode = scanMode, scan = scanNum,
                                             **scan_labels)
            
                if self.deisoreduce and self.centroid:
                    scan = deisotope_reduce_scan(scan, **self.deisoreduce_MS2_args)  
                scan = [x for x in scan if x[0] > self.min_mz]
                assert charge, title         
                
                yield scan, title, mz, charge
                self.scans_written += 1





