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

import math
import os
import re

from collections import defaultdict

import matplotlib

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.backends.backend_svg import FigureCanvasSVG
from matplotlib.ticker import ScalarFormatter

from multiplierz.mass_biochem import generate_labels, mz as calc_mz


# Dictionary of color codes for the different ions. Keeping it in
# one place makes it simpler to change. Doubly-charged ions are a
# lighter version of the singly-charged color. 'MH' is precursor,
# 'mult' is multiple different ions, '?' is unknown ion type
ion_color_dict = {  'b': '#0000FF',   'b++': '#9696FF',
                    'y': '#FF0000',   'y++': '#FF9696',
                    'c': '#0000FF',   'c++': '#9696FF',
                  'z+1': '#FF0000', 'z+1++': '#FF9696',
                  'z+2': '#FF0000', 'z+2++': '#FF9696',
                   'MH': '#CDCD00',  'mult': '#EE30A7',
                    '?': '#FFA500'}

# This dictionary takes the place of a bunch of if-else stuff which
# was used to determine coordinates
ion_dict = {'b': (0.95, '_'), 'b++': (0.97, '_'),
            'y': (0.85, '_'), 'y++': (0.83, '_')}

# This one is the special case for the ETD-TRAP
# Do we need both? If we put the b-ion entries into this one, it does double duty.

# the question is--does an ETD-TRAP have b ions that should be ignored,
# and does a non-ETD have c/z ions that should be ignored?

# Answer: Mascot will never return such ions to us, and when we're
# doing this ourselves we can specify which ions to pay attention to.
# So, we can merge the dictionaries. But that will require some decisions
# on spacing--need all the locations visible.

ETDTRAP_ion_dict = {  'c': (0.95, '_'),   'c++': (0.97, '_'),
                      #'y': (0.87, '_'),   'y++': (0.85, '_'),
                    'z+1': (0.85, '_'), 'z+1++': (0.83, '_')}

def add_pep_text(fig, sequence, ion_matches=None,
                 phospho=None, oxid=None, silac=None, instrument='ESI-TRAP'):
    '''Draws the sequence on top of an ms2 plot. This function
    takes an existing figure, which should be sized so that
    there is room at the top for these graphics.'''

    phospho = set(phospho or [])
    oxid = set(oxid or [])
    silac = set(silac or [])

    w,h = fig.get_size_inches()

    font_size = min(float(w) * 72.0 / (len(sequence)+2), # 72 points to an inch, w is in inches
                    17.28) # 17.28 points is the 'x-large' font size

    font_step = 1.0 / (len(sequence) + 2)
    font_offset = 1.5 * font_step

    xcoords = [(font_offset + i * font_step) for i in range(len(sequence))]

    draw_ticks = len(sequence) >= 10
    tick_font_size = font_size / 2.5 if font_size > 10.0 else font_size / 2.0

    for i,AA in enumerate(sequence):
        if i in phospho:
            color = 'green'
        elif i in oxid:
            color = 'blue'
        elif i in silac:
            color = 'grey'
        else:
            color = 'black'

        # draw the amino acid in the right color
        t = fig.text(xcoords[i],
                     0.9,
                     AA.upper(),
                     color=color,
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=font_size)

        # if we have ion matches, draw this AA's ions. we use _ and ~ characters
        if ion_matches:
            if draw_ticks and i > 0:
                # draw b-ion tick
                if i % 5 == 0:
                    fig.text(xcoords[i] - font_step / 2.0,
                             0.925,
                             '%d' % i,
                             color='black',
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=tick_font_size)
                # draw y-ion tick
                if (len(sequence) - i) % 5 == 0:
                    fig.text(xcoords[i] - font_step / 2.0,
                             0.875,
                             '%d' % (len(sequence) - i),
                             color='black',
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=tick_font_size)

            for j in ion_matches[i]:
                # get the y-coordinate, the color, and the character to print
                if instrument == 'ETD-TRAP' and j in ETDTRAP_ion_dict:
                    (ycoord, c) = ETDTRAP_ion_dict[j]
                    colr = ion_color_dict[j]
                elif instrument != 'ETD-TRAP' and j in ion_dict:
                    (ycoord, c) = ion_dict[j]
                    colr = ion_color_dict[j]
                else:
                    continue
                fig.text(xcoords[i],
                         ycoord,
                         c,
                         color=colr,
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontsize=font_size)

    return xcoords


def translate_labels(labels, seq_len):
    '''Tranlates the labels from mascot into the label list that
    our image-maker wants'''

    included = set(('b', 'y', 'c', 'z+1'))

    # complicated regex matches ion, location, and (optional) charge and theoretical mass/mass error
    ion_re = re.compile(r'(.+)\((\d+)\).*?([\+\-]{2})?(?: \[.+\])?$')

    ion_info = [set() for i in range(seq_len)]
    ion_coords = [defaultdict(list) for i in range(seq_len)]

    for mass,lbl in labels:
        m = ion_re.match(lbl)
        if m:
            # soon: going to include every ion we saw, even though
            # multiplierz doesn't currently deal with most of them
            for i in included:
                if m.group(1).startswith(i):
                    #i = m.group(1)
                    j = int(m.group(2))
                    if i[0] in ('a','b','c'):
                        ion_info[j-1].add(i + (m.group(3) or ''))
                        ion_coords[j-1][mass].append(lbl)
                    elif i[0] in ('x','y','z'):
                        ion_info[-j].add(i + (m.group(3) or ''))
                        ion_coords[-j][mass].append(lbl)
                    break

    return ion_info, ion_coords


def get_color(labels):
    '''Takes a list of labels and returns a color for each point,
    depending on whether the labels are all of one set, or mixed'''

    # regex matches ion (but not '0' or '*') and optional charge
    lbl_re = re.compile(r'(?:(MH)|(.+?))(?(2)[\*0]?\(\d+\).*?)([\+\-]{2,})?')
    lbl_set = list(set(lbl_re.match(lbl).groups() for lbl in labels))
    if len(lbl_set) == 1:
        if lbl_set[0][0] == 'MH':
            return ion_color_dict['MH'] # yellow for precursor

        i = lbl_set[0][1] + (lbl_set[0][2] or '')

        if i in ion_color_dict:
            return ion_color_dict[i]
        else:
            return ion_color_dict['?'] # orange for unknown ion label
    else:
        return ion_color_dict['mult'] # purple for multiple ion types--fairly uncommon


def make_ms1_im(save_file, mz, xy, scan_mode, pm_scanDot, title=None, fmt='PNG', im_size=(8.0,6.0)):
    '''Plots an MS scan, centered on the target mz'''

    fig = Figure()

    _make_ms1(fig, mz, xy, scan_mode, pm_scanDot, title=title)

    dpi = int(sum(im_size)/14.0 * 100.0)
    if fmt == 'PDF':
        FigureCanvasPdf(fig).print_pdf('.'.join([save_file, fmt.lower()]))
    else:
        FigureCanvasAgg(fig).print_figure('.'.join([save_file, fmt.lower()]), dpi=dpi, format=fmt)


def make_xic_im(save_file, mz, xic_time, xic_int, scan_dot, bin_times, bin_ints, title=None, fmt='PNG', im_size=(8.0,6.0)):
    '''Plots a XIC, with labels for the primary MSMS scan and any neighboring scans
    in the same time and mz area'''

    fig = Figure()

    _make_xic(fig, mz, xic_time, xic_int, scan_dot, bin_times, bin_ints, title=title)

    dpi = int(sum(im_size)/14.0 * 100.0)
    if fmt == 'PDF':
        FigureCanvasPdf(fig).print_pdf('.'.join([save_file, fmt.lower()]))
    else:
        FigureCanvasAgg(fig).print_figure('.'.join([save_file, fmt.lower()]), dpi=dpi, format=fmt)


def make_ms2_im(save_file, scan, scan_mode, peptide, labels=None, ion_list=('b', 'y'),
                charge='', score='', title=None, fmt='PNG', im_size=(8.0,6.0),
                label_top=50, filter_labels=True, y_range = None, x_range = None,
                figure_labels = {}, tolerance = None, **settings):
    '''Creates an image for the MS/MS scan with theoretical fragments overlayed and sequence on top
    scan = [(mz, int)...]
    scan_mode = 'p' or 'c'
    peptide = PEPpTIDE-[100.0] (note that mods should be included)
    labels = [(mz, label)...] or None (labels will be generated)
    ion_list = ('b', 'y') (this list can be much longer if desired)
    charge = 2
    score = 50.0
    '''
    # initialize defaults and update with optional arguments
    _settings = dict(show_theor_mz=True, ms2_mz_figs=2, show_mass_error=False,
                     mass_error_figs=2, mass_error_units='ppm')
    _settings.update(settings)

    fig = Figure()

    #_make_ms2(fig, scan, scan_mode, peptide, labels,
              #ion_list, charge, score, label_top,
              #title=title, filter_labels=filter_labels) #, pretty_print)
    
    a, peptide_annotations = _make_ms2(fig, scan, scan_mode, peptide, labels,
              ion_list, charge, score, label_top,
              title=title, filter_labels=filter_labels, y_range = y_range, x_range = x_range, tolerance= tolerance, **_settings)
    
    # print a
    out = open(save_file + '.txt', mode = 'w')
    for m in a:
        for ion in m[2]:
            print >> out,  ion, ':\t', m[0:2]
    out.close()

    
    if figure_labels:
        ax = fig.get_axes()[0]
        bot, top = ax.get_ylim()
        left, right = ax.get_xlim()
        hm = (right - left)*0.05
        vm = (top - bot)*0.05 
        text = 'Precursor %s m/z\nScan %s' % (figure_labels['Precursor MZ'], figure_labels['Scan Number'])
        ax.text(right-hm, top-vm, text,
                verticalalignment = 'top', horizontalalignment = 'right')
    
    #for label, value in figure_labels.items():
        #if label == 'Precursor MZ':
            #text = 'Precursor %s m/z' % value
            #ax.text(left+hm, top+vm, text, verticalalignment = 'top', horizontalalignment = 'left')
        #elif label == 'Scan Number':
            #text = 'Scan %s' % value
            #ax.text(right+hm, top+vm, text, verticalalignment = 'top', horizontalalignment = 'right')
        #else:
            #raise Exception, "Tyopes?"
    
         

    dpi = int(sum(im_size)/14.0 * 100.0)
    if fmt == 'PDF':
        FigureCanvasPdf(fig).print_pdf('.'.join([save_file, fmt.lower()]))
    else:
        FigureCanvasAgg(fig).print_figure('.'.join([save_file, fmt.lower()]), dpi=dpi, format=fmt)


def make_mascot_ms2_im(save_file, web_img, peptide, charge, score, ion_info=None,
                       phospho=None, oxid=None, instrument='ESI-TRAP', fmt='PNG', im_size=(8.0,6.0)):
    '''Displays the Mascot MS2 image with labeled peaks, along with our own sequence
    image of color-coded mods and ions'''

    fig = Figure()

    axes = fig.add_axes([0.125,  0.1,  0.775,  0.8])

    axes.imshow(matplotlib.image.imread(web_img), origin='lower', interpolation='nearest')
    axes.set_title("Charge: %s    Score: %s" % (charge, score))

    axes.get_xaxis().set_visible(False)
    axes.get_yaxis().set_visible(False)

    # the negative position is to crop off a blue line at the bottom
    new_pos = [0.0, -0.02, 1.0, 0.77]
    axes.set_position(new_pos)

    add_pep_text(fig, peptide, ion_info, phospho, oxid, instrument)

    dpi = int(sum(im_size)/14.0 * 100.0)
    if fmt == 'PDF':
        FigureCanvasPdf(fig).print_pdf('.'.join([save_file, fmt.lower()]))
    else:
        FigureCanvasAgg(fig).print_figure('.'.join([save_file, fmt.lower()]), dpi=dpi, format=fmt)


def make_mirror_im(save_file, scan1, scan2, title="Mirror Image", fmt='PNG', im_size=(8.0,6.0)):
    '''Creates a mirror plot of two scans--i.e. plots scan2 upside-down
    on the same x-axis as scan1, allowing them to be compared easily'''
    # Flip scan 2
    scan2 = [(mz, -1*intensity) for mz,intensity in scan2]

    #Re-format save_file string so that matplot likes the syntax
    save_file = re.sub('\\\\','\\\\\\\\',save_file)

    x_scan1 = [x1 for (x1,y1) in scan1]
    y_scan1 = [y1 for (x1,y1) in scan1]

    x_scan2 = [x1 for (x1,y1) in scan2]
    y_scan2 = [y1 for (x1,y1) in scan2]

    fig = Figure()

    axes = fig.add_axes([0.125,  0.1,  0.775,  0.8])

    axes.plot(x_scan1, y_scan1, c='r')
    axes.plot(x_scan2, y_scan2, c='k')

    axes.set_title(title)

    dpi = int(sum(im_size)/14.0 * 100.0)
    if fmt == 'PDF':
        FigureCanvasPdf(fig).print_pdf('.'.join([save_file, fmt.lower()]))
    else:
        FigureCanvasAgg(fig).print_figure('.'.join([save_file, fmt.lower()]), dpi=dpi, format=fmt)


def make_venn(fig, A, B, AB, A_label=None, B_label=None, title=None, eps=0.001):
    '''Plot a proportional 2-set Venn diagram. A and B are the sizes of the two sets,
    AB is the size of the intersection, and eps is an error margin for the proportional
    placement. E.g. if eps is 0.01 then the areas of the plot will be accurate to ~1%.

    A lower eps will give a more accurate plot at the expense of longer running time.
    The method uses a bisecting search algorithm to find the right proportions.'''

    def get_wAB(rA, rB, d):
        '''Given the radii of the two circles and the offset, calculate
        the area of intersection.'''

        alpha = 2 * math.acos((d**2 + rA**2 - rB**2) / (2*rA*d))
        beta = 2 * math.acos((d**2 + rB**2 - rA**2) / (2*rB*d))

        wAB = 0.5 * rA**2 * (alpha - math.sin(alpha)) + 0.5 * rB**2 * (beta - math.sin(beta))

        return wAB

    def calc_venn(A, B, AB, eps=0.0001):
        '''Given the sizes of the two sets and their intersection, figure
        out the spacing within ~eps error.'''

        def bin_search(d0, d1):
            '''binsecting search algorithm to find the right distance.'''
            dh = 0.5 * (d0 + d1)

            w = get_wAB(rA, rB, dh)

            if abs((w - wAB) / wAB) < eps:
                return dh
            elif w < wAB:
                return bin_search(d0, dh)
            else:
                return bin_search(dh, d1)

        max_w = max((A,B))

        # normalize to 1.0
        wA = 100.0 * A / max_w
        wB = 100.0 * B / max_w
        wAB = 100.0 * AB / max_w

        rA = math.sqrt(wA / math.pi)
        rB = math.sqrt(wB / math.pi)

        # minimum and maximum possible answers
        dmin,dmax = sorted((abs(rA - rB), rA + rB))

        if wAB == 0.0:
            d = dmax
        elif wAB == 100.0:
            d = dmin
        else:
            d = bin_search(dmin, dmax)

        return rA, rB, d

    self.figure.clear()
    # big axes to used up as much space as we can
    axes = self.figure.add_axes([0.05, 0.05, 0.9, 0.85])

    rA,rB,d = calc_venn(A, B, AB, eps)

    max_r = max((rA, rB))

    # half the distance, for placing the circles
    d2 = 0.5 * d

    # offset for a label
    a_off = 0.86 * rA
    # offset for b label
    b_off = 0.86 * rB

    # midpoint of a
    a_mid = (-rA - rB) / 2
    # midpoint of b
    b_mid = (rB + rA) / 2
    # midpoint of intersection
    ab_mid = (rA - rB) / 2

    # x lim so that there's just enough space
    x_lim = 1.05 * (d2 + max_r)
    # y lim with padding for label
    y_lim = 1.1 * max_r

    c1 = Circle((-d2,0), rA, alpha=0.2, fc='red')
    c2 = Circle((d2,0), rB, alpha=0.2, fc='blue')

    axes.add_patch(c1)
    axes.add_patch(c2)

    axes.set_xlim(-x_lim, x_lim)
    axes.set_ylim(-y_lim, y_lim)

    # label text
    if A_label:
        axes.annotate(A_label, (-a_off - d2, a_off), ha='center')
    if B_label:
        axes.annotate(B_label, (b_off + d2, b_off), ha='center')

    # size text
    for v,x in zip((A-AB, B-AB, AB), (a_mid, b_mid, ab_mid)):
        if v:
            axes.annotate(str(v), (x, 0), ha='center', va='center')

    # set aspect to equal to prevent distortion
    axes.set_aspect('equal')
    # turn off the axis labels, they're useless for a Venn
    axes.set_axis_off()

    if title:
        axes.set_title(title)


def _make_ms1(fig, mz, xy, scan_mode, pm_scanDot, title=None, half_window=4.2):
    '''Plots an MS scan, centered on the target mz.

    This is an internal method, it expects a matplotlib Figure instance'''

    fig.clear()
    fig.set_facecolor('w')

    if scan_mode == 'c':
        plot_xy = []
        for sc in xy:
            # add data point with zeros before and after. the x stays the same,
            # which means the lines are all vertical (looks nicer)
            plot_xy.extend(((sc[0],0), sc, (sc[0],0)))
    else:
        plot_xy = xy

    axes = fig.add_axes([0.125,  0.1,  0.775,  0.8])

    axes.plot([i[0] for i in plot_xy],
              [i[1] for i in plot_xy],
              c='k')

    if pm_scanDot:
        axes.set_xlim((pm_scanDot[0] - half_window, pm_scanDot[0] + half_window))
        try:
            axes.set_ylim((0.0, max(i[1] for i in xy if abs(i[0] - pm_scanDot[0]) <= half_window) * 1.1))
        except ValueError:
            axes.set_ylim((0.0, 1.0))
        axes.axvline(x=pm_scanDot[0], ymin=0, ymax=1, ls='--', c='b')


    if title is None:
        if isinstance(mz, float):
            mz = '%.3f' % mz

        title = '%s m/z' % mz

    if title:
        axes.set_title(title)

    axes.set_xlabel('m/z')
    axes.set_ylabel('Relative Abundance')

    axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

    return (xy,)


def _make_xic(fig, mz, xic_time, xic_int, scan_dot, bin_times, bin_ints, title=None):
    '''Plots a XIC, with labels for the primary MSMS scan and any neighboring scans
    in the same time and mz area.

    This is an internal method, it expects a matplotlib Figure instance'''

    #assert len(scan_dot) == 2
    if len(scan_dot) != 2 and scan_dot != None:
        if len(scan_dot) >= 1 and len(scan_dot[0]) == 2:
            scan_dot = scan_dot[0]
        else:
            scan_dot = None

    fig.clear()
    fig.set_facecolor('w')

    axes = fig.add_axes([0.125,  0.1,  0.775,  0.8])

    axes.plot(xic_time, xic_int, '--rs', linewidth=2, markeredgecolor='k',
              markerfacecolor='g', markersize=5)

    if len(bin_times) > 0:
        axes.plot(bin_times, bin_ints, 'yo', markersize=10)

    # plot the scan dot last, so it stays on top of other markers
    if scan_dot is not None:
        axes.plot([scan_dot[0]], [scan_dot[1]], 'b^', markersize=10)

    if title is None:
        if isinstance(mz, float):
            mz = '%.3f' % mz

        title = '%s m/z' % mz

    if title:
        axes.set_title(title)

    axes.set_xlabel('Time (min)')
    axes.set_ylabel('Abundance')

    axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

    return (zip(xic_time, xic_int) + ([scan_dot] if scan_dot else []),)


def _make_ms2(fig, scan, scan_mode, peptide, labels=None, ion_list = None,
              charge=None, score=None, label_top=50, filter_labels=True,
              pretty_print=False, title=None, y_range = None, x_range = None, tolerance = None, **settings):
    '''Creates an image for the MS/MS scan with theoretical fragments labeled
    and the sequence on top, with observed ions annotated

    scan = [(mz, int)...]
    ion_list = ('b', 'y')

    This is an internal method, it expects a matplotlib Figure instance'''
    # scan_mode is unused; if you want to go and fix all the things it would break,
    # take it out!

    # initialize defaults and update with optional arguments
    _settings = dict(show_theor_mz=True, ms2_mz_figs=2, show_mass_error=False,
                     mass_error_figs=2, mass_error_units='ppm')
    _settings.update(settings)

    allTheIons = ('b', 'y', 'a', 'c', 'x', 'z',
                  'b++', 'y++', 'a++', 'c++', 'x++', 'z++',
                  'a0', 'b0', 'c0', 'x0', 'y0',
                  'a*', 'b*', 'y*',
                  'MH') # Respecting guidelines from docstring of fragment().
    if ion_list == None:
        ion_list = allTheIons

    if peptide:
        clean_re = re.compile(r'-?(\[(.+?)\]|([cdiop]|s\d{1,2}))-?')

        clean_pep = clean_re.sub('', peptide)

        if labels is None: # create theoretical fragments and do some labeling
            # passing in the top 50 intensities for labeling
            labels = generate_labels(sorted(scan, key=lambda mi: mi[1], reverse=True)[:label_top],
                                     peptide, ion_list, charge, tolerance = tolerance, **_settings) # default tolerance is 0.6 Da

        ion_info, ion_coords = translate_labels(labels, len(clean_pep))

        label_dict = defaultdict(list)
        for m,lbl in labels:
            label_dict[m].append(lbl)

        phospho = []
        oxid = []
        silac = []

        if set(ion_list).intersection(('c', 'z+1')):
            instrument = 'ETD-TRAP'
        else:
            instrument = ''

        # in case someone has explicitly named these mods, we need to replace them
        oxphossilac_re = re.compile('\[(o)xidation\]|\[(p)hospho\]|\[(s)ilac\]', flags=re.IGNORECASE)
        # and we'll replace all the other bracket-based mods because they might contain
        # lower-case characters
        oxphossilac_pep = re.sub('\[.+?\]', '', oxphossilac_re.sub(lambda g: (g.group(1) or g.group(2)).lower(),
                                                         peptide))
        
        # Terminal peptides will also still have their -, invalidating mod location!  Fixing that.
        if oxphossilac_pep[0] == '-': oxphossilac_pep = oxphossilac_pep[1:]

        # check for o and p (oxidation and phosphorylation) in the peptide string
        for i,aa in enumerate(oxphossilac_pep):
            if aa == 'o':
                oxid.append(len(clean_re.sub('', oxphossilac_pep[:i])))
            elif aa == 'p':
                phospho.append(len(clean_re.sub('', oxphossilac_pep[:i])))
            elif aa == 's':
                silac.append(len(clean_re.sub('', oxphossilac_pep[:i])))            

        phospho = frozenset(phospho)
        oxid = frozenset(oxid)
        silac = frozenset(silac)

    fig.clear()
    fig.set_facecolor('w')

    x = [x1 for x1,y1 in scan]
    y = [y1 for x1,y1 in scan]

    if peptide:
        axes = fig.add_axes([0.125, 0.1, 0.775, 0.65])
    else:
        axes = fig.add_axes([0.125, 0.1, 0.775, 0.8])

    ## stem plot: unfortunately slower, but more consistent
    markerline, stemlines, baseline = axes.stem([x1 for x1,y1 in scan if y1 > 0],
                                                [y1 for x1,y1 in scan if y1 > 0],
                                                linefmt='k-',
                                                basefmt='k-',
                                                markerfmt='k ')
                            

    if peptide:
        (x_lbl,y_lbl,lbls) = zip(*[(x1,y1,label_dict[x1])
                                   for x1,y1 in scan if x1 in label_dict]) or ((),) * 3
        if x_lbl:
            axes.scatter(x_lbl, y_lbl, s=5, marker='o', linewidth=0.5, facecolors='w', zorder=2)

            for s in stemlines:
                if s.get_xdata()[0] in label_dict:
                    s.set_linewidth(0.7)
                    s.set_color(get_color(label_dict[s.get_xdata()[0]]))
                else:
                    s.set_linewidth(0.5)

    if title is None and (charge or score):
        if isinstance(charge, float):
            charge = '%d' % charge

        if isinstance(score, float):
            score = '%.2f' % score

        # this is deprecated behavior, but I'm not sure how to deal with it gracefully
        title = ''.join(('Charge: %s+' % charge if charge else '',
                         '    ' if charge and score else '',
                         'Score: %s' % score if score else ''))

    if title:
        axes.set_title(title)

    axes.set_xlabel('m/z')
    axes.set_ylabel('Abundance')

    if not peptide:
        axes.set_ylim(ymin=0.0, ymax=max(y) * 1.1)

        axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

        return (None, None)

    # this stuff needs tweaking.
    if pretty_print:
        x_off = -0.01
        y_off = 0.02
        y_chr = 0.023
        last_xy = (-100,0)
        last_lbl = -100

    if filter_labels:
        #lbl_re = re.compile(r'(([byc]|z\+1)\(\d+\)(-\d+)?) .*')
        lbl_re = re.compile(r'(([byc]|z\+1)\(\d+\)(\++)?(-\d+)?).*')

    for x1,y1 in scan:
        if y1 == 0:
            continue

        if pretty_print:
            # translate data points into display pixels
            x2,y2 = axes.transAxes.inverted().transform(axes.transData.transform([x1,y1]))

        if x1 in label_dict:
            if filter_labels == 'sensible':
                lbl = ', '.join([x.split()[0] for x in label_dict[x1]])
            elif filter_labels:
                lbl = ', '.join(lbl_re.match(lbl).group(1) for lbl in label_dict[x1] if lbl_re.match(lbl))
            else:
                lbl = ', '.join(label_dict[x1])

            if not lbl:
                continue

            if pretty_print:
                if ((x2 - last_xy[0]) < 0.017
                    and last_off > y2
                    and y2 + y_off + y_chr * len(lbl) > last_xy[1]):
                    y2 = last_off
                    arrowprops = {'arrowstyle': '-',
                                  'color': 'g',
                                  'linewidth': 0.5}
                else:
                    arrowprops = None

                axes.annotate(lbl, (x1,y1),
                              (x2, y2),
                              #(x2 + x_off, y2 + y_off),
                              arrowprops=arrowprops,
                              textcoords='axes fraction',
                              rotation='vertical',
                              fontsize='x-small',
                              horizontalalignment='center',
                              verticalalignment='bottom')
            else:
                axes.annotate(lbl, (x1,y1),
                              #(x2, y2),
                              (0, 2),
                              textcoords='offset points',
                              rotation='vertical',
                              fontsize='x-small',
                              horizontalalignment='center',
                              verticalalignment='bottom')
        else:
            lbl = ''

        if pretty_print:
            last_xy = (x2,y2)
            last_off = y2 + y_off + y_chr * len(lbl)

    ion_sets = []

    for coord_dict in ion_coords:
        ion_sets.append([(x1,y1,get_color(coord_dict[x1]))
                         for x1,y1 in scan if x1 in coord_dict])

    annotations = dict()
    for i,ion_set in enumerate(ion_sets):
        if not ion_set:
            continue
        c = axes.scatter([x1[0] for x1 in ion_set],
                         [y1[1] for y1 in ion_set],
                         s=20, marker='o', facecolors='none',
                         edgecolors=[c[2] for c in ion_set],
                         linewidths=2.0, zorder=4, visible=False)
        annotations[i] = c

    # set the y-limit to give some space above the top peak for labels
    if y_range:
        axes.set_ylim(ymin = y_range[0], ymax = y_range[1])
    else:
        axes.set_ylim(ymin=0.0, ymax=max(y) * 1.3)
        
    if x_range:
        axes.set_xlim(xmin = x_range[0], xmax = x_range[1])

    axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

    peptide_xes = add_pep_text(fig, clean_pep, ion_matches=ion_info,
                               phospho=phospho, oxid=oxid, silac=silac, instrument=instrument)

    peptide_annotations = [(x,annotations.get(i, None)) for i,x in enumerate(peptide_xes)]

    return zip(x_lbl, y_lbl, lbls), peptide_annotations




#if __name__ == '__main__':
    #print "TEST MODE"
    #from multiplierz.mzAPI import mzFile
    #import matplotlib.pyplot as pyt
    #datafile = r'\\glu2\PipelineStuff\testUser\copdTestCys\2014-10-02-Amgen-COPD-Set2-15-800.raw'
    #data = mzFile(datafile)
    #scan = data.scan(16829)
    ##fig = Figure()
    #foo = make_ms2_im("foobar3.pdf", [x for x in scan if x[0] < 1000], None, '[229.162]-AAAAGLGHPApSPGGpSEDGPPGSEEEDAAREGTPGSPGRGR', 
                    #None, ('b', 'y'), charge = 6, y_range = (0, 10000), x_range = (300, 1000))
    #print "Done!"