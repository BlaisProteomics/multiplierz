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


__author__ = 'Jignesh Parikh, James Webber, William Max Alexander'
__version__ = '2.0.0'

__all__ = ['mzAPI', 'mzTools', 'mzReport',
           'myHome', 'myData', 'logger_message', 'SettingsFile']

import logging
import os, sys



# Based on matplotlib's _get_home function (in matplotlib/__init__.py)
def _get_home():
    """
    Find user's home directory if possible. Otherwise raise an error.
    """

    path = ''

    try:
        path=os.path.expanduser("~")
    except:
        pass

    if not os.path.isdir(path):
        for evar in ('HOME', 'USERPROFILE', 'TMP'):
            try:
                path = os.environ[evar]
                if os.path.isdir(path):
                    break
            except:
                pass

    if path:
        return path
    else:
        raise RuntimeError('No home environment variable found')


# these directories should be cross-platform:

# user's home directory
myHome = _get_home()

# multiplierz data folder
myData = os.path.join(myHome, ".multiplierz")
myTemp = os.path.join(myData, "TEMP")
SettingsFile = os.path.join(myData, 'settings.txt')
modFile = os.path.join(myData, 'mods.txt')

verbose_mode = True
def vprint(thing):
    """
    Multiplierz verbosity-sensitive print function, for use 
    internally.  Prints only if multiplierz.verbose_mode == True.
    """
    global verbose_mode # Slightly speeds up lookup, supposedly.
    if verbose_mode:
        print thing
        


class legacy_logging(object):
    # The old logging module assumed everything was in the old mz-scripting
    # box, and stomped over all the logger-settings.
    def __init__(self):
        #self.logger = logging.getLogger('multiplierz_legacy_logger')
        pass
    
    def __call__(self, level = 30, message = 'Foo'):
        vprint(message)
        #self.logger.log(level, str(message))
        
    
    def set_level(self, level):
        pass

logger_message = legacy_logging()

def get_module_logger(module_name):
    log = logging.getLogger(module_name)
    log.addHandler(logging.NullHandler())
    return log



def load_mods():
    """
    A list of modifications for use in Comet, XTandem, and etc.
    
    List is (modification, specificity, mass) tuples.
    """
    
    assert os.path.exists(modFile), 'mods.txt not found!'
    
    mods = []
    with open(modFile, 'r') as file:
        for line in file:
            if not line.strip() or line.strip()[0] == '#':
                continue
            
            mod, sites, mass = line.split()
            mods.append((mod, sites, float(mass)))
    
    return mods
    
fastaFiles = os.path.join(myData, 'fastafiles.txt')
def fastaList():
    fastas = []
    with open(fastaFiles, 'r') as file:
        for line in file:
            fastas.append(line.strip())
    return fastas
        




# The first time multiplierz is loaded, it will set up the
# userDirectory/.multiplierz directory, and populate it with
# some data files.  These are in particular the data that
# users are supposed to be able to edit or update, such as
# the settings file or list of peptide modifications.

def initialSettings():
    # Can't be in settings because that module can't be imported until
    # a settings file exists.
    defaultSettingsFile = """
    # Multiplierz Settings File
    # Lines beginning with a '#' are considered comments.
    
    # This is used only for HTTP requests made by pep2gene, and not required.
    user email=not_set
    
    
    logger verbosity=30
    format default=.xls
    image width=8.0
    image height=6.0
    
    # peak viewer
    XIC_view time=0.5
    XIC_gen time=1.0
    XIC_gen mz=0.02
    MS1_view mz=4.0
    
    # report viewer
    sig_figs ms1_mz=2
    sig_figs ms1_int=0
    sig_figs ms2_mz=2
    sig_figs ms2_int=0
    sig_figs xic_time=1
    sig_figs xic_int=1
    ion_labels theor=1
    ion_labels error=0
    ion_labels figs=2
    ion_labels units=ppm
    
    mzServer use=ask
    
    # Comet
    comet directory=C:\Path\To\Comet.exe
    
    # XTandem
    xtandem directory=C:\Path\To\XTandem.exe
    
    # Mascot
    mascot server=Undefined
    mascot version=2.4
    mascot ms2=True
    mascot security=False
    mascot var_mods=True
    mascot max hits=100000
    mascot ion cutoff=5
    mascot bold red=False
    mascot show input query=True
    mascot rank one only=True
    mascot pep quant=False
    """    
    ptr = open(os.path.join(myData, 'settings.txt'), 'w')
    ptr.write(defaultSettingsFile)
    ptr.close()
    
def deployUnimod():
    import shutil
    unimodFile = os.path.join(os.path.dirname(__file__), 'unimod.sqlite')
    if not os.path.exists(unimodFile):
        unimodFile = None
        basedir = os.path.dirname(sys.executable)
        for files, subdirs, path in os.walk(basedir):
            if 'unimod.sqlite' in files:
                unimodFile = os.path.join(path, 'unimod.sqlite')
                break
            
    if not unimodFile:
        print "WARNING: No unimod.sqlite found in %s (%s)" % (basedir, sys.executable)    
    else:
        shutil.copy(unimodFile, os.path.join(myData, 'unimod.sqlite'))    
    

def initialMods():
    DefaultModificationFile = """
    # Format is (Modificiation Name) (AA Site Specificity) (Mass)

    Oxidation	M	15.9949
    Phosphorylation	STY	79.966331
    """   
    mods = open(os.path.join(myData,'mods.txt'), 'w')
    mods.write(DefaultModificationFile)
    mods.close()    
    
def initializeFasta():
    fastas = open(os.path.join(myData,'fastafiles.txt'), 'w')
    fastas.close()
    
    

if not os.path.exists(myData):
    print "Multiplierz directory (%s) not found; creating it." % myData
    os.mkdir(myData)
    

requiredFiles = [('settings.txt', initialSettings),
                 ('unimod.sqlite', deployUnimod), 
                 ('mods.txt', initialMods), 
                 ('fastafiles.txt', initializeFasta)]
requiredSubdirs = ['TEMP', # Still used by old mzReport imaging functions.
                   'pyCometDatabases'] # Probably obsolete from fastafiles.txt?


for subdir in requiredSubdirs:
    if not os.path.exists(os.path.join(myData, subdir)):
        print "Required multiplierz data directory %s not found!  Creating it." % subdir
        os.mkdir(os.path.join(myData, subdir))
for filename, initializer in requiredFiles:
    if not os.path.exists(os.path.join(myData, filename)):
        print "Required multiplierz data file %s not found!  Creating it." % filename
        initializer()

from multiplierz.settings import settings
protonMass = 1.0072764 # Used at various points.
