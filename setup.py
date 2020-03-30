from setuptools import setup, find_packages, dist

from codecs import open
from os import path




install_requires = ['numpy', 'matplotlib', 'openpyxl', 'xlrd',
                    'xlwt', 'requests', 'beautifulsoup4', 'pythonnet',
                    'pypiwin32; platform_system == "Windows"',
                    'comtypes; platform_system == "Windows"']



README = """ 
**multiplierz** is a Python software library and associated GUI \
desktop environment for managing proteomic mass spectrometry workflows and \
data analysis. Using the mzAPI interface to native instrument data formats, \
multiplierz is provides a complete toolset for a variety of methods for \
peptide identification, quantitation, and experimental reporting.

More information can be found on [the multiplierz Github \
page](https://github.com/BlaisProteomics/multiplierz) .

"""

setup(name = 'multiplierz',
      version = '2.2.1',
      description = 'The MultiplierZ proteomics package',
      long_description = README,
      author = 'William Max Alexander (et al.)',
      author_email = 'williamM_alexander@dfci.harvard.edu',
      classifiers = ['Development Status :: 4 - Beta',
                     'Intended Audience :: Science/Research',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
                     'Programming Language :: Python',
                     ],
      keywords = 'biology bioinformatics proteomics spectrometry',
      packages = find_packages(),
      package_data = {'multiplierz.mzAPI':['agilentdlls/*',
                                           't2ddlls/*',
                                           'wiffdlls/*',
                                           'rawdlls/*',
                                           'brukerlib/*'],
                      'multiplierz':['unimod.sqlite', '_msparser.pyd',]},
      include_package_data=True,
      install_requires = install_requires
      )



        
