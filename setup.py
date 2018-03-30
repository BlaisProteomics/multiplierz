from setuptools import setup, find_packages, dist

from codecs import open
from os import path


import platform
if 'windows' in platform.platform().lower():
    install_requires = ['numpy', 'comtypes', 'matplotlib', 'pypiwin32',
                          'openpyxl', 'xlrd', 'xlwt', 'requests'], # Removed 'lxml'.
else:
    print "Preparing Linux mode installation."
    install_requires = ['numpy', 'comtypes', 'matplotlib',
                        'openpyxl', 'xlrd', 'xlwt', 'requests']

setup(name = 'multiplierz',
      version = '2.0.12',
      description = 'The multiplierz proteomics package',
      author = 'William Max Alexander (et al.)',
      author_email = 'williamM_alexander@dfci.harvard.edu',
      classifiers = ['Development Status :: 4 - Beta',
                     'Intended Audience :: Science/Research',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
                     'Programming Language :: Python :: 2.7',
                     ],
      keywords = 'biology bioinformatics proteomics spectrometry',
      packages = find_packages(),
      package_data = {'multiplierz.mzAPI':['mzAPI/agilentdlls/*',
                                           'mzAPI/t2ddlls/*',
                                           'mzAPI/wiffdlls/*',
                                           'mzAPI/rawdlls/*'],
                      'multiplierz':['unimod.sqlite', '_msparser.pyd',]},
      include_package_data=True,
      install_requires = install_requires
      )



        
