from setuptools import setup, find_packages, dist

from codecs import open
from os import path


setup(name = 'multiplierz',
      version = '2.0.1',
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
                                           'mzAPI/wiffdlls/*'],
                      'multiplierz':['unimod.sqlite', '_msparser.pyd', 'mzAPI/MSFileReader_x86_x64_v3.0SP3.exe']},
      include_package_data=True,
      install_requires = ['numpy', 'comtypes', 'matplotlib', 'pypiwin32',
                          'openpyxl', 'xlrd', 'xlwt', 'requests'], # Removed 'lxml'.
      )



        