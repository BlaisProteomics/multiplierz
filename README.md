# Multiplierz
The *Multiplierz* Proteomics Library


*multiplierz* is a Python software library and associated GUI desktop environment for managing proteomic mass spectrometry workflows and data analysis.  Using the mzAPI interface to native instrument data formats, *multiplierz* is provides a complete toolset for a variety of methods for peptide identification, quantitation, and experimental reporting.

See [the wiki for installation instructions and documentation](https://github.com/MaxAlex/multiplierz/wiki/Installation).

The mzDesktop GUI application can be found at [our Sourceforge page](https://sourceforge.net/projects/multiplierz/).

***

### Key Features

#### Native, Open-Source Python Codebase

*multiplierz* is written for and in Python 2.7, a mature and user-friendly language that is becoming a standard for scientific data analysis, and is fully open-source.  This allows researchers with programming experience full freedom in customizing and improving the *multiplierz* library for individual, novel use cases.  We expect *multiplierz* capabilities to grow alongside the rapid developments in the field of computational proteomics.

#### mzAPI Data Interface

Typical proteomic bioinformtics workflows use tools such as [ProteoWizard](http://proteowizard.sourceforge.net/) to first convert native instrument data into a text-based XML format such as mzML or mzXML; for high-throughput experiments, this entails substantial demands in terms of data processing, storage, and organization.  Instead, *multiplierz* offers fully transparent access to native instrument data through mzAPI, a standardized interface to Thermo, AB Sciex, and Agilent machine data formats, as well as the universal mzML format.

#### mzSearch Database Search Interface

A key step in any proteomic workflow is peptide/protein identification, a sophisticated algorithmic problem solved by applications such as Mascot, Comet, and X! Tandem.  To support routine bioinformatic workflows it is necessary to automate the operation of these programs, so *multiplierz* offers rich, programmatic access to the applications mentioned, with full support for the wide array of parameters requested by each, so that searches can be customized and initiated through elegant Python code.

#### mzDesktop GUI Interface

Not all biologists are computer scientists, so *multiplierz* provides a fully-graphical interface to most of its capabilities, including management of MGF and FASTA files, protein coverage visualization, and database search coordination.


#### Chromatographic Feature Detection, Quantitation Analysis, Intact Protein Charge Deconvolution, and more...

***

### News

* 8/22/2017: The paper on *multiplierz*'s 2.0 release has been [published in Proteomics](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201700091/full).

* 8/01/2017: [mzStudio](https://github.com/BlaisProteomics/mzStudio), an interactive proteomics data browser built using *multiplierz*, has been [published in Proteomes](http://www.mdpi.com/2227-7382/5/3/20/html).

* 5/22/2017: The 2.0 release edition of *multiplierz* is pending publication!  We will update here as features are added, bugs are fixed, etc.

***

### Citation

If you use *multiplierz* in an academic publication, please cite us:

**Alexander, William M., et al. "multiplierz v2. 0: a Python‐based ecosystem for shared access and analysis of native mass spectrometry data." Proteomics (2017).**

***

### Contact the Author

Questions related to use and modification of the *multiplierz* library should be referred to W. Max Alexander at williamM_alexander@dfci.harvard.edu.
