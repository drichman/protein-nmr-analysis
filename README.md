## Overview

A collection of tools to help process and analyze protein NMR data. Archives code from my thesis work and papers in the Garcia-Moreno Lab at Johns Hopkins 2011-2015. Not a cohesive Python package.

## Contents

### modules

Import or copy into data analysis notebooks or scripts

- **nmrfn** Functions for parsing Sparky files from different types of experiments and some ancillary functions
- **hdx_analysis** Functions for analyzing hydrogen-deuterium exchange 
- **rd_analysis** Functions for analyzing relaxation dispersion
- **noe_plot** Functions for calculating and plotting heteronuclear NOE results

### notebooks

IPython Notebooks from a subset of the analysis tasks performed in my work; for the best example of a notebook that documents and automates the quantitative work of a paper, see 2014-12-06_l25k_l125k_analyses_and_figs.ipynb

### scripts

Specific data analyses:

- (folder) **hdx-example** Self-contained example (script, module, and input data) of HDX decay fitting and analysis; good for pedagogical purposes
- (folder) **rd-example** Self-contained example (script and input data) of estimating the conformational exchange contribution to relaxation dispersion (RD) from the low- and high-CPMG frequency endpoints of RD
- (folder) **deprecated-proteins-2014** Early scripts to visualize relaxation data before I organized code in functions or used IPython Notebooks; archival for Richman et al. *Proteins* 2014

Assorted tools:

- **nmrpipe_batchprocess.bash** For batch-processing (namely running fid.com and ft2.com repeatedly) through a series of spectra such as relaxation timepoints
- **guardd_formats.py** Converts a Sparky "rh" table to the input format for GUARDD
- **nmr_clean.py** Cleans non-essential Bruker files and processed data from directories
- **transport.bash** A way to fetch a script from a central directory, eg if you store a template fid.com somewhere (but this is not strictly NMR related)