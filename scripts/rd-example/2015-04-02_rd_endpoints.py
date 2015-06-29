# Python script for estimating Rex from assigned RD endpoints spectra
# Dan Richman, 2015 April
"""
Requires Sparky peak lists ("lt" type) from the "many loops" and "few loops" 
spectra (columns turned on in Sparky: Assignment, w1, w2, Data Height; see
included example files for more detail).

Just a reminder:
many loops : high CPMG freq : much refocusing : small conf exch decay (low Rex)
few loops : low CPMG freq  : little refocusing : big conf exch decay (high Rex)

Output: plot of Rex estimates vs sequence
Usage: python 2015-04-02_rd_endpoints.py

(Depends on Python 3, pandas, numpy, and matplotlib.)
"""

import pandas as pd
from numpy import log
import matplotlib.pyplot as plt

##### Definitions of functions that are used below

def est_endpts(few_loops, many_loops, ct_time):
    """Calculate the intensity ratio of CPMG extremes and make Rex estimate.
    
    Parameters
    ----------
    few_loops : Series with peak height data vs residue
    many_loops : Series with peak height data vs residue
    ct_time : float representing constant-time relaxation period in ms
    
    Returns
    -------
    Series containing Rex estimates indexed by residues
    """
    ratio = few_loops/many_loops #low_peak_intensity/high_peak_intensity
    rex_est = -log(ratio)/(ct_time*1e-3)
    return rex_est

def sparky_lt(f):
    """Read Sparky "lt" list, return parameters indexed by residues.
    
    Parameters
    ----------
    f : string filename
    
    Returns
    -------
    Series of heights as integers, indexed by amino acid labels
    """
    data = pd.read_table(f, delim_whitespace=True, squeeze=True,
        index_col='Assignment').drop('Height',axis=1).dropna()
    if '?-?' in data.index:
        data.drop('?-?', inplace=True)
    resnames = [i.rstrip('NH-') for i in data.index]
    data.index = resnames
    data.columns=['wN','wH','height']
    return data

ab_up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
ab_lo = ab_up.lower()

def int_index(df):
    ints = [int(i.strip(ab_up+ab_lo+'-')) for i in df.index]
    df.index = ints
    return df.sort_index()

##### End function defs

##### Script workflow follows
# Adjust filenames, CT R2 decay length, and plot parameters as necessary

many = int_index(sparky_lt('l25k_600_1008.06.list')['height'])
few = int_index(sparky_lt('l25k_600_34.92.list')['height'])

ct = 28.8 # constant-time R2 decay period in milliseconds

rex = est_endpts(few, many, ct)

plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size']= 9

rex.plot(kind='bar', figsize=(18,5));

plt.ylim(ymax = 25, ymin = -2) # Change these limits to suit Rex value range,
# or comment out this line to restore autoscaling y axis.

plt.xticks(rotation=75)
plt.xlabel('Residue');
plt.ylabel('$R_{ex}$ ($s^{-1}$)');

plt.grid(False)

plt.savefig('rex_estimates.png')