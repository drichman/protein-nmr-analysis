# Module with functions relevant to data parsing and analysis
# Dan Richman, 2014 March

import pandas as pd
import numpy as np


# Parsers and their dependencies    
    
def parse_rd(f):
    """Returns DataFrame indexed by timepoints as floats (in seconds)"""
    data = sparky_rh(f).T
    time = data.index.to_series().str.match('([0-9]\.[0-9]+).*').str.get(0)
    data.index = time.astype('float')
    new_ind = []
    for i in data.index:
        if i == 0:
            new_ind.append(i)
        else:
            new_ind.append(1/(4*i*1e-3)) #  Convert to "CPMG freq" (Hz)
    data.index = new_ind
    return data.sort_index()
    
def parse_hx(f, t0='manual'):
    """Returns DataFrame indexed by datetime timepoints
    
    t0='manual' if timepoints in seconds are already encoded in Sparky rh file       
    columns; otherwise should be a timestamp formatted as '%Y-%m-%d_%H:%M'
    """
    data = sparky_rh(f).T
    if t0 == 'manual':
        data.index = manual_timepoints(data.index)
    else:
        data.index = auto_timepoints(data.index, t0)
    return data

def manual_timepoints(time_list):
    """Parses timepoints encoded as seconds in Sparky rh column names"""
    time_strings = time_list.to_series().str.match('([0-9]+)\/.*').str.get(0)
    return time_strings.values.astype('float')
    # () first set of match instructions, [0-9] contains digits character class.
    # + repeats character class at least once, \ ensures / is matched literally.
    # . matches anything, * repeats preceding item zero or more times.

def auto_timepoints(time_list, t0):
    """Parses timepoints encoded as dates/times in Sparky rh column names"""
    time_strings = [i.split('/')[0] for i in time_list]
    times = pd.to_datetime(time_strings, format='%Y-%m-%d_%H:%M')
    t0 = pd.to_datetime(t0, format='%Y-%m-%d_%H:%M')
    delta_t = [i-t0 for i in times]
    return [i.total_seconds() for i in delta_t]

def parse_noe(f):
    """Returns DataFrame with MultiIndex columns to reflect NOE off or on"""
    data = sparky_rh(f)
    cols = [i.split('/')[0] for i in data.columns]
    multi_cols = [ [i.rstrip('0123456789') for i in cols], cols]
    data.columns = pd.MultiIndex.from_arrays(multi_cols)
    return data.sort_index()

aa_chars = 'AEGHIKLMNPRSTUVY-?'

def sparky_rh(f):
    """Basic parser of Sparky rh file into DataFrame.
    
    Parameters
    ----------
    f : string filename, often passed from parse_hx or parse_rd call
    
    Returns
    -------
    DataFrame indexed by residue; columns are some identifier eg timepoints
    """
    data = pd.read_table(f, delim_whitespace=True,
        index_col='Assignment').iloc[:,1:].drop(['SD'], axis=1)
    int_index = [int(i.strip(aa_chars)) for i in data.index]
    data.index = int_index
    return data
                
def sparky_lt(f):
    """Read Sparky heights list; return Series of heights indexed by residues.
    
    Parameters
    ----------
    f : string filename
    
    Returns
    -------
    Series of heights as integers, indexed by integer amino acid labels
    """
    data = pd.read_table(f, delim_whitespace=True, squeeze=True,
        index_col='Assignment').ix[1:,1:2]
    int_index=[int(i.strip(aa_chars)) for i in data.index]
    data.index = int_index
    return data

def sphere_file(f):
    """Determine the reference sphere file for a given pH and parse it.
    
    Parameters
    ----------
    f : filename
    
    Returns
    -------
    DataFrame with residues and reference k values, indexed by residues
    """
    return pd.read_fwf(f, colspecs=[(0, 3), (8, 17)], index_col=0, 
        squeeze=True)


# Calculation functions

def decay(t, a, k):
    return a*np.exp(-k*t)
