# NMR heteronuclear NOE analysis tools
# Dan Richman, 2013-2014

import pandas as pd
import numpy as np
import nmrfn
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

def noe_analysis():
    """Heteronuclear NOE analysis workflow example; runs as main."""

    files = ['f34k_7p2_noe_heights.txt', 'dphs_7p4_noe_heights.txt']
    proteins = ['f34k','dphs']

    heights = [nmrfn.parse_noe(f) for f in files]
    heights = pd.concat(heights, axis=1, keys=proteins)

    results = [noe_calcs(heights[i]) for i in proteins]
    results = pd.concat(results, axis=1, keys=proteins)

    fig = plt.figure(figsize=(10,6))
    val = plt.subplot(211)
    dif = plt.subplot(212)
    
    seqmin = 6
    seqmax = 144
    noemin = 0.5
    noemax = 0.95

    dy = results.dphs.noe - results.f34k.noe
    dif_plot(val, dy.index, dy.ix[8:141], seqmin, seqmax)
    noe_plot(proteins, results, seqmin, seqmax, noemin, noemax, xax=False)

    s1 = plt.Rectangle((0, 0), 1, 1, fc=proteins['f34k'])
    s2 = plt.Rectangle((0, 0), 1, 1, fc=proteins['dphs'])
    val.legend([s1, s2], ['F34K', '∆+PHS'], loc=8)

    title = 'HN NOE: F34K pH 7.2 compared with ∆+PHS pH 7.4'
    savefile = 'f34k_7p2_noe_analysis.pdf'
    save_plot(title, savefile)

def noe_calcs(df):
    """Calc heteronuclear NOE from single "on" and "off" pair of experiments.
    
    Parameters
    ----------
    df : DataFrame prepared by parse_noe in nmrfn module
         
    Returns
    -------
    DataFrame indexed by residue, columns are NOE value and standard deviation
    """
    noe_status = ['off','on']
    off_mean, on_mean = [df[i].mean(1) for i in noe_status]
    off_sdev, on_sdev = [df[i].std(1) for i in noe_status]
    noe_val = on_mean/off_mean
    noe_err = np.sqrt( (noe_val**2) * ( (on_sdev/on_mean)**2 + \
        (off_sdev/off_mean)**2 ) )
    return pd.concat([noe_val,noe_err], axis=1, keys=['noe','err'])

def noe_plot(dic, df, seqmin, seqmax, noemin, noemax, xax=True):
    """Builds subplot of NOE values from multiple datasets.

    Parameters
    ----------
    dic : dictionary of dataset names and corresponding colors
    df : DataFrame with NOE calculation results
    seqmin, seqmax : integers to control x axis length (residues)
    noemin, noemax : floats to control y axis
    xax : boolean to turn x axis display on (True) or off (False)
    """
    for i in proteins:
        x = df[i].noe.index
        y = df[i].noe
        e = df[i].err
        val.errorbar(x, y, yerr=e, fmt='o', markersize=4, color=dic[i])
    val.set_xlim(seqmin, seqmax)
    val.xaxis.set_major_locator(MaxNLocator(15))
    val.set_ylim(noemin, noemax)
    val.set_ylabel('NOE')
    plt.setp(val.get_xticklabels(), visible=xax)

def dif_plot(sp, x, y, seqmin, seqmax):
    """Builds subplot of differences between NOE values from two datasets.
    
    Parameters
    ----------
    sp : plt.subplot object
    x,y : DataFrame.index and DataFrame objects, or arrays
    seqmin, seqmax : integers to control x axis length (residues)
    """
    sp.bar(x, y, color='gray')
    sp.set_xlim(seqmin, seqmax)
    sp.xaxis.set_major_locator(MaxNLocator(15))
    sp.set_xlabel('Residue')
    sp.set_ylabel('Difference')
    sp.grid(False)
    sp.axhline(y=0, linewidth=1, color='black')

def save_plot(title, savefile):
    fig.suptitle(title)
    fig.savefig(savefile)

if __name__ == '__main__':
    noe_analysis()