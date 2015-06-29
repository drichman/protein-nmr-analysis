# NMR relaxation dispersion analysis functions
# Dan Richman, 2013-2014

import pandas as pd
import numpy as np
import nmrfn
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from os import mkdir, system

def rd_analysis(): 
    """For scripting purposes this function should be reorganized.
    Demonstrates workflow replaced by IPython Notebooks. Runs as main.
    
    Returns
    -------
    Happy scientists
    """
    
    curve_files = ['L25K_6p3_rd_600_full_heights.txt',
        'L25K_6p3_rd_800_full_heights.txt']
    endpt_files = ['l25k_6p3_600_rd_l9=2.txt','l25k_6p3_600_rd_l9=32.txt',
        'l25k_6p3_800_rd_l9=2.txt','l25k_6p3_800_rd_l9=15.txt']

    # Total R2 relaxation period (ms):
    CT_B1 = 38.4
    CT_B2 = 28.8
        
    # Load data; 2 field strengths, 2 data types (full curves, endpoints):
    # Full dispersion curve peak intensities:
    B1curve, B2curve = [nmrfn.parse_rd(f) for f in curve_files]
    B1curve['index']=B1curve.index # Eliminating triplicate measurement...
    B1curve = B1curve.drop_duplicates(cols='index') #...
    B1curve = B1curve.drop('index', axis=1) #...done.

    # Curve endpoints peak intensities:
    B1few, B1many, B2few, B2many = [nmrfn.parse_lt(f) for f in endpt_files]
    
    # Compute Rex estimates from endpoints data
    rexB1endpts = est_endpts(B1few, B1many, CT_B1)
    rexB2endpts = est_endpts(B2few, B2many, CT_B2)

    # Convert full dispersion curve peak intensities to R2 values
    B1curveR2 = calc_R2series(B1curve, CT_B1)
    B2curveR2 = calc_R2series(B2curve, CT_B2)
    
    # Using highest-frequency CPMG R2 as intrinsic R2:
    B1intR2 = B1curveR2.iloc[0,:]
    B2intR2 = B2curveR2.iloc[0,:]
    
    # Compute Rex estimates from fitting the full R2 curves:
    kex0 = 1e3
    Rex0 = 1.0
    rexB1fit = fit_curves(B1curveR2, B1intR2, kex0, Rex0)
    rexB2fit = fit_curves(B2curveR2, B2intR2, kex0, Rex0)

    rex_output(rexB1fit, rexB2fit, rexB1endpts, rexB2endpts)
    alpha_output(rexB1fit, rexB2fit, rexB1endpts, rexB2endpts)
           
    plot_r2_cpmg(B1curveR2, 'R2_series_600_MHz')
    plot_r2_cpmg(B2curveR2, 'R2_series_800_MHz')
    
    #bin_rex(rex_vs_seq, output='rex_est_bins.txt')    

def rex_output(rexB1fit, rexB2fit, rexB1endpts, rexB2endpts):
    rex = pd.concat([rexB1fit, rexB2fit, rexB1endpts, rexB2endpts], axis=1)
    rex.columns = ['Rx1','Rx1er','Rx2','Rx2er','Rx1ep','Rx2ep']
    rex.to_csv('rex.txt', sep='\t', float_format='%.2f')

def alpha_output(rexB1fit, rexB2fit, rexB1endpts, rexB2endpts):
    alpha_fits = alpha(rexB1fit, rexB2fit)
    alpha_endpts = alpha(rexB1endpts, rexB2endpts)
    alpha_both = pd.concat([alpha_fits, alpha_endpts], axis=1)
    alpha_out = alpha_both.loc[:,['Rex','Height']]
    alpha_out.columns = ['fits_a','endpts_a']
    alpha_out.to_csv('alpha.txt', sep='\t', float_format='%.2f')

def alpha(rexB1, rexB2, B1=14.1, B2=18.8):
    """
    Parameters
    ----------
    rexB1, rexB2 : series containing Rex values indexed by residues
    B1, B2 : floats, magnetic fields in Tesla, B1 < B2
    
    Returns
    -------
    Series of alpha values indexed by residues
    """
    B_ratio = (B2+B1)/(B2-B1)
    alpha = B_ratio*(rexB2 - rexB1)/(rexB2 + rexB1)
    return alpha.sort_index()

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
    rex_est = -np.log(ratio)/(ct_time*1e-3)
    return rex_est

def fit_curves(R2data,R2int,kex0,Rex0):
    """Perform dispersion fitting for whole set of residues.
    
    Parameters
    ----------
    R2data : DataFrame with R2 vs tau_cp (as indices), columns are residues
    
    Returns
    -------
    fits : DataFrame indexed by residues, columns: Rex, Rex error
    """
    fits=[fit_disp(R2data[res], R2int.ix[res], kex0, Rex0) for res in R2data]
    fits = pd.concat(fits, axis=1, keys=R2data.columns).T
    return fits

def fit_disp(s, r, kex0, Rex0):
    """Fit dispersion function to R2 vs tau_cp data and extract Rex."""
    tcp = s.index.values.astype('float')
    y = s-r
    try:
        popt, pcov = curve_fit(disp, xdata=tcp, ydata=y, p0=[kex0,Rex0])
        Rex = popt[1]
        Rex_er = np.sqrt(pcov.diagonal()[1])
    except RuntimeError:
        print("Error in curve_fit, residue ", str(s.name))
        Rex = np.nan
        Rex_er = np.nan
    return pd.Series([Rex, Rex_er], index=['Rex','Rex_er'])

def disp(tcp, kex, Rex):
    return Rex*(1 - (2/(kex*tcp))*np.tanh(kex*tcp/2))

def calc_R2series(data, ct_time):
    """Calculate the R2 series for each residue.
    
    Parameters
    ----------
    data : DataFrame of peak heights indexed by tau_cp (ms), cols are residues
    ct_time : float representing constant-time relaxation period in ms
    
    Returns
    -------
    DataFrame of R2 values indexed by CPMG frequency, each residue as a column
    """
    data_series = data.iloc[1:,:]
    data_I0 = data.iloc[0,:]
    calcs = [calc_r2(data_series[res], data_I0[res], ct_time*1e-3) 
        for res in data_series]
    r2_series = pd.concat(calcs, axis=1, keys=data.columns)
    return r2_series

def calc_r2(s, I0_val, ct_time_in_s):
    """Calculate the R2 values in a relaxation dispersion series.
    
    Parameters
    ----------
    s : Series
    ct_time_in_s : float representing constant-time relaxation period (seconds)
    
    Returns
    -------
    Series containing R2 values indexed by CPMG periods in ms
    """
    new = pd.Series(index=s.index) # new is so that s doesn't get overwritten
                                   # ("Returning a view vs a copy" in pandas)
    for i in s.index:
        new.ix[i] = np.log(I0_val/s.ix[i])/ct_time_in_s
    return new

def plot_r2_cpmg(R2data, directory):
    """Plot R2 vs CPMG frequency data to inspect."""
    R2data = R2data.sort_index(axis=1)
    for res in R2data:
        plt.figure()
        ax = R2data[res].plot(style='o')
        ax.set_ylabel('$R_2$ (1/s)')
        ax.set_xlabel('CPMG freq (Hz)')
        plt.savefig(str(res)+'.png')
    mkdir(directory)
    system('mv *.png '+directory)

def plot_rex(data, output):
    """Plot Rex estimates vs sequence.
    
    Parameters
    ----------
    heights_file : string directory name
    output : string file name with suffix eg .pdf

    Returns
    -------
    None
    """
    fig, rex = plt.subplots(figsize=(10, 5))
    rex.bar(data.index, data, align='center', color='gray')
    rex.set_ylabel('$R_{ex}$ estimated (1/s)')
    rex.set_xlabel('Residue')
    fig.savefig(output)

def bin_rex(data, output):
    bins = data.groupby(pd.cut(data, bins=[-2.5,5,10,15,25]))
    with open(output, 'w') as fw:
        for b, g in bins:
            fw.write(str(b)+'\n')
            fw.write(str(g)+'\n')
    #bins.sort_index().to_csv(output, sep='\t')

def r2_cpmg_inspection():
    """Tool to visualize R2 vs CPMG freq computed through Ananya's workflow.
       Rarely used (nmr-analysis/rd_analysis.py replaced Ananya's workflow).
       But this can be executed if necessary.
       
       Because of its rare use, this function contains imports only it uses:"""
    from os import listdir, getcwd
    from matplotlib.backends.backend_pdf import PdfPages

    files = [i for i in listdir(getcwd()) if i.endswith('_R2.out')]

    fig = plt.figure(figsize=(8.5,11))
    #fig.suptitle('R_2 vs CPMG freq, DPHS L25K pH 6.3, 25C',fontsize=12)
    fz = 8

    for i, filename in enumerate(sorted(files, \
        key = lambda x: int(x.split('_')[0].translate(None, \
        "AaDdEeGgHhIiKkLlMmNnoPpQRrSsTtUuVvYy")))):
        resname = filename.split('_')[0]
        resdata = resname
        resdata =pd.read_table(filename,sep='   ',names=['nu','r2'],header=None)
        sub = fig.add_subplot(8,4,i+1)
        sub.plot(resdata.nu,resdata.r2,'bo')
        plt.setp(sub.get_xticklabels(),visible=False)
        plt.setp(sub.get_yticklabels(),fontsize=fz)
        sub.text(0.85,0.85,resname,\
            horizontalalignment='right',verticalalignment='top',\
            transform=sub.transAxes,fontsize=fz)
        if sub.is_last_row():
            sub.set_xlabel('CPMG freq (1/s)',fontsize=fz)
            pl.setp(sub.get_xticklabels(),fontsize=fz,visible=True,rotation=45)
        if sub.is_first_col():
            sub.set_ylabel('$R_2$ (1/s)',fontsize=fz)
    
    plt.tight_layout(pad = 6.0, h_pad = 0.1, w_pad = 0.0)
    
    plotfile = PdfPages('R2_vs_CPMG_plots.pdf')
    plotfile.savefig(fig)
    plotfile.close()
        
if __name__ == '__main__':
    rd_analysis()