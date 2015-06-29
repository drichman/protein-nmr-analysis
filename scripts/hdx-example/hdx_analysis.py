# NMR hydrogen-deuterium exchange analysis demonstration script
# Dan Richman, 2013-2014

from os import listdir, getcwd
import pandas as pd
import numpy as np
import nmrfn
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def analysis_scheme(peak_heights_files=None):
    """What you want to do. Runs as main.
    
    Parameters
    ----------
    peak_heights_files : None by default, enables reading files from directory
    
    Returns
    =======
    Happy scientists.
    """
    
    if peak_heights_files is None:
        peak_heights_files = \
            [i for i in listdir(getcwd()) if i.endswith('_heights.txt')]
    t0_list = ['2013-12-26_16:38', '2013-12-27_15:45', '2013-12-28_10:29']
    results = multi_pH(dict(zip(peak_heights_files, t0_list))).sort_index()

    calc_slopes(results) # also includes plot_slopes
    bin_dg(results)
    formatter(results).to_csv('hdx_rates.txt', sep='\t', float_format='%.2e')
    
    #pH_vals = ['ph3p9','ph4p5','ph6p7']
    #plot_fits(results, pH_vals)
    #plot_pcterr(results, pH_vals, avg, err)?

def multi_pH(files):
    fits = [single_pH(filename, t0) for filename, t0 in files.items()]
    pH = list(map(extract_pH, files)) # as strings, remember
    return pd.concat(fits, keys=pH)

def single_pH(filename, t0, ref_filename=None):
    """Perform fitting of Series in one file and calculate derived quantities.
    
    Parameters
    ----------
    filename : file with peak height vs time data at one pH*
    ref_filename : file with SPHERE reference exchange rates at one pH
            
    Returns
    -------
    DataFrame indexed by residue, with decay rate, error, 
        and derived quantities for each residue
    """
    data = nmrfn.parse_hx(filename, t0)
    fits = [fit_decay(data[res], a0=4e7, k0=1e-6) for res in data]
    fits = pd.concat(fits, axis=1, keys=data.columns).T

    if ref_filename is None:
        # Make a handy default.
        ref_filename = extract_pH(filename) + '_sphere'
    ref = nmrfn.sphere_file(ref_filename)

    R = 1.986e-3 # Universal gas constant in kcal K^(-1) mol^(-1)
    T = 298.0 # Temperature in K

    fits['dg'] = -R*T*np.log(fits['k_obs']/ref) # Log here is natural log
    fits['dger'] = R*T*fits['k_error']/fits['k_obs']
    #fits['lpf'] = np.log10(ref/fits['k_obs']) # Log_10(protectionfactor)
    return fits

from re import findall
extract_pH = lambda x: findall('[0-9]p[0-9]', x)[0]

def fit_decay(s, a0, k0):
    """Fit a Series to a decaying exponential.
    
    Parameters
    ----------
    s : Series which may contain NaNs
    a0 : guess
    k0 : guess
    
    Returns
    -------
    Series containing k_obs and k_error
    """
    h = s.dropna()
    y = np.asarray(h.values, dtype=float)
    t = np.asarray(h.index.values, dtype=float)
    index = ['dec_frac', 'k_obs', 'k_error', 'er_frac']
    if h.count() < 3:
        print('Res ', s.name, ' has < 3 points, so no fit')
        return pd.Series([np.nan, np.nan, np.nan, np.nan], index=index)
    dec_frac = h.tail().mean()/h.head().mean()
    popt, pcov = curve_fit(nmrfn.decay, xdata = t, ydata = y, p0 = [a0, k0])
    k_obs = popt[1]
    k_error = np.sqrt(pcov.diagonal()[1])
    er_frac = k_error/k_obs
    return pd.Series([dec_frac, k_obs, k_error, er_frac], index=index)

def calc_slopes(fits):
    slope = (fits.ix['4p5']['dg'] - fits.ix['3p9']['dg'])/(4.5-3.9)
    err = (fits.ix['4p5']['dger'] + fits.ix['3p9']['dger'])/(4.5-3.9)
    slope.to_csv('hdx_slopevals.txt', sep='\t', float_format='%.2f')
    plot_slopes(slope, err)
    
    slopebins = slope.groupby(pd.cut(slope, bins=[-0.5,0.5,1,2,5]))
    with open('hdx_dg_slope_bins.txt', 'w') as fw:
        for b, g in slopebins:
            fw.write(str(b)+'\n')
            fw.write(str(g)+'\n')

def bin_dg(fits):
    dgtable = fits['dg'].unstack(level=0)
    dgtable.to_csv('hdx_dg.txt', sep='\t', float_format='%.2f')
    
    dgbins = fits['dg'].groupby(pd.cut(fits['dg'], bins = [0,3,6,9,12]))
    dgbins_unstack = dgbins.apply(lambda x: x.unstack(level=0))
    dgbins_unstack.to_csv('hdx_dg_bins.txt', sep='\t', float_format='%.2f')

def formatter(fits):
    #    fits = fits[['dec_ratio','er_frac','dg','dger','k_obs','k_error']]
    fits['dec_frac'] = fits['dec_frac'].map(lambda x: '%.2f' % x)
    fits['er_frac'] = fits['er_frac'].map(lambda x: '%.2f' % x)
    fits['dg'] = fits['dg'].map(lambda x: '%.2f' % x)
    fits['dger'] = fits['dger'].map(lambda x: '%.2f' % x)
    return fits

def plot_fits(df, pH_vals, output='hdx_rates_plot.pdf'):
    fig, axarr = plt.subplots(len(pH_vals), sharex=True, figsize=(8.5,11))
    ymin = 1e-9
    ymax = 1e-2
    xmin = 10
    xmax = 143
    for i,n in enumerate(axarr):
        a = pH_vals[i]
        n.bar(df.ix[a]['k_obs'].index, df.ix[a]['k_obs'],
                     yerr=df.ix[a]['k_error'], log=1, align='center', color='gray')
        n.set_ylim(ymin, ymax)
        n.set_xlim(xmin, xmax)
        
        n.xaxis.set_major_locator(MaxNLocator(15))
        
        n.set_ylabel('%s $k_{obs}$ (1/s)' % (a))
    
    axarr[2].set_xlabel('Residue')

    fig.savefig(output)

def plot_pcterr(df, pH_vals, avg, std, output='hdx_k_err_plot.pdf'):
    fig, axarr = plt.subplots(len(pH_vals), sharex=True, figsize=(8.5,11))
    ymin = 0
    ymax = 200
    xmin = 10
    xmax = 143
    for i,n in enumerate(axarr):
        a = pH_vals[i]
        n.bar(df.ix[a]['pcterr'].index, df.ix[a]['pcterr'],
              align='center', color='gray')
        n.set_ylim(ymin, ymax)
        n.set_xlim(xmin, xmax)
        
        n.xaxis.set_major_locator(MaxNLocator(15))
        
        n.axhline(y=100, linewidth=1, color='black')
                                
        n.set_ylabel('%s pct err on $k_{obs}$' % (a))
    
    axarr[2].set_xlabel('Residue')

    fig.savefig(output)
    
def plot_slopes(slope, err, output='hdx_dg_slopes_plot.pdf'):
    fig, ddg = plt.subplots(figsize=(10, 5))
    ddg.grid(True)
    ddg.bar(slope.index, slope, yerr=err, align='center', color='gray')
    dphsddg = 2
    #xmin = 6; xmax = 143; ymin = -0.5; ymax = 8
    #axlim = [xmin, xmax, ymin, ymax]

    ddg.axhline(y=dphsddg, linewidth=3, color='red')
    
    iongrps = ['E10','D19','D21','E57','E67','D73','E75','D83','D95','E101',\
    'H121','E122','E129','E135']

    for i in iongrps:
        loc = int(i.strip('DEH'))
        ddg.annotate(i,xy=(loc,3.5),xytext=(loc,4.5), fontsize=8, rotation=90,\
            arrowprops=dict(facecolor='black',width=0.5,headwidth=2,\
            shrink=0.05))

    ddg.set_ylabel('$\Delta\Delta G_{local}$/$\Delta$pH (kcal/mol/pH)')
    ddg.set_xlabel('Residue')

    fig.savefig(output)    

def hdx_simulation(f, t0='manual'):
    """Simulation of H-D exchange decay to assess accuracy of fitting.
    Generates data from exponential, adds random noise, fits noisy series.
    Prints results as a table to inspect visually.
    
    Parameters
    ==========
    f : string filename (Sparky rh table), provides realistic time values
    t0 : string, default 'manual' applies if timepoints in seconds are already
         encoded in the Sparky rh file columns; otherwise should be a timestamp
         in the format '%Y-%m-%d_%H:%M'
    """
    t_series = nmrfn.parse_hx(f, t0)

    k_input_names = ['5e-05']*4 + ['1e-05']*4 + ['5e-06']*4 + ['1e-06']*4 + \
        ['5e-07']*4 + ['1e-07']*4 + ['5e-08']*4
    noise_scale_names = ['0.01', '0.05', '0.1', '0.2']*7
    arrays = [k_input_names, noise_scale_names]
    tuples = zip(*arrays)
    index = pd.MultiIndex.from_tuples(tuples, names=['k_in', 'noiz'])
    columns=['rat_avg','k_avg','ker_avg','k_std','er_rat','dg_sd']
    df = pd.DataFrame(columns=columns, index=index)
    
    k_in = {'5e-05':5e-5, '1e-05':1e-5, '5e-06':5e-6, '1e-06':1e-6, 
        '5e-07':5e-7, '1e-07':1e-7, '5e-08':5e-8}
    noise_in = {'0.01':0.01, '0.05':0.05, '0.1':0.1, '0.2':0.2}
    a = 7e7
    for kkey, kval in k_in.items():
        cln = nmrfn.decay(t_series, a, kval) # clean data
        for nkey, nval in noise_in.items():
            ratios = []; k_fits = []; k_errs = []
            while len(k_fits) < 1000:
                noisy = cln + np.random.normal(loc=0,scale=nval*a,size=len(cln))
                ratios.append(noisy[-6:-1].mean()/noisy[0:4])
    
                a0 = 4e7; k0 = 1e-6
                popt, pcov = curve_fit(nmrfn.decay, \
                                    xdata=t_series, ydata=noisy, p0=[a0,k0])
    
                k_fits.append(popt[1])
                k_errs.append( np.sqrt(pcov.diagonal()[1]) )
    
            ratios = np.array(ratios)
            k_fits = np.array(k_fits)
            k_errs = np.array(k_errs)
    
            k_avg = k_fits.mean(); k_er = k_errs.mean(); k_sd = k_fits.std()
            RT = 1.986e-3 * 298.0
            df.ix[kkey, nkey]['rat_avg'] = '{0:.2f}'.format(ratios.mean())
            df.ix[kkey, nkey]['k_avg'] = '{0:.2e}'.format(k_avg)
            df.ix[kkey, nkey]['ker_avg'] = '{0:.2e}'.format(k_er)
            df.ix[kkey, nkey]['er_rat'] = '{0:.2f}'.format(k_er/k_avg)
            df.ix[kkey, nkey]['k_std'] = '{0:.2e}'.format(k_sd)
            df.ix[kkey, nkey]['dg_sd'] = '{0:.2f}'.format(RT*k_sd/k_avg)
    
    print(df)
    
    #print 'Fitted k and error: %s +/- %s' % (kfit, kerr)
    
    #plt.ylabel('Peak height (arbitrary units)')
    #plt.xlabel('Time (s)')
    #plt.errorbar(t_series, noisy, fmt='ro') # noisy data
    #plt.plot(t_series,data_helpers.decay(t_series,afit,kfit), 'r-',lw=2)
    #plt.plot(t_series,data_helpers.decay(t_series,a_in,k_in), 'b-')
    #plt.show()
    
if __name__ == '__main__':
    analysis_scheme()