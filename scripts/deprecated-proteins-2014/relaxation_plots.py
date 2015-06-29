#!/usr/bin/python

from os import listdir, getcwd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

c = {'ph3p3':'r', 'ph4p7':'g', 'ph7p4':'b'}
psz = 4 # scatter plot marker size

raw = plt.figure(figsize=(8.5,11)) # figsize args are w,h in inches
r1 = raw.add_subplot(311)
r2 = raw.add_subplot(312, sharex=r1)
noe = raw.add_subplot(313, sharex=r1)

def parse_file(f):
    """Read text file with relax. values; return DataFrame indexed by sequence.
    
    Parameters
    ----------
    f : string filename
    
    Returns
    -------
    DataFrame indexed by sequence, with relaxation values or NaNs if no data for
        that residue
    """
    data = pd.read_table(f, delim_whitespace=True, usecols = [0,5,6,7,8,9,10], \
        names = ['res','r1','r1r','r2','r2r','noe','noer'], \
        header = None, index_col = 'res')
    data = data.reindex(index=range(1,150))

extract_pH = lambda x: x.split('_')[0]

def build_data(raw_file=None):
    """Build a DataFrame containing all data to be used in plotting.
    
    Parameters
    ----------
    raw_file : Specifies file, or default None enables scanning the directory
    
    Returns
    -------
    DataFrame containing relaxation parameters, indexed by protein sequence, 
        multi-indexed by pH
    """
    if raw_file is None:
        raw_file = [i for i in listdir(getcwd()) if i.endswith('rlx_raw')]
    raw = [parse_file(filename) for filename in raw_file]
    pH_val = map(extract_pH, raw_file)
    raw = pd.concat(raw, keys=pH_val)

r1.errorbar(raw.index, raw.r1, yerr=raw.r1r,fmt='o', markersize=psz, color=c[m])
r2.errorbar(raw.index, raw.r2, yerr=raw.r2r,fmt='o', markersize=psz, color=c[m])
noe.errorbar(raw.index, raw.noe,yerr=raw.noer,fmt='o',markersize=psz,color=c[m])

seqmin = 6
seqmax = 143
r1.set_xlim(seqmin, seqmax)
r1.xaxis.set_major_locator(MaxNLocator(15))
plt.setp(r1.get_xticklabels(), visible = False)
plt.setp(r2.get_xticklabels(), visible = False)
noe.set_xlabel('Residue')

r1min = 0.95
r1max = 1.55
r1.set_ylim(r1min, r1max)
r2min = 7
r2max = 17
r2.set_ylim(r2min, r2max)
noemin = 0.6
noemax = 0.95
noe.set_ylim(noemin, noemax)
r1.yaxis.set_major_locator(MaxNLocator(8, prune ='both'))
r2.yaxis.set_major_locator(MaxNLocator(8, prune ='both'))
noe.yaxis.set_major_locator(MaxNLocator(9, prune ='both'))
r1.set_ylabel('$R_1$ (1/s)')
r2.set_ylabel('$R_2$ (1/s)')
noe.set_ylabel('NOE')

r1.grid(True, which = 'major')
r2.grid(True, which = 'major')
noe.grid(True, which = 'major')

plt.tight_layout(pad = 8.0, h_pad = 0.5)

# Parameters for dynamics results (S2, tau_local, and Rex) plots
s2min = 0.5
s2max = 1.0
tau_topmin = 0.2
tau_topmax = 3.0
tau_botmin = 0.0
tau_botmax = 0.055
rexmin = 0.0
rexmax = 6.0

dyn = plt.figure(figsize=(8.5,11))
s2 = dyn.add_subplot(411)
tau_top = dyn.add_subplot(412, sharex=s2)
tau_bot = dyn.add_subplot(413, sharex=s2)
# http://matplotlib.1069221.n5.nabble.com/a-break-in-the-y-axis-td13343.html
rex = dyn.add_subplot(414, sharex=s2)

dynf = [i for i in listdir(getcwd()) if i.startswith('dynamics_')]
for f in dynf:
    m = f.split('_')[1]
    n = m
    n = pd.read_fwf(f, colspecs = [(0, 4), (29, 34), (35, 44), (45, 50)], \
        names=['res','s2','tau','rex'], index_col = 'res')
    s2.plot(n.index, raw.s2, 'o', markersize = psz, color = c[m])
    tau_top.plot(n.index, raw.tau, 'o', markersize = psz, color = c[m])
    tau_bot.plot(n.index, raw.tau, 'o', markersize = psz, color = c[m])
    rex.plot(n.index, raw.rex, 'o', markersize = psz, color = c[m])

s2.set_xlim(seqmin, seqmax)
s2.xaxis.set_major_locator(MaxNLocator(15))
tau_top.xaxis.tick_top()
tau_top.tick_params(labeltop='off')
tau_top.spines['bottom'].set_visible(False)
tau_bot.spines['top'].set_visible(False)
tau_bot.xaxis.tick_bottom() 
plt.setp(s2.get_xticklabels(), visible = False)
plt.setp(tau_bot.get_xticklabels(), visible = False)

s2.set_ylim(s2min, s2max)
tau_top.set_ylim(tau_topmin, tau_topmax)
tau_bot.set_ylim(tau_botmin, tau_botmax)
rex.set_ylim(rexmin, rexmax)
s2.yaxis.set_major_locator(MaxNLocator(6, prune ='both'))
tau_top.yaxis.set_major_locator(MaxNLocator(6, prune ='upper'))
tau_bot.yaxis.set_major_locator(MaxNLocator(6, prune ='lower'))
rex.yaxis.set_major_locator(MaxNLocator(7, prune ='both'))
s2.set_ylabel('$S_2$')
tau_top.set_ylabel('$\\tau_{local}$ (ns)')
tau_bot.set_ylabel('$\\tau_{local}$ (ns)')
rex.set_ylabel('$R_ex$ (1/s)')
rex.set_xlabel('Residue')

d = .015 # Length of diagonal break lines in tau axis
kwargs = dict(transform=tau_top.transAxes, color='k', clip_on=False)
tau_top.plot((-d,+d),(-d,+d), **kwargs)      # Top-left diagonal 
tau_top.plot((1-d,1+d),(-d,+d), **kwargs)    # Top-right diagonal 
kwargs.update(transform=tau_bot.transAxes)   # Switch to the bottom axes 
tau_bot.plot((-d,+d),(1-d,1+d), **kwargs)    # Bottom-left diagonal 
tau_bot.plot((1-d,1+d),(1-d,1+d), **kwargs)  # Bottom-right diagonal 

s2.grid(True, which = 'major')
tau_top.grid(True, which = 'major')
tau_bot.grid(True, which = 'major')
rex.grid(True, which = 'major')

plt.tight_layout(pad = 8.0, h_pad = 0.5)

plotfile = PdfPages('plots.pdf')
plotfile.savefig(raw)
plotfile.savefig(dyn)
plotfile.close()