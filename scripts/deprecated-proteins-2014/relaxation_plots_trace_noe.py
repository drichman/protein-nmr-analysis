#!/usr/bin/python

from os import listdir, getcwd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

c = {'ph3p3':'r', 'ph4p7':'g', 'ph7p4':'b'}

raw = plt.figure(figsize=(8.5,7.5)) # figsize args are w,h in inches
G = gridspec.GridSpec(2, 3)
noe1 = plt.subplot(G[0,0])
noe2 = plt.subplot(G[0,1], sharey=noe1)
noe3 = plt.subplot(G[0,2], sharey=noe1)
noe4 = plt.subplot(G[1,0], sharey=noe1)
noe5 = plt.subplot(G[1,1], sharey=noe1)
noe6 = plt.subplot(G[1,2], sharey=noe1)

rawf = [i for i in listdir(getcwd()) if i.endswith('rlx_raw')]
for f in rawf:
    m = f.split('_')[0]
    n = m # n can be assigned as DataFrame while m is preserved is string
    n = pd.read_table(f, delim_whitespace=True, usecols = [0,5,6,7,8,9,10], \
        names = ['res','r1','r1r','r2','r2r','noe','noer'], \
        header = None, index_col = 'res')
    data = n.reindex(index=range(1,150))
    noe1data = data.iloc[10:21,:]
    noe2data = data.iloc[30:41,:]
    noe3data = data.iloc[54:62,:]
    noe4data = data.iloc[61:71,:]
    noe5data = data.iloc[80:96,:]
    noe6data = data.iloc[102:112,:]
    psz = 5 # Scatter plot marker size
    noe1.plot(noe1data.index, noe1data.noe, 'o-', markersize=psz, color=c[m])
    noe2.plot(noe2data.index, noe2data.noe, 'o-', markersize=psz, color=c[m])
    noe3.plot(noe3data.index, noe3data.noe, 'o-', markersize=psz, color=c[m])
    noe4.plot(noe4data.index, noe4data.noe, 'o-', markersize=psz, color=c[m])
    noe5.plot(noe5data.index, noe5data.noe, 'o-', markersize=psz, color=c[m])
    noe6.plot(noe6data.index, noe6data.noe, 'o-', markersize=psz, color=c[m])
    
seqmin = 6
seqmax = 143
noe1.set_xlim(12, 21)
noe2.set_xlim(32, 41)
noe3.set_xlim(57, 62)
noe4.set_xlim(63, 71)
noe5.set_xlim(82, 96)
noe6.set_xlim(104, 112)
noe1.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe2.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe3.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe4.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe5.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe6.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
plt.setp(noe2.get_yticklabels(), visible = False)
plt.setp(noe3.get_yticklabels(), visible = False)
plt.setp(noe5.get_yticklabels(), visible = False)
plt.setp(noe6.get_yticklabels(), visible = False)
noe5.set_xlabel('Residue')

noemin = 0.63
noemax = 0.88
noe1.set_ylim(noemin, noemax)
noe1.yaxis.set_major_locator(MaxNLocator(8, prune='both'))
noe1.set_ylabel('NOE')
noe4.set_ylabel('NOE')

plt.tight_layout(h_pad=4.0)

plotfile = PdfPages('plots_noe_main.pdf')
plotfile.savefig(raw)
plotfile.close()