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
r2 = plt.subplot(G[1,:])

rawf = [i for i in listdir(getcwd()) if i.endswith('rlx_raw')]
for f in rawf:
    m = f.split('_')[0]
    n = m # n can be assigned as DataFrame while m is preserved is string
    n = pd.read_table(f, delim_whitespace=True, usecols = [0,5,6,7,8,9,10], \
        names = ['res','r1','r1r','r2','r2r','noe','noer'], \
        header = None, index_col = 'res')
    data = n.reindex(index=range(1,150))
    noe1data = data.iloc[10:20,:]
    noe2data = data.iloc[30:40,:]
    noe3data = data.iloc[60:70,:]
    psz = 5 # Scatter plot marker size
    r2.plot(data.index, data.r2, '-', markersize=psz, color=c[m])
    noe1.plot(noe1data.index, noe1data.noe, 'o-', markersize=psz, color=c[m])
    noe2.plot(noe2data.index, noe2data.noe, 'o-', markersize=psz, color=c[m])
    noe3.plot(noe3data.index, noe3data.noe, 'o-', markersize=psz, color=c[m])


seqmin = 6
seqmax = 143
noe1.set_xlim(12, 20)
noe2.set_xlim(32, 40)
noe3.set_xlim(62, 70)
r2.set_xlim(seqmin, seqmax)
noe1.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe2.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
noe3.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
r2.xaxis.set_major_locator(MaxNLocator(15))
plt.setp(noe2.get_yticklabels(), visible = False)
plt.setp(noe3.get_yticklabels(), visible = False)
#noe2.set_xlabel('Residue')
r2.set_xlabel('Residue')

r2min = 7
r2max = 17
r2.set_ylim(r2min, r2max)
noemin = 0.64
noemax = 0.85
noe1.set_ylim(noemin, noemax)
r2.yaxis.set_major_locator(MaxNLocator(8, prune='both'))
noe1.yaxis.set_major_locator(MaxNLocator(8, prune='both'))
r2.set_ylabel('$R_2$ (1/s)')
noe1.set_ylabel('NOE')

plt.tight_layout() # Possible argument: pad=

plotfile = PdfPages('plots_trace_main.pdf')
plotfile.savefig(raw)
plotfile.close()