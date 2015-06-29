#!/usr/bin/python

from os import listdir, getcwd
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

c = {'ph3p3':'r', 'ph4p7':'g', 'ph7p4':'b'}
psz = 4 # scatter plot marker size

raw = pl.figure(figsize=(8.5,11)) # figsize args are w,h in inches
r1 = raw.add_subplot(311)
r2 = raw.add_subplot(312, sharex=r1)
noe = raw.add_subplot(313, sharex=r1)

rawf = [i for i in listdir(getcwd()) if i.endswith('rlx_raw')]
for f in rawf:
    m = f.split('_')[0]
    n = m # n can be assigned as DataFrame while m is preserved is string
    n = pd.read_table(f, delim_whitespace=True, usecols = [0,5,6,7,8,9,10], \
        names = ['res','r1','r1r','r2','r2r','noe','noer'], \
        header = None, index_col = 'res')
    data = n.reindex(index=range(1,150))
    r1.plot(data.index, data.r1, '-', markersize=psz, color=c[m])
    r2.plot(data.index, data.r2, '-', markersize=psz, color=c[m])
    noe.plot(data.index, data.noe, '-', markersize=psz, color=c[m])
    print m
    r1pe = 100*data.r1r/data.r1
    r2pe = 100*data.r2r/data.r2
    noepe = 100*data.noer/data.noe
    print 'R1 percent error avg = %f' % r1pe.mean()
    print 'R2 percent error avg = %f' % r2pe.mean()
    print 'NOE percent error avg = %f' % noepe.mean()

seqmin = 6
seqmax = 143
r1.set_xlim(seqmin, seqmax)
r1.xaxis.set_major_locator(MaxNLocator(15))
pl.setp(r1.get_xticklabels(), visible = False)
pl.setp(r2.get_xticklabels(), visible = False)
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

r1.grid(False, which = 'major')
r2.grid(False, which = 'major')
noe.grid(False, which = 'major')

pl.tight_layout(pad = 8.0, h_pad = 0.5)

plotfile = PdfPages('plots_trace.pdf')
plotfile.savefig(raw)
plotfile.close()