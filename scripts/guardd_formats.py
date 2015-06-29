# Script to process and convert NMR relaxation dispersion data to GUARDD format
# Dan Richman, 2014
# (GUARDD copyright Ian Kleckner, http://code.google.com/p/guardd/)

import argparse
import nmrfn
from rd_analysis import calc_R2series

# Handling command line arguments
cmdln = argparse.ArgumentParser(description='Prepare input for GUARDD')
cmdln.add_argument('outtype', help='r2 or height')
cmdln.add_argument('rh', help='Sparky rh file name to read in')
cmdln.add_argument('seqlist', help='Protein sequence file')
cmdln.add_argument('fout', help='Filename for resulting GUARDD script')
cmdln.add_argument('CT', type=float, help='Const-time relaxation period in ms')
cmdln.add_argument('name', help='Name that dataset should have in GUARDD')
cmdln.add_argument('ax', help='Relaxation nucleus', choices=['15N','13C'])
cmdln.add_argument('b0', type=float, \
    help='Field (MHz). Brew2:600.53, Brew:600.13, 800:799.7, 500:499.9')
cmdln.add_argument('temp', type=float, help='Temperature in C')
cmdln.add_argument('coh', help='Coherence during relaxation: SQ:TRUE, MQ:FALSE')
cmdln.add_argument('atom', help='Relaxation atom', \
    choices=['NH','C','CO','\delta_1'])
args = cmdln.parse_args()

# Function to print dataset information
def guardd_dataset(name, ax, b0, temp, CT, coh, f):
    text = \
'''DATASET
NAME\t{0}
AX\t{1}
B0\t{2}
TEMPC\t{3}
TCPMG\t{4}
SQX\t{5}
SETSPECS
'''
    print(text.format(name, ax, b0, temp, round(CT*1e-3, 4), coh), file = f)

# Function to print curve information and R2 and error vs freq
def guardd_r2(s, f):
    text = \
'''INDEX\t{0}
ATOM\t{1}
RESIDUE\tpending
OBS\tVCPMG\tR2\tERROR'''
    print(text.format(s.name, args.atom), file = f)
    for o, t in enumerate(s.index.values, start=1):
        line = '{0}\t{1}\t{2}'
        print(line.format(o, round(t,2), round(s.ix[t],2)), file = f)
    print('ADDDATA\n', file = f)

# Function to print curve information and normalized peak intensity vs freq
def guardd_height(s, f):
    text = \
'''INDEX\t{0}
ATOM\t{1}
RESIDUE\tpending
OBS\tVCPMG\tINTENSITY'''
    print(text.format(s.name, args.atom), file = f)
    norm = s.values[0]
    for o, t in enumerate(s.index.values, start=1):
        line = '{0}\t{1}\t{2}'
        print(line.format(o, round(t,2), round(s.iloc[o-1]/norm,3)), file = f)
    print('ADDDATA\n', file = f)

# Actual output procedure
out = open(args.fout, 'w')
guardd_dataset(args.name,args.ax,args.b0,args.temp,args.CT,args.coh,out)
heights = nmrfn.parse_rd(args.rh)

if args.outtype == 'r2':
    r2_vals = calc_R2series(heights, args.CT)
    for res in r2_vals:
        guardd_r2(r2_vals[res], out)
if args.outtype == 'height':
    for res in heights:
        guardd_height(heights[res], out)

print('SEQUENCEFILE\t{0}'.format(args.seqlist), file = out)
out.close()