#! /Users/drichman/miniconda3/envs/py2/bin/python2.7

__title__ = "nmrclean.py"
__author__ = "Carlos Castaneda <carlos.a.castaneda@gmail.com>"
__version__ = "2009-07-07"
__usage__ = "Usage: nmr_clean.py directory"
__description__ = \
"""
nmr_clean.py 

Removes intermediate NMR processing files, including *.fid, *.ft1, *.ft2, *.ft3,
and directories including fid*/, ft*/, lp*/.  This program assumes that
NMRPipe is the software used for data processing.
    
MODIFIED by Dan Richman to delete pdata, audita.txt, cpdprg3, scon2, spnam*, .temp, .ased, and .par

The program will begin deleting files from 'directory' and descend into all subdirectories
and recursively delete these files.
"""

import os, shutil, sys


def clean_nmr_data_directories(input_root):

	deletable_ext = ('.ucsf','.fid','.ft1','.ft2','.ft3','.temp','.ased','.par','audita.txt')
	deletable_start = ('spnam','scon','cpdprg')
	deletable_dirs = ('fid','ft','lp','pdata')
	
	print 'NMRClean: Removing unnecessary NMR files in:'
	for root, dirs, files in os.walk(input_root):
		print root
		for f in files:
			if f.endswith(deletable_ext):
				os.remove(os.path.join(root,f))
       			elif f.startswith(deletable_start):
                		os.remove(os.path.join(root,f))
		for d in dirs:
			if d in deletable_dirs:
				shutil.rmtree(os.path.join(root,d))

if __name__ == '__main__':
	try:
		input_dir = sys.argv[1]
	except IndexError:
		print __usage__
		sys.exit()
	if not os.path.isdir(input_dir):
		print "\"%s\" does not exist!" % input_dir
		sys.exit()
	clean_nmr_data_directories(input_dir)
