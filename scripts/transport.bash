#!/bin/bash
# Copies code, eg an NMR script, into your working directory

dir="/Users/drichman/GDrive/code/"
name=$1

if [ $# -eq 0 ]; then
	ls -l $dir
	echo GOTTA PICK ONE - usage 'transport <name>'
else
	if [ -e $dir$name ]; then
		cp $dir$name ${PWD}
	else
		ls -l $dir
		echo GOTTA PICK ONE THAT EXISTS
	fi
fi
