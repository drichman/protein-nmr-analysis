#!/bin/bash
# DER 2012
# For a processing a batch of NMR directories (eg relaxation series)
# using the same .com and ft2.com. First, make and run fid.com and
# ft2.com in one of the experiment directories. Then go to the parent
# directory and call this script.

find */fid.com -exec cp {} . \;
find */ft2.com -exec cp {} . \;
ls -d * | sed 's/\(.*\)/cp fid.com \1\/fid.com/' | csh
ls -d * | sed 's/\(.*\)/cp ft2.com \1\/ft2.com/' | csh
rm fid.com
rm ft2.com

for i in $( ls -d * ); do
    cd $i
    fid.com
    ft2.com
    pipe2ucsf test.ft2 ../$i.ucsf
    cd ../
done

# An alternative to this loop is
# find -name fid.com (or ft2.com) -execdir {} \;
# but this seems to require removing . from $PATH, using   
# export PATH=`echo $PATH | sed 's/.://'`  
# but I am reluctant to do this all the time, so for now am sticking with the
# inelegant loop.
