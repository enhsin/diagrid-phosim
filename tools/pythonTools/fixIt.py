#!/opt/rocks/bin/python

from __future__ import with_statement 
import os, sys
import string
import re
import subprocess

myFileName = sys.argv[1]
pbsDir = '/home/nms/pt1.2imsimTrunk/'

with file(myFileName, 'rU') as myFile:
    for line in myFile:
        pbsFileName = line.strip()
        if not pbsFileName or pbsFileName.startswith("#"):
            continue
        pbsFileName = pbsFileName.strip()

        with file(pbsFileName, 'rU') as myPbs:
            print 'Working on file %s.' %(myPbs)
            myStr = myPbs.read()
            newString = string.replace(myStr,"../ancillary/sky/sed_flat.txt","flatSED/sed_flat.txt")
            #newString = string.replace(myStr,"../ancillary/sky/sed_flat.txt","flatSED/sed_flat.txt", 1) #replaces just one occurrance
            #mFile = file(pbsPath, 'w')
            mFile = file(pbsFileName, 'w')
            mFile.write(newString)

            print 'Fixed %s.' %(myPbs)
            


    
    
