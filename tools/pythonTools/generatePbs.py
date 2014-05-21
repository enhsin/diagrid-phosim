#!/share/apps/lsst_gcc440/Linux64/external/python/2.5.2/bin/python

"""
Brief:   Python script to generate parameter files and PBS scripts for
         each sensor job.  Script calls fullFocalplanePbs.py for this purpose.
         Can be run standalone or called by generateVisitPbs.py.

Usage:   python generatePbs.py [options]
Options: -f: Name of file containing trim file list
         -p: Name of policy file
         -e: Name of additional parameter file
         
Date:    May 03, 2010
Author:  Nicole Silvestri, U. Washington, nms21@uw.edu
Updated: March 11, 2011

Notes:   1. The trimfiles.lis can contain either the relative or absolute paths to the
            trimfiles. (eg. data/catalogs/trim1234_567890_0.gz)

         2. You should modify imsimPbsPolicy.paf to the appropriate settings for
            your account before running this script.
"""

import os, re, sys
import fullFocalplanePbs
import time, datetime
import shutil
import subprocess
import string
from optparse import OptionParser
import lsst.pex.policy as pexPolicy
#from lsst.sims.catalogs.generation.db import jobDB

def parse_options(argv):
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-f", "--file",    action="store", type="string", dest="file",    help="Name of file containing the trim files", default="trimfiles.lis")
    parser.add_option("-p", "--policy",  action="store", type="string", dest="policy",  help="Imsim PBS Policy File Name", default="imsimPbsPolicy.paf")
    parser.add_option("-e", "--extraid", action="store", type="string", dest="extraidFile",  help="Extra parameter file name.", default="")
    (options, argv) = parser.parse_args()
    return options

def generateParameters(file, policy, extraidFile):

    """

    Run the fullFocalplanePbs.py script, populating it with the
    correct user and cluster job submission information from an LSST
    policy file. 
    
    """

    # Read trimfiles from your file list
    myFiles = '%s' %(file)
    files = open(myFiles).readlines()

    for trimfiles in files:
        trimfiles = trimfiles.strip()
        trimfile = os.path.basename(trimfiles)
        basename, extension = os.path.splitext(trimfiles)
        #if os.path.isfile(trimfiles):
        if extension == '.gz':
            print 'Unzipping trimfile: ', basename
            cmd = 'gunzip %s' %(trimfiles)
            subprocess.check_call(cmd, shell=True)
            print 'Running fullFocalplanePbs.py on: ', basename
            fullFocalplanePbs.main(basename, policy, extraidFile)
##             cmd = 'gzip %s' %(basename)
##             subprocess.check_call(cmd, shell=True)
        else:
            print 'Running fullFocalplanePbs.py on: ', trimfiles
            fullFocalplanePbs.main(trimfiles, policy, extraidFile)
##             cmd = 'gzip %s' %(basename)
##             subprocess.check_call(cmd, shell=True)
    return

if __name__ == "__main__":
    options     = parse_options(sys.argv)
    file        = options.file
    policy      = options.policy
    extraidFile = options.extraidFile
    
    generateParameters(file, policy, extraidFile)
 

