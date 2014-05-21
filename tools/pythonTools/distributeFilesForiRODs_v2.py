""" distributeFilesForiRODS
    Given a filename and a top level directory move a file to the
    correct file structure and create directories if necessary

    ajc@astro.washington.edu 4/6/2010

    5/24/10: added iRODS registration - NMS  
    """

import os
import sys
import shutil
import glob
import subprocess
import time, datetime
from optparse import OptionParser

def parse_options(argv):
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage=usage)
    #parser.add_option("-t", "--topLevelDir", action="store", type="string", dest="topLevelDir",  help="Directory on local storage for iRODS.", default="/share/pogo1/iRODS/PT1/")
    #parser.add_option("-f", "--filename", action="store", type="string", dest="filename",  help="Director to imsim_*.gz files.", default=None)
    parser.add_option("-g", "--gIrodsDir", action="store", type="string", dest="gIrodsDir",  help="Directory on global iRODS storage.", default="/ImSimGrid/ImSimData/global/")
    parser.add_option("-l", "--lIrodsDir", action="store", type="string", dest="lIrodsDir",  help="Directory on local iRODS storage.", default="/share/pogo1/iRODS/PT1/")
    (options, argv) = parser.parse_args()
    return options


def parseFileName(fileName):
    """ Return the path to the file

    Expected filename format is imsim_85770141_f0_R01_S22_C17_E001.fits.gz
    Expected filename format is eimage_85770141_f0_R01_S22_E001.fits.gz
    """
    
    #remove directory structure and .fits 
    stubs = os.path.basename(fileName).partition(".")
    
    try:
        imageType, obshistid, filter, raft, sensor, channel, exposure = stubs[0].split("_")
    except:
        raise ValueError("Unable to unpack filename %s" % fileName)
    
    path = "imSim/raw/v%08d/%s/%s/%s/%s/%s" % (int(obshistid),filter,exposure, raft, sensor, channel)

##     try:
##         imageType, obshistid, filter, raft, sensor, exposure = stubs[0].split("_")
##     except:
##         raise ValueError("Unable to unpack filename %s" % fileName)

##     path = "imSim/raw/v%08d/%s/%s/%s/%s" % (int(obshistid), filter, exposure, raft, sensor)
    
    return path 

def distributeFilesForiRODs(topLevelDir, fileName, gIrodsDir, lIrodsDir):
    """Move an image file to the iRODS directory structure"""

    path = parseFileName(fileName)
    
    # generate directories. Race condition can exist if the directory
    # is create between checking existance and makedirs
    dir =  os.path.join(topLevelDir,path)
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except OSError:
            pass

    # move file
    try:
        shutil.move(fileName, dir)
    except:
        raise OSError("Unable to move filename %s" % fileName)

##     # register with irods
##     localDir = os.path.join(lIrodsDir, dir)
##     print "Local path to file:", localDir
##     globalDir = os.path.join(gIrodsDir, dir)
##     print "Global path to file:", globalDir
##     try:
##         cmd = "iput -f -R ncsa-disk1 %s/%s %s/%s" %(localDir, filename, globalDir, filename)
##         # subprocess.checkZ_call(cmd, shell=True)
        
##         # write a logfile to record transferred files
##         tempDate = datetime.date.today()
##         sDate = str(tempDate)
##         year, mo, day = sDate.split('-')
##         logFile = '%s_%s_%s_%s%s%s_transferred.iRodsLog' %(obshistid, filter, exposure, year, mo, day) 
##         with file(logFile, 'w') as myFile:
##             myFile.write('%s/%s' %(localDir, filename))
##     except:
##         # write a logfile to record failed file transfers
##         logFile = '%s_%s_%s_%s%s%s_failed.iRodsLog' %(obshistid, filter, exposure, year, mo, day) 
##         with file(logFile, 'w') as myFile:
##             myFile.write('%s' %(cmd))
    
def main():

    usage = "usage: %prog toplevelDir filename (or wildcard) -iRODS"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    try:
        topLevelDir = args[0]
        # assume that system will expand wild cards into the args list
        count = 0
        for fileName in args[1:]:
            print "Filename: ", fileName
            distributeFilesForiRODs(topLevelDir, fileName, gIrodsDir, lIrodsDir)
            count += 1
            print count
    except ValueError, (ErrorMessage):
        print "Parse Error %s" % ErrorMessage
    except OSError, (ErrorMessage):
        print "IO Error %s", ErrorMessage
    except IndexError, (ErrorMessage):
        print "IndexError: USAGE distributeFilesForiRODS toplevelDir filename (or wildcard)"
        
if __name__ == '__main__':
    options  = parse_options(sys.argv)
    #topLevelDir = options.topLevelDir
    #filename = options.filename
    gIrodsDir = options.gIrodsDir
    lIrodsDir = options.lIrodsDir
    main()
