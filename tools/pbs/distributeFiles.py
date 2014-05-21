""" distributeFiles
    Given a filename and a top level directory move a file to the
    correct file structure and create directories if necessary

    Author: ajc@astro.washington.edu - 4/6/2010

    Modified to work with python imsim scripts - 6/3/10 NMS  
    """

import os
import sys
import shutil
import glob
import subprocess
import time, datetime
from optparse import OptionParser

def parseImsimFileName(fileName):
    """ Return the path to the file
    Expected filename format is imsim_85770141_f0_R01_S22_C17_E001.fits.gz
    """
    
    #remove directory structure and .fits 
    stubs = os.path.basename(fileName).partition(".")

    #map filter number to filter character
    filtmap = {"f0":"fu", "f1":"fg", "f2":"fr", "f3":"fi", "f4":"fz", "f5":"fy"}
    
    try:
        imageType, obshistid, filter, raft, sensor, channel, exposure = stubs[0].split("_")
        filter = filtmap[filter]
        path = "imSim/PT1.2/raw/v%08d-%s/%s/%s/%s/" % (int(obshistid),filter,\
                                                 exposure, raft, sensor)
        filename = "imsim_%s_%s_%s_%s_%s.fits.gz"%(obshistid,raft,sensor,channel,exposure)
       
    except:
        raise ValueError("Unable to unpack filename %s" % fileName)
    
    #path = "imSim/raw/v%08d/%s/%s/%s/%s/%s" % (int(obshistid),filter,exposure, raft, sensor, channel)
    

    return path, filename

def parseEimageFileName(fileName):
    """ Return the path to the file
    Expected filename format is eimage_85770141_f0_R01_S22_E001.fits.gz
    """
        
    #remove directory structure and .fits 
    stubs = os.path.basename(fileName).partition(".")
    
    #map filter number to filter character
    filtmap = {"f0":"fu", "f1":"fg", "f2":"fr", "f3":"fi", "f4":"fz", "f5":"fy"}
    
    try:
        imageType, obshistid, filter, raft, sensor, exposure = stubs[0].split("_")
        filter = filtmap[filter]
        path = "imSim/PT1.2/eimage/v%08d-%s/%s/%s/" % (int(obshistid),filter,\
                                                 exposure, raft)
        filename = "eimage_%s_%s_%s_%s.fits.gz"%(obshistid,raft,sensor,exposure)
    except:
        raise ValueError("Unable to unpack filename %s" % fileName)

    #path = "imSim/raw/v%08d/%s/%s/%s/%s" % (int(obshistid), filter, exposure, raft, sensor)

    return path, filename

def distributeImsims(topLevelDir, fileName):
    """Move an image file to the MSS directory structure"""

    path, filename = parseImsimFileName(fileName)
    
    # generate directories. Race condition can exist if the directory
    # is create between checking existance and makedirs
    dir =  os.path.join(topLevelDir,path)
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except OSError:
            pass

    # map filter number to filter character
    filtmap = {"f0":"fu", "f1":"fg", "f2":"fr", "f3":"fi", "f4":"fz", "f5":"fy"}

    # move file
    try:
        shutil.move(fileName, os.path.join(dir, filename))
    except:
        raise OSError("Unable to move filename %s" % fileName)

def distributeEimages(topLevelDir, fileName):
    """Move an image file to the MSS directory structure"""

    path, filename = parseEimageFileName(fileName)
    
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
        shutil.move(fileName, os.path.join(dir, filename))
    except:
        raise OSError("Unable to move filename %s" % fileName)

def main():

    usage = "usage: %prog transferPath/base/directory image/path/imageName baseName (imsim or eimage)"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    topLevelDir = args[0]
    fileName = args[1]
    baseName = args[2]

    if baseName == 'imsim':
        distributeImsims(topLevelDir, fileName)
    else:
        distributeEimages(topLevelDir, fileName)

    
if __name__ == '__main__':
    main()
