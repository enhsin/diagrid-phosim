import os, re, sys
import subprocess


def showThumbs(start, end, inFile):

    """
    
    Script takes a starting and ending file number and plots
    thumbnails of each file in the given list of files.

    """


    
    myFiles = open(inFile).readlines()
    
    numFiles = len(myFiles)
    print 'Number of Files: ', numFiles
    print 'Starting File Number: ', start
    print 'Ending File Number: ', end
    count = 1
    cmd = '/share/home/nms/bin/ds9 -cmap SLS -scale mode 95 '
    for files in myFiles:
        files = files.strip()
        if count >= int(start) and count <= int(end):
            #cmd = cmd+'-tile '+files+ ' -zoom to 2 -pan to %f %f wcs fk5 -wcs align yes ' %(objRaDeg, objDecDeg)
            cmd = cmd+' -tile '+files
        count += 1
    #print cmd    
    subprocess.check_call(cmd, shell=True)

if __name__== "__main__":
    
    if not len(sys.argv) == 4:
        print "usage: python ds9Thumb.py start end inFile"
        quit()
        
    start = sys.argv[1]
    end = sys.argv[2]
    inFile = sys.argv[3]
    
    showThumbs(start, end, inFile)   
