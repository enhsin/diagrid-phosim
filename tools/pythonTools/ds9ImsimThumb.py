import os, re, sys
import subprocess


def showThumbs(start, end, inFile):

    """
    
    Script takes a starting and ending file number and plots
    thumbnails of each file in the given list of files.

    """
    
    # 871731611 - imsim_88827141_f4_R22_S11_C12_E000.fits.gz
    # Random star near center of amplifier 12 image coordinates
    # ra = 21:04:12.345, dec = -05:21:01.07
    objRaDeg = 21.0*15.0 + 4.0/60.0*15.0 + 12.345/3600.0*15.0
    print 'Obj RA (degrees):', objRaDeg
    objDecDeg = -5.0 - 21.0/60.0 - 1.07/3600.0
    print 'Obj DEC (degrees):', objDecDeg
    
    imageFiles = open(inFile).readlines()
    
    numFiles = len(imageFiles)
    print 'Number of Files: ', numFiles
    print 'Starting File Number: ', start
    print 'Ending File Number: ', end
    #start = 1
    #end = 10
    count = 1
    #cmd = 'ds9 -cmap SLS -scale mode 95 '
    cmd = 'ds9 -scale histequ '
    for files in imageFiles:
        files = files.strip()
        if count >= int(start) and count <= int(end):
            #cmd = cmd+'-tile '+files+ ' -zoom to 2 -pan to %f %f wcs fk5 -wcs align yes ' %(objRaDeg, objDecDeg)
            #cmd = cmd+'-tile '+files+ ' -zoom to 1 '
            #cmd = cmd +'-tile ' +files + ' -zoom to 8 -pan to 316.05 -5.367 wcs fk5 '
            cmd = cmd +'-tile ' +files + ' -zoom to 8 -pan to 316.042 -5.3619 wcs fk5 '
            count += 1
##         cmd = cmd +files + ' '
##         count += 1
    subprocess.check_call(cmd, shell=True)

if __name__== "__main__":
    
    if not len(sys.argv) == 4:
        print "usage: python ds9Thumb.py start end inFile"
        quit()
        
    start = sys.argv[1]
    end = sys.argv[2]
    inFile = sys.argv[3]
    
    showThumbs(start, end, inFile)   
