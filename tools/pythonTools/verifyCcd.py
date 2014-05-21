from __future__ import with_statement
import os, sys

def checkCcds(myfile, path2imsim):
    """
    
    Script takes a list of visit directories and a path to create a
    list of missing directories/files from the distributed imsim
    directory structure and also creats a resubmit list with paths to
    the location of the pbs files.

    Script only checks for the existance of the file - does not
    validtae the file's integrity.
    
    myfile: file with list of visit dirs as
    v[obshistid]-f[filter-letter] (eg.v886469580-fz)
    
    path2imsim: the absolute path to the imSim directory
    location.

    """
    visitDirs = open(myfile).readlines()
    
    for visits in visitDirs:
        visits = visits.strip()
        vobshist, ffilt = visits.split('-')
        obshistid = vobshist.replace('v', '')
        filter = ffilt.replace('f', '')
    
        visit = '%s-f%s' %(obshistid, filter)
        filename = '%s_missing.txt' %(visit)
        resubmitList = 'resubmit_%s_missing.txt' %(visit)
    
        rxlist = ['0','1','2','3','4']
        rylist = ['0','1','2','3','4']
        sxlist = ['0','1','2']
        sylist = ['0','1','2']
        pairlist = ['0','1']
        ampxlist = ['0','1']
        ampylist = ['0','1','2','3','4','5','6','7']
        rawDir = '%s/imSim/PT1.2/raw/v%s' %(path2imsim, visit)
        eimageDir = '%s/imSim/PT1.2/eimage/v%s' %(path2imsim, visit)
        for ex in pairlist:
            for rx in rxlist:
                for ry in rylist:
                    for sx in sxlist:
                        for sy in sylist:
                            if rx+ry == '00':
                                print 'Skipping Invalid Raft Id:', rx+ry
                                break
                            elif rx+ry == '04':
                                print 'Skipping Invalid Raft Id:', rx+ry
                                break
                            elif rx+ry == '40':
                                print 'Skipping Invalid Raft Id:', rx+ry
                                break
                            elif rx+ry == '44':
                                print 'Skipping Invalid Raft Id:', rx+ry
                                break
                            else:                      
                                snap = 'E00'+ex
                                raft = 'R'+rx+ry
                                ccd = 'S'+sx+sy
                                ccdRawDir = os.path.join(rawDir, snap, raft, ccd)
                                ccdEimageDir = os.path.join(eimageDir, snap, raft)
                                eimage = 'eimage_%s_%s_%s_%s.fits.gz' %(obshistid, raft, ccd, snap)
                                pbsPath = '%s/%s/run%s/pbs_%s_%s_%s_%s.pbs' %(path2imsim, visit, obshistid, obshistid, raft, ccd, snap)
                                # Check the eimages
                                if os.path.isdir(ccdEimageDir):
                                    eimagePath = os.path.join(ccdEimageDir, eimage)
                                    if os.path.isfile(eimagePath):
                                        continue
                                    else:
                                        with file(filename, 'a') as myFile:
                                            myFile.write('%s is missing. \n' %(eimage))

                                      #  with file(resubmitList, 'a') as myFile:
                                      #      myFile.write('%s \n' %(pbsPath))
                                else:
                                    with file(filename, 'a') as myFile:
                                        myFile.write('%s is missing.\n' %(ccdEimageDir))
                                    
                                  #  with file(resubmitList, 'a') as myFile:
                                  #      myFile.write('%s \n' %(pbsPath))
                                # Check the amplifier images
                                if os.path.isdir(ccdRawDir):
                                    for ax in ampxlist:
                                        for ay in ampylist:
                                            amp = 'C'+ax+ay
                                            ampimage = 'imsim_%s_%s_%s_%s_%s.fits.gz' %(obshistid, raft, ccd, amp, snap)
                                            ampPath = os.path.join(ccdRawDir, ampimage)
                                            if os.path.isfile(ampPath):
                                                continue
                                            else:
                                                with file(filename, 'a') as myFile:
                                                    myFile.write('%s is missing. \n' %(ampimage))
                                                    
                                                with file(resubmitList, 'a') as myFile:
                                                    myFile.write('%s \n' %(pbsPath))

                                else:
                                    with file(filename, 'a') as myFile:
                                        myFile.write('%s is missing.\n' %(ccdRawDir))
                            
                                    with file(resubmitList, 'a') as myFile:
                                        myFile.write('%s \n' %(pbsPath))
                                    

if __name__ == "__main__":
    
    if not len(sys.argv) == 3:
        print "usage: python verifyCcd.py myfile path2imsim"
        quit()

    myfile = sys.argv[1]
    path2imsim = sys.argv[2]

    checkCcds(myfile, path2imsim)
