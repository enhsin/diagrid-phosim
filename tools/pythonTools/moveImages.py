#!/opt/rocks/bin/python

from __future__ import with_statement
import os, sys, time
import shutil
import subprocess
from optparse import OptionParser

def parse_options(argv):
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-o", "--obsid", action="store", type="string", dest="obsid",  help="Givethe obsId of the CCDs to check", default=None)
    (options, argv) = parser.parse_args()
    return options

def checkCcds(obsId):
  
    rxlist = ['0','1','2','3','4']
    rylist = ['0','1','2','3','4']
    sxlist = ['0','1','2']
    sylist = ['0','1','2']
    pairlist = ['0','1']
    ampxlist = ['0','1']
    ampylist = ['0','1','2','3','4','5','6','7']
    errors = []
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
                            id = 'R'+rx+ry + '_' + 'S'+sx+sy + '_' +'E00'+ex
                            ccd = obsid + '_' + id
                            eimage = ccd + '/eimage_' + ccd + '.fits.gz'
                            efitsFile = 'eimage_' + ccd + '.fits.gz'
                            
                            # Verify that the CCD directory exists.
                            if os.path.isdir(ccd):
                                print 'CCD: ', ccd
                              
                                # Verify that the eimage exists
                                # If it exists, move it to the base directory
                            
                                if os.path.isfile(eimage):
                                    print 'eimage :', eimage
                                    print 'moving %s to %s' %(eimage, ccd)
                                    shutil.move('%s' %(eimage), '.')
                                    
                                # Verify that the individual amp images exists.
                                # If the amp exists, move it to the base directory.
                                for ax in ampxlist:
                                    for ay in ampylist: 
                                        ampid = 'R'+rx+ry+'_'+'S'+sx+sy+'_'+'C'+ax+ay+'_'+'E00'+ex+'.fits.gz'
                                        amp = ccd + '/imsim_' + obsid + '_' + ampid
                                        if os.path.isfile(amp):
                                            print 'AMP :', amp
                                            print 'moving %s to %s' %(amp, ccd)
                                            shutil.move('%s' %(amp),'.')
                                            
                                # Remove all of the now empty directories
                                cmd = os.rmdir('%s' %(ccd))
                                
                                  
if __name__ == "__main__":
    options  = parse_options(sys.argv)
    obsid = options.obsid

    checkCcds(obsid)
