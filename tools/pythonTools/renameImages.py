#!/opt/rocks/bin/python

from __future__ import with_statement
import os, sys, time
import shutil
import subprocess
from optparse import OptionParser

def parse_options(argv):
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-o", "--obsid", action="store", type="string", dest="obsid",  help="Givethe obsId of the CCDs to fix", default=None)
    parser.add_option("-f", "--filter", action="store", type="string", dest="filter",  help="Give the filterId of the CCDs to fix (f#)", default=None)
    (options, argv) = parser.parse_args()
    return options

def checkCcds(obsId, filter):
  
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
                            #newCcd = obsid + '_' + filter + '_' + id
                            #eimage = newCcd + '/eimage_' + ccd + '.fits.gz'
                            #newEimage = newCcd + '/eimage_' + newCcd + '.fits.gz'
                            eimage = ccd + '/eimage_' + filter + '_' + ccd + '.fits.gz'
                            newEimage = ccd + '/eimage_' + obsid + '_' + filter + '_' + id + '.fits.gz' 
                            
                            # Verify that the CCD directory exists
                            # If it exists, rename it
                            if os.path.isdir(ccd):
##                                 print 'CCD: ', ccd
##                                 print 'moving %s to %s' %(ccd, newCcd)
##                                 shutil.move('%s' %(ccd), '%s'%(newCcd))
                                
                                # Verify that the eimage exists
                                # If it exists, rename it
                                if os.path.isfile(eimage):
                                    print 'eimage :', eimage
                                    print 'moving %s to %s' %(eimage, newEimage)
                                    shutil.move('%s' %(eimage), '%s'%(newEimage))
                                    
                                # Verify that the individual amp images exists.
                                # If the amp exists, rename it.
                                for ax in ampxlist:
                                    for ay in ampylist: 
                                        ampid = 'R'+rx+ry+'_'+'S'+sx+sy+'_'+'C'+ax+ay+'_'+'E00'+ex+'.fits.gz'
                                        #amp = newCcd + '/imsim_' + obsid + '_' + ampid
                                        amp = ccd + '/imsim_' + filter +'_' + obsid + '_' + ampid
                                        #newAmp = newCcd + '/imsim_' + obsid + '_' + filter + '_' + ampid
                                        newAmp = ccd + '/imsim_' + obsid + '_' + filter + '_' + ampid
                                        if os.path.isfile(amp):
                                            print 'AMP :', amp
                                            print 'moving %s to %s' %(amp, newAmp)
                                            shutil.move('%s' %(amp),'%s' %(newAmp))
                                    
if __name__ == "__main__":
    options  = parse_options(sys.argv)
    obsid = options.obsid
    filter = options.filter

    checkCcds(obsid, filter)
