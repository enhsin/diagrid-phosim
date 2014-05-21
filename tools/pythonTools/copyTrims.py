import os, sys
import shutil

def copyFiles(trimlist):

    trims = open(trimlist).readlines()
    trimDir = '/share/pogo3/krughoff/allVarPT1.2/done_trans'
    runDir = '/share/pogo3/nms/pt1.2productionTrims/2770/clouds'

    x = 0
    for trim in trims:
        trim = trim.strip()
        trimfile = os.path.join(trimDir, trim)
        shutil.copy(trimfile, runDir)
        x += 1
        print 'Copied number %i, file %s.' %(x, trimfile)        

if __name__ == "__main__":

    if not len(sys.argv) == 2:
        print "usage: python copyTrims.py trimlist"
        quit()

    trimlist = sys.argv[1]

    copyFiles(trimlist)
