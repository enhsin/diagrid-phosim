import os, sys
import subprocess


def unzip(filelist, fileDir):

    filelist = open(filelist).readlines()
    x = 0
    for files in filelist:
        files = files.strip()
        file = os.path.join(fileDir,files) 
        cmd = 'tar xzfPC %s %s' %(file, fileDir) 
        subprocess.check_call(cmd, shell=True)
        x += 1
        print 'Unzipped number %i, file %s.' %(x, file)
        os.remove(file)
        print 'Removed %s', file

    return

if __name__ == "__main__":

    if not len(sys.argv) == 3:
        print "usage: python unzip.py filelist filedir"
        quit()

    filelist = sys.argv[1]
    fileDir = sys.argv[2]

    unzip(filelist, fileDir)
