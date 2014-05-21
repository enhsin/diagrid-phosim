import os
import numpy
import math

filePath = '/share/pogo3/krughoff/TuesdayTrims' 
file = 'trim2544_85501858_short_meta.trim'
filename = os.path.join(filePath, file)

myvals = numpy.loadtxt(filename, delimiter=' ', skiprows=23, usecols=(2,3), unpack=False)

ra = myvals[:,0]    	      
decl = myvals[:,1]

savePath = '/share/home/nms/pt1.2imsimTrunk'
newname = 'trim2544_85501858.reg'
filePathName = os.path.join(savePath, newname)
print 'Writing file: %s' %(newname)
f = open(filePathName, 'w')

for i in range(0,len(ra)-1):
    f.write('fk5; point %sd %sd \n' % (ra[i], decl[i]))  
f.close()
