#!/usr/bin/env python

import sys
## import matplotlib
## matplotlib.use('Agg')
## import aplpy
import numpy
import pyfits
import pylab as py
#import img_scale

import lsst.testing.pipeQA         as pipeQA
import lsst.testing.pipeQA.figures as qaFig
from lsst.sims.catalogs.generation.db import jobDB

def main():
    
    stateEnum = {"RUNNING":0, "FINISHED":1, "FAILED":2}
    jobid = jobDB.JobId(id=9999, owner='nms')
    js = jobDB.JobState(jobid=jobid)
    states = js.showStates()
    
    #Example state key:
    #885448631_R:3,0 S:1,1:0_JS
  
    chips = {}

    for k in states.keys():
        if not (k[-2:] == "JS"):
            continue
        obsid = k[0:9]
        chipname = k[10:23]
        raftid = k[12:15]
        sensid = k[18:21]
        chips[chipname] = {"raftname":raftid, "sensorname":sensid,
                           "status":stateEnum[states[k]]}
        
    ##############################
    # Adding a mapped FPA figure
    # - The mapped fpa will appear on the right 'nav' panel
    # - when you click on a ccd, you'll select data to display for that ccd

    tsMapFpa = pipeQA.TestSet(group="%s" %(obsid),
            label="FPA-imsim")
    #tsMapFpa.addTest("test3", 1, [2, None], "Top level focalplane level test")
    #tsMapFpa.addTest("test4", 1, [0, 2], "Another top level test")
    tsMapFpa.addTest("%s" %(obsid), 1, [0, 1], "Top Level Test for visit %s." %(obsid))
    
    camInfo = pipeQA.LsstSimCameraInfo()
    fpaFig = qaFig.FpaQaFigure(camInfo)

    # Fill the fpaFig.data attribute with the values to display.
    # This will color-code ccds in the focal plane according to these values
    data = fpaFig.data
    map  = fpaFig.map
    #print chips.keys()
    for raft in sorted(data.keys()):
        ccdDict = data[raft]
        for ccd, value in ccdDict.items():
            key = ccd+":0"
            #print key
            #Above the key was constructed to have the form "R:rx,ry S:sx,sy:snap" 
            #which is how the fpaFig.map keys are constructed
            if not chips.has_key(key):
                continue
            #Get enumerated status
            data[raft][ccd] = chips[key]["status"]
            # this time, add a mouse-over string (no spaces) for the map
            map[raft][ccd]  = "data-from-" + ccd

            # make a test plot for this ccd.  Right now it is just a sine wave
            # with noise, but in the future we could add a grey scale image of
            # the sensor, or similar.
##             qafig = qaFig.QaFigure()
##             fig = qafig.fig


##             myimg = pyfits.getdata('imsim_885335891_R01_S00_C00_E000.fits')
##             py.imshow(myimg)

            
##             areaLabel = fpaFig.getAreaLabel(raft, ccd) # areaLabels are handled by FpaQaFigure
## ##             #tsMapFpa.addFigure(qafig, "sine-wave.png", "Sine wave from "+ ccd, areaLabel=areaLabel)

            
## ##             tsMapFpa.addFigure(qafig, 'myfig.png', "Grayscale image for"+ ccd, areaLabel=areaLabel)
## ##             # We can also add tests at the ccd level.  I can imagine having a
## ##             # test for file size here.

##             testName = 'EimageSizeTest'
##             value = 5.6 #Mb
##             range = [4.4, 6.4] #Mb
  
##             test1 = 

##             tsMapFpa.addTest(testName, value, range, "Test the file size of the eimage.", areaLabel=areaLabel)
            
##             tsMapFpa.addTest("test", 1, [0, 2], "Here will be a test for file size", areaLabel=areaLabel)
            
    # now tell fpaFig to make the focalplane matplotlib figure, and add the fpaFigure to the TestSet
    fpaFig.makeFigure(vlimits=[-1.0, 1.0], cmap="gray", cmapUnder="#0000ff", cmapOver="#ff0000")
    tsMapFpa.addFigure(fpaFig, camInfo.name+".png", "FPA for camera: "+camInfo.name, navMap=True)
    

if __name__ == '__main__':
    main()
    
