#!/usr/bin/env python

import sys
import numpy

import lsst.testing.pipeQA         as pipeQA
import lsst.testing.pipeQA.figures as qaFig
from lsst.sims.catalogs.generation.db import jobDB

#############################################################
#
# Main body of code
#
#############################################################

# This tutorial assumes you've read through dispQaTutorial1.py and dispQaTutorial2.py !
def main():
    #There are currently only two states that I know of
    stateEnum = {"RUNNING":0, "FINISHED":1}

    #below is basically a total hack to try to get something to work.  When
    #the state keys are a little more discriptive and state logging is done
    #on a per visit basis, I think this will require less nuttyness.

    jobid = jobDB.JobId(id=77, owner='anon')
    js = jobDB.JobState(jobid=jobid)
    states = js.showStates()
    #Example state key:
    #77_8883317101308JS
    test_obsid = '88545064'
    chips = {}
    rafts = ["0,1", "0,2", "0,3", "1,0", "1,1", "1,2", "1,3", "1,4",\
             "2,0", "2,1", "2,2", "2,3", "2,4", "3,0", "3,1", "3,2",\
             "3,3", "3,4", "4,1", "4,2", "4,3"]
    sensors = ["0,0", "0,1", "0,2", "1,0", "1,1", "1,2", "2,0", "2,1", "2,2"]
    senid = 0
    raftid = 0 
    for k in states.keys():
        if not (k[-2:] == "JS"):
            continue
        obsid = k[3:11]
        if not obsid == test_obsid:
            continue
        chipid = int(k[13:-2])
        snap = int(k[12])
        chips["R:"+rafts[raftid]+" "+"S:"+sensors[senid]+":0"] = {
                "raftname":rafts[raftid], "sensorname":sensors[senid],
                "status":stateEnum[states[k]]}
        senid += 1
        senid %= 9
        if senid == 0:
            raftid += 1
        raftid %= 21

    ##############################
    # Adding a mapped FPA figure
    # - The mapped fpa will appear on the right 'nav' panel
    # - when you click on a ccd, you'll select data to display for that ccd
    # let's add this to a TestSet
    tsMapFpa = pipeQA.TestSet(group="imsimstatus",
            label="mapped-FPA-figure-imsim")
    tsMapFpa.addTest("test3", 1, [2, None], "Top level focalplane level test")
    tsMapFpa.addTest("test4", 1, [0, 2], "Another top level test")
    
    camInfo = pipeQA.LsstSimCameraInfo()
    fpaFig = qaFig.FpaQaFigure(camInfo)

    # Fill the fpaFig.data attribute with the values to display.
    # This will color-code ccds in the focal plane according to these values
    data = fpaFig.data
    map  = fpaFig.map
    print chips.keys()
    for raft in sorted(data.keys()):
        ccdDict = data[raft]
        for ccd, value in ccdDict.items():
            key = ccd+":0"
            print key
            #Above the key was constructed to have the form "R:rx,ry S:sx,sy:snap" 
            #which is how th fpaFig.map keys are constructed
            if not chips.has_key(key):
                continue
            #Get enumerated status
            data[raft][ccd] = chips[key]["status"]
            # this time, add a mouse-over string (no spaces) for the map
            map[raft][ccd]  = "data-from-" + ccd

            # make a test plot for this ccd.  Right now it is just a sine wave
            # with noise, but in the future we could add a grey scale image of
            # the sensor, or similar.
            qafig = qaFig.QaFigure()
            fig = qafig.fig
            ax = fig.add_subplot(111) #1 column, first row, first column
            x = 2.0*numpy.pi*numpy.arange(100)/100
            ax.plot(x, numpy.sin(x)+numpy.random.normal(0.0, 0.1, len(x)), 'r-')
            ax.set_ylim([-1.0, 1.0])
            
            areaLabel = fpaFig.getAreaLabel(raft, ccd) # areaLabels are handled by FpaQaFigure
            tsMapFpa.addFigure(qafig, "sine-wave.png", "Sine wave from "+ ccd, areaLabel=areaLabel)
            # We can also add tests at the ccd level.  I can imagine having a
            # test for file size here.
            tsMapFpa.addTest("test4", 1, [0, 2], "Here will be a test for file size", areaLabel=areaLabel)
            
    # now tell fpaFig to make the matplotlib figure, and add the fpaFigure to the TestSet
    fpaFig.makeFigure(vlimits=[-1.0, 1.0], cmap="gray", cmapUnder="#0000ff", cmapOver="#ff0000")
    tsMapFpa.addFigure(fpaFig, camInfo.name+".png", "FPA for camera: "+camInfo.name, navMap=True)

if __name__ == '__main__':
    main()
    
