#!/usr/bin/env python
##
## @package phosim.py
## @file phosim.py
## @brief python methods called by phosim script
##
## @brief Created by:
## @author John R. Peterson (Purdue)
##
## @brief Modified by:
## @author Emily Grace (Purdue)
## @author Nathan Todd (Purdue)
## @author En-Hsin Peng (Purdue)
## @author Glenn Sembroski (Purdue)
## @author Jim Chiang (SLAC)
## @author Jeff Gardner (Google)
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## The intention is that script is the only file necessary to run phosim for
## any purpose-- ranging from a single star on your laptop to full fields
## on large-scale computing clusters.  There is no physics in this
## script.  Its only purpose is to move files to the correct place.
## This file is called by the phosim script

import os
import subprocess
import sys, glob, optparse, shutil
import multiprocessing

## print the usage
def usage():
     script=os.path.abspath(__file__)
     os.system(script+' x --help')

##jobChip is a function that run an individual chip for a single exposure
def jobChip(observationID, cid, eid, filt, outputDir, binDir, instrDir, instrument='lsst', run_e2adc=True):
    fid = observationID + '_' + cid + '_' + eid
    segfile = instrDir+'/segmentation.txt'
    runProgram("raytrace < raytrace_"+fid+".pars", binDir)
    removeFile('raytrace_'+fid+'.pars')
    if run_e2adc:
        runProgram("e2adc < e2adc_"+fid+".pars", binDir)
        removeFile('e2adc_'+fid+'.pars')
        for line in open(segfile):
            aid=line.split()[0]
            if cid in line and aid != cid:
                rawImage=instrument+'_a_'+observationID+'_f'+filt+'_'+aid+'_'+eid+'.fits.gz'
                shutil.move(rawImage,outputDir+'/'+rawImage)
    eImage=instrument+'_e_'+observationID+'_f'+filt+'_'+cid+'_'+eid+'.fits.gz'
    shutil.move(eImage,outputDir+'/'+eImage)

## runProgram function calls each of the phosim programs using subprocess.call
#  it raises an exception and aborts if the return code is non-zero.
def runProgram(command, binDir=None, argstring=None):
    myCommand = command
    if binDir is not None:
        myCommand = os.path.join(binDir, command)
    if argstring is not None:
        myCommand += argstring
    if subprocess.call(myCommand, shell=True) != 0:
        raise RuntimeError("Error running %s" % myCommand)

## removeFile deletes files (if they do not exist, it will catch the OSError exception and silently proceed.)
def removeFile(filename):
    try:
         os.remove(filename)
    except OSError:
         pass

## assignPath figures out the whole path
def assignPath(path,phosimDir):
    if path == 'none':
        return 'none'
    full_path = phosimDir + '/' + path
    if os.path.exists(full_path):
        return full_path
    elif os.path.exists(path):
        return os.path.abspath(path)
    raise RuntimeError('%s does not exist.' % path)

## checkPaths makes sure the paths make sense
def checkPaths(opt,phosimDir):
    for x in ['outputDir','workDir','binDir','dataDir','sedDir','imageDir','extraCommands']:
        exec('opt.%s=assignPath(opt.%s,phosimDir)' % (x,x))

    #Now find the instrument directory
    instrDir = os.path.join(opt.dataDir, opt.instrument)
    if not os.path.exists(instrDir):
        raise RuntimeError('%s does not exist.' % instrDir)
    opt.instrDir=instrDir

    opt.instrument = instrDir.split("/")[-1]
    if len(opt.instrument)==0:
        opt.instrument = instrDir.split("/")[-2]


## PhosimFocalplane is a class for handling phosim files and directories for one focalplane.
#    In order to perform all of the preprocessing steps, use doPreproc().  The
#    raytrace steps can then be scheduled (exactly how depends on the 'grid' option)
#    using scheduleRaytrace().  Once the computation is finished, intermediate files
#    can be removed using cleanup().  Thus, an example workflow would be as follows:
#    focalplane = PhosimFocalplane(phosimDir, outputDir, workDir, binDir, dataDir, instrDir, opt.grid, grid_opts)
#        focalplane.doPreproc(instanceCatalog, extraCommands, sensor)
#        focalplane.scheduleRaytrace(instrument, run_e2adc, keep_screens)
#        focalplane.cleanup(keep_screens)
class PhosimFocalplane(object):

     ## Constructor for PhosimFocalplane
     #  grid:      'no', 'condor', 'cluster'
     #  grid_opts: A dictionary to supply grid options.  Exactly which options
     #  depends on the value of 'grid':
     #  'no':      'numproc' = Number of threads used to execute raytrace.
     #  'condor':  'universe' = Condor universe ('vanilla', 'standard', etc)
     #  'cluster': 'script_writer' = callback to generate raytrace batch scripts
     #             'submitter' = optional callback to submit the job
     def __init__(self, phosimDir, opt, grid_opts={}):

        self.phosimDir = phosimDir
        self.outputDir = opt.outputDir
        self.workDir = opt.workDir
        self.binDir = opt.binDir
        self.dataDir = opt.dataDir
        self.instrDir = opt.instrDir
        self.sedDir = opt.sedDir
        self.imageDir = opt.imageDir
        self.flatdir = False
        self.tarfile = False
        self.extraCommands = None
        self.instanceCatalog = None
        self.userCatalog = None
        self.chipID = None
        self.runFlag = None
        self.devtype = None
        self.devvalue = None
        self.grid = opt.grid
        self.grid_opts = grid_opts
        self.execEnvironmentInitialized = False
        if self.grid == 'condor':
            assert 'universe' in self.grid_opts
            self.flatdir=True if self.grid_opts['universe'] == 'vanilla' else False
        elif self.grid == 'diagrid':
            self.flatdir=True
            self.tarfile=True

     ## doPreproc is a method to run all of the non-chip steps.
     def doPreproc(self, instanceCatalog, extraCommands, sensor):
        self.loadInstanceCatalog(instanceCatalog, extraCommands)
        os.chdir(self.workDir)
        self.writeInputParamsAndCatalogs()
        self.generateAtmosphere()
        self.generateInstrumentConfig()
        self.trimObjects(sensor)

     ## Parses the instance catalog
     def loadInstanceCatalog(self, instanceCatalog, extraCommands):
        self.instanceCatalog = instanceCatalog
        self.extraCommands = extraCommands
        defaultCatalog=open(os.path.join(self.phosimDir,'default_instcat')).readlines()
        self.userCatalog=open(instanceCatalog).readlines()
        for line in defaultCatalog+self.userCatalog:
             lstr=line.split()
             if "obshistid" in line:
                  self.observationID=lstr[1]
             elif "moonra" in line:
                  self.moonra=lstr[1]
             elif "moondec" in line:
                  self.moondec=lstr[1]
             elif "sunalt" in line:
                  self.solaralt=lstr[1]
             elif "moonalt" in line:
                  self.moonalt=lstr[1]
             elif "dist2moon" in line:
                  self.moondist=lstr[1]
             elif "moonphase" in line:
                  self.phaseang=lstr[1]
             elif "mjd" in line:
                  self.tai=lstr[1]
             elif "seeing" in line:
                  self.constrainseeing=lstr[1]
             elif "rottelpos" in line:
                  self.spiderangle=lstr[1]
             elif "Azimuth" in line or "azimuth" in line:
                  self.azimuth=lstr[1]
             elif "Altitude" in line or "altitude" in line:
                  self.altitude=lstr[1]
             elif "rotskypos" in line:
                  self.rotationangle=lstr[1]
             elif "Unrefracted_RA" in line or "rightascension" in line:
                  self.pointingra=lstr[1]
             elif "Unrefracted_Dec" in line or "declination" in line:
                  self.pointingdec=lstr[1]
             elif "SEED" in line  or "seed" in line:
                  self.obsseed=lstr[1]
             elif "date" in line:
                  self.monthnum=lstr[1].split('/')[1]
             elif "filter" in line:
                  self.filt=lstr[1]
             elif "VISTIME" in line or "vistime" in line:
                  self.vistime=float(lstr[1])
             elif "NSNAP" in line or "nsnap" in line:
                  self.nsnap=int(float(lstr[1]))
             elif "MINSOURCE" in line or "minsource" in line:
                  self.minNumSources=int(float(lstr[1]))
             elif "CAMCONFIG" in line or "camconfig" in line:
                  self.camconfig=int(float(lstr[1]))
             elif "DOMEINT" in line or "domeint" in line:
                  self.domeint=float(lstr[1])
             elif "DOMEWAV" in line or "domewav" in line:
                  self.domewav=float(lstr[1])
             elif "TELCONFIG" in line or "telconfig" in line:
                  self.telconfig=int(float(lstr[1]))
             elif "TEMPERATURE" in line or "temperature" in line:
                  self.temperature=lstr[1]
             elif "TEMPVAR" in line or "tempvar" in line:
                  self.tempvar=lstr[1]
             elif "PRESSURE" in line or "pressure" in line:
                  self.pressure=lstr[1]
             elif "PRESSVAR" in line or "pressvar" in line:
                  self.pressvar=lstr[1]
             elif "OVERDEPBIAS" in line or "overdepbias" in line:
                  self.overdepbias=lstr[1]
             elif "CCDTEMP" in line or "ccdtemp" in line:
                  self.ccdtemp=lstr[1]
             elif "ALTVAR" in line or "altvar" in line:
                  self.altvar=lstr[1]
             elif "CONTROL" in line or "control" in line:
                  self.control=lstr[1]
             elif "ACTUATOR" in line:
                  self.actuator=line.split("ACTUATOR ")[1]
             elif "actuator" in line:
                  self.actuator=line.split("actuator ")[1]

        self.eventfile=0
        self.throughputfile=0
        self.centroidfile=0
        self.opdfile=0
        if extraCommands != 'none':
             for line in open(extraCommands):
                  lstr=line.split()
                  if "extraid" in line:
                       self.extraid=lstr[1]
                       self.observationID=self.observationID+self.extraid
                  if "eventfile" in line:
                       self.eventfile=int(float(lstr[1]))
                  if "throughputfile" in line:
                       self.throughputfile=int(float(lstr[1]))
                  if "centroidfile" in line:
                       self.centroidfile=int(float(lstr[1]))
                  if "opdfile" in line:
                       self.opdfile=int(float(lstr[1]))

     ## writeInputParamsAndCatalogs encapsulate the two instance catalog processing functions.
     def writeInputParamsAndCatalogs(self):
        self.writeInputParams()
        self.writeCatalogList()

     ## writeInputParams takes some of the parsed input parameters out of the
     # instance catalog and puts them in a single file.  We mainly have to read this
     # in and write out the same file because the instance catalog format cannot change
     #rapidly due to compatibility with opSim and catSim.
     def writeInputParams(self):
        self.inputParams='obs_'+self.observationID+'.pars'
        pfile=open(self.inputParams,'w')
        pfile.write("obshistid %s\n" % self.observationID)
        pfile.write("moonra %s\n" % self.moonra)
        pfile.write("moondec %s\n" % self.moondec)
        pfile.write("solaralt %s\n" % self.solaralt)
        pfile.write("moonalt %s\n" % self.moonalt)
        pfile.write("moondist %s\n" % self.moondist)
        pfile.write("phaseang %s\n" % self.phaseang)
        pfile.write("tai %s\n" % self.tai)
        pfile.write("constrainseeing %s\n" % self.constrainseeing)
        pfile.write("spiderangle %s\n" % self.spiderangle)
        pfile.write("azimuth %s\n" % self.azimuth)
        pfile.write("altitude %s\n" % self.altitude)
        pfile.write("altvar %s\n" % self.altvar)
        pfile.write("temperature %s\n" % self.temperature)
        pfile.write("tempvar %s\n" % self.tempvar)
        pfile.write("pressure %s\n" % self.pressure)
        pfile.write("pressvar %s\n" % self.pressvar)
        pfile.write("overdepbias %s\n" % self.overdepbias)
        pfile.write("ccdtemp %s\n" % self.ccdtemp)
        pfile.write("actuator %s" % self.actuator)
        pfile.write("control %s\n" % self.control)
        pfile.write("rotationangle %s\n" % self.rotationangle)
        pfile.write("pointingra %s\n" % self.pointingra)
        pfile.write("pointingdec %s\n" % self.pointingdec)
        pfile.write("obsseed %s\n" % self.obsseed)
        pfile.write("monthnum %s\n" % self.monthnum)
        pfile.write("filter %s\n" % self.filt)
        pfile.write("vistime %g\n" % self.vistime)
        pfile.write("camconfig %s\n"% self.camconfig)
        pfile.write("seddir %s\n" % self.sedDir)
        pfile.write("imagedir %s\n" % self.imageDir)
        pfile.write("datadir %s\n" % self.dataDir)
        pfile.write("instrdir %s\n" % self.instrDir)
        pfile.write("bindir %s\n" % self.binDir)
        pfile.write("telconfig %d\n" % self.telconfig)
        pfile.write("minsource %d\n" % self.minNumSources)
        pfile.write("domelight %g\n" % self.domeint)
        pfile.write("domewave %g\n" % self.domewav)
        if self.flatdir:
             pfile.write("flatdir 1\n")
        if self.tarfile:
             pfile.write("tarfile 1\n")
        pfile.close()

     ## writeCatalogList simply makes a list of possible
     # sub-catalogs (using the includeobj option) or lists of
     # objects put in the instance catalog.  The former is useful
     # for 1000s of objects, whereas the latter is useful for entire
     # focalplanes (millions).  Hence we support both of these options.
     def writeCatalogList(self):
        assert self.instanceCatalog
        assert self.userCatalog
        l=0
        objectCatalog=open('objectcatalog_'+self.observationID+'.pars','w')
        for line in self.userCatalog:
             if "object" in line:
                  objectCatalog.write(line)
                  l+=1
        objectCatalog.close()
        ncat=0
        catalogList=open('catlist_'+self.observationID+'.pars','w')
        if l>0:
             catalogList.write("catalog %d objectcatalog_%s.pars\n" % (ncat,self.observationID))
             ncat=1
        else:
             removeFile('objectcatalog_'+self.observationID+'.pars')
        catDir = os.path.dirname(self.instanceCatalog)
        for line in self.userCatalog:
            if "includeobj" in line:
                path = os.path.join(catDir, line.split()[1])
                if not os.path.isabs(catDir):
                    path = os.path.join("..", path)
                catalogList.write("catalog %d %s\n" % (ncat, path))
                ncat+=1
        catalogList.close()

     ## generateAtmosphere runs the atmosphere program
     def generateAtmosphere(self):
          assert self.inputParams
          assert self.extraCommands
          inputParams='obsExtra_'+self.observationID+'.pars'
          pfile=open(inputParams,'w')
          pfile.write(open(self.inputParams).read())
          if self.extraCommands!='none':
              pfile.write(open(self.extraCommands).read())
          pfile.close()
          runProgram("atmosphere < "+inputParams, self.binDir)
          removeFile(inputParams)

     ## generateInstrumentConfig runs the instrument program
     def generateInstrumentConfig(self):
          assert self.inputParams
          assert self.extraCommands
          inputParams='obsExtra_'+self.observationID+'.pars'
          pfile=open(inputParams,'w')
          pfile.write(open(self.inputParams).read())
          if self.extraCommands!='none':
             for line in open(self.extraCommands):
                  lstr=line.split()
                  if "dlsm" in line:
                       pfile.write(line)
          pfile.close()
          runProgram("instrument < "+inputParams, self.binDir)
          removeFile(inputParams)

     ## trimObjects runs the trim program
     #  Note this is overly complicated because we want to allow the trimming
     #  on grid computing to be done in groups to reduce the I/O of sending
     #  the entire instance catalog for every chip.  This complex looping
     #  isn't necessary for desktops.
     def trimObjects(self, sensors):

          self.initExecutionEnvironment()

          camstr="%03d" % int(float(bin(self.camconfig).split('b')[1]))
          if self.camconfig==0:
               camstr='111'
          fp=open(self.instrDir+"/focalplanelayout.txt").readlines()
          chipID=[]
          runFlag=[]
          devtype=[]
          devvalue=[]

          #Go through the focalplanelayout.txt filling up the arrays
          for line in fp:
               lstr=line.split()
               addFlag=0
               if "Group0" in line and camstr[2]=='1': addFlag=1
               elif "Group1" in line and camstr[1]=='1': addFlag=1
               elif "Group2" in line and camstr[0]=='1': addFlag=1
               if addFlag==1:
                    chipID.append(lstr[0])
                    runFlag.append(1)
                    devtype.append(lstr[6])
                    devvalue.append(float(lstr[7]))

          # See if we limit ourselves to a specific set of chipID (seperated by "|").
          if sensors != 'all':
               lstr = sensors.split('|')
               for i in range(len(chipID)): runFlag[i]=0
               for j in range(len(lstr)):
                    for i in range(len(chipID)):
                         if lstr[j]==chipID[i]:
                              runFlag[i]=1
                              break

          lastchip=chipID[-1]
          chipcounter1=0
          chipcounter2=0
          tc=0
          i=0
          trimJobID=[]
          for cid in chipID:
               if chipcounter1==0:
                    jobName='trim_'+self.observationID+'_'+str(tc)
                    inputParams=jobName+'.pars'
                    pfile=open(inputParams,'w')

               pfile.write('chipid %d %s\n' % (chipcounter1,cid))
               chipcounter1+=1
               if runFlag[i]==1:
                    chipcounter2+=1
               if chipcounter1==9 or cid==lastchip:   #Do groups of 9 to reduce grid computing I/O
                    trimJobID.append('none') 
                    pfile.write(open('obs_'+self.observationID+'.pars').read())
                    if self.flatdir:
                         for line in open('catlist_'+self.observationID+'.pars'):
                              lstr=line.split()
                              pfile.write('%s %s %s\n' % (lstr[0],lstr[1],lstr[2].split('/')[-1]))
                    else:
                         pfile.write(open('catlist_'+self.observationID+'.pars').read())
                    pfile.close()
                    if chipcounter2>0:
                         if self.grid in ['no', 'cluster']:
                              runProgram("trim < "+inputParams, self.binDir)
                         elif self.grid == 'condor':
                              nexp=self.nsnap if devtype[i]=='CCD' else int(self.vistime/devvalue[i])
                              condor.writeTrimDag(self,jobName,tc,nexp)
                         elif self.grid == 'diagrid': 
                              nexp=self.nsnap if devtype[i]=='CCD' else int(self.vistime/devvalue[i]) 
                              trimJobID[tc]=diagrid.writeTrimDag(self,jobName,tc,nexp) 
                         else:
                              sys.stderr.write('Unknown grid type: %s' % self.grid)
                              sys.exit(-1)
                    if self.grid in ['no', 'cluster'] or (self.grid in ['condor','diagrid'] and chipcounter2==0):
                         removeFile(inputParams)
                    chipcounter1=0
                    chipcounter2=0
                    tc+=1
               i=i+1
          self.chipID = chipID
          self.runFlag = runFlag
          self.devtype = devtype
          self.devvalue = devvalue
          self.trimJobID = trimJobID

     #scheduleRaytrace sets up the raytrace & e2adc jobs and also figures out the
     #numbers of exposures to perform.
     def scheduleRaytrace(self, instrument='lsst', run_e2adc=True, keep_screens=False):

        assert self.extraCommands
        chipcounter1=0
        tc=0
        counter=0
        jobs=[]
        i=0
        seg=open(self.instrDir+'/segmentation.txt').readlines()
        observationID = self.observationID

        for cid in self.chipID:
            if self.runFlag[i]==1:
                numSources=self.minNumSources
                if self.grid in ['no', 'cluster']:
                    numSources=len(open('trimcatalog_'+observationID+'_'+cid+'.pars').readlines())
                    numSources=numSources-2
                if numSources>=self.minNumSources:
                    nexp=self.nsnap if self.devtype[i]=='CCD' else int(self.vistime/self.devvalue[i])
                    ex=0
                    while ex<nexp:
                        eid="E%03d" % (ex)
                        fid=observationID + '_' + cid + '_' + eid
                        pfile=open('image_'+fid+'.pars','w')
                        pfile.write("chipid %s\n" % cid)
                        pfile.write("exposureid %d\n" % ex)
                        pfile.write("nsnap %d\n" % nexp)
                        pfile.close()

                        # PHOTON RAYTRACE
                        pfile=open('raytrace_'+fid+'.pars','w')
                        pfile.write(open('obs_'+observationID+'.pars').read())
                        pfile.write(open('atmosphere_'+observationID+'.pars').read())
                        pfile.write(open('optics_'+observationID+'.pars').read())
                        pfile.write(open('chip_'+observationID+'_'+cid+'.pars').read())
                        pfile.write(open('image_'+fid+'.pars').read())
                        if self.extraCommands!='none':
                            pfile.write(open(self.extraCommands).read())
                        if self.grid in ['no', 'cluster']:
                            pfile.write(open('trimcatalog_'+observationID+'_'+cid+'.pars').read())
                        pfile.close()

                        # ELECTRONS TO ADC CONVERTER
                        if run_e2adc:
                            pfile=open('e2adc_'+fid+'.pars','w')
                            pfile.write(open('obs_'+observationID+'.pars').read())
                            pfile.write(open('readout_'+observationID+'_'+cid+'.pars').read())
                            pfile.write(open('image_'+fid+'.pars').read())
                            pfile.close()

                        if self.grid == 'no':
                            p=multiprocessing.Process(target=jobChip,
                                                      args=(observationID,cid,eid,self.filt, self.outputDir,
                                                            self.binDir, self.instrDir),
                                                      kwargs={'instrument': instrument, 'run_e2adc': run_e2adc})
                            jobs.append(p)
                            p.start()
                            counter+=1
                            if counter==self.grid_opts.get('numproc', 1):
                                for p in jobs:
                                    p.join()
                                counter=0
                                jobs=[]
                        elif self.grid == 'cluster':
                            if self.grid_opts.get('script_writer', None):
                                self.grid_opts['script_writer'](observationID, cid, eid, self.filt,
                                                                   self.outputDir, self.binDir, self.dataDir)
                            else:
                                sys.stderr.write('WARNING: No script_writer callback in grid_opts for grid "cluster".\n')
                            if self.grid_opts.get('submitter', None):
                                self.grid_opts['submitter'](observationID, cid, eid)
                            else:
                                sys.stdout.write('No submitter callback in self.grid_opts for grid "cluster".\n')
                        elif self.grid == 'condor':
                            condor.writeRaytraceDag(self,cid,eid,tc,run_e2adc)
                        elif self.grid == 'diagrid': 
                            diagrid.writeRaytraceDag(self,cid,eid,tc,run_e2adc)

                        removeFile('image_'+fid+'.pars')
                        ex+=1

            chipcounter1+=1
            if chipcounter1==9:
                tc+=1
                chipcounter1=0

            if self.grid in ['no', 'cluster']:
                if os.path.exists('trimcatalog_'+observationID+'_'+cid+'.pars'):
                    removeFile('trimcatalog_'+observationID+'_'+cid+'.pars')
            removeFile('readout_'+observationID+'_'+cid+'.pars')
            removeFile('chip_'+observationID+'_'+cid+'.pars')
            i+=1

        removeFile('obs_'+observationID+'.pars')
        if not keep_screens:
             removeFile('atmosphere_'+observationID+'.pars')
        removeFile('optics_'+observationID+'.pars')
        removeFile('catlist_'+observationID+'.pars')

        if self.grid == 'no':
            for p in jobs:
                p.join()
        elif self.grid == 'condor':
            condor.submitDag(self)
        elif self.grid == 'diagrid': 
	    diagrid.submitDax(self)
        os.chdir(self.phosimDir)
        return

     ## Generic methods for handling execution environment
     def initExecutionEnvironment(self):
        if self.execEnvironmentInitialized:
            return
        if self.grid == 'condor':
            self.initCondorEnvironment()
        elif self.grid == 'diagrid': 
            self.initDiagridEnvironment() 
        elif self.grid == 'cluster':
            self.initClusterEnvironment()
        self.execEnvironmentInitialized = True

     ## general method to delete files at end
     def cleanup(self,keep_screens):
        if self.grid in ['no', 'cluster']:
            os.chdir(self.workDir)
            removeFile('objectcatalog_'+self.observationID+'.pars')
            removeFile('tracking_'+self.observationID+'.pars')
            if not keep_screens:
                 removeFile('airglowscreen_'+self.observationID+'.fits.gz')
                 for f in glob.glob('atmospherescreen_'+self.observationID+'_*') :
                      removeFile(f)
                 for f in glob.glob('cloudscreen_'+self.observationID+'_*') :
                      removeFile(f)
            else:
                 f='atmosphere_'+self.observationID+'.pars'
                 shutil.move(f,self.outputDir+'/'+f)
                 f='airglowscreen_'+self.observationID+'.fits.gz'
                 shutil.move(f,self.outputDir+'/'+f)
                 for f in glob.glob('atmospherescreen_'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
                 for f in glob.glob('cloudscreen_'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            if self.eventfile==1:
                 f='output.fits'
                 shutil.move(f,self.outputDir+'/'+f)
            if self.throughputfile==1:
                 for f in glob.glob('throughput_*'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            if self.centroidfile==1:
                 for f in glob.glob('centroid_*'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            if self.opdfile==1:
                 f='opd.fits'
                 shutil.move(f,self.outputDir+'/'+f)
	    	    
            os.chdir(self.phosimDir)

     ## Condor method to setup directories
     def initCondorEnvironment(self):
        sys.path.append(self.phosimDir+'/condor')
        global condor
        import condor
        condor.initEnvironment(self)

     def initDiagridEnvironment(self): 
        sys.path.append(self.phosimDir+'/diagrid') 
        global diagrid 
        import diagrid 
        diagrid.initEnvironment(self) 

     ## Cluster methods
     def initClusterEnvironment(self):
        pass

##main function
def main():

     phosimDir=os.path.split(os.path.abspath(__file__))[0]

     parser = optparse.OptionParser(usage='%prog instance_catalog [<arg1> <arg2> ...]')
     parser.add_option('-c','--command',dest="extraCommands",default="none")
     parser.add_option('-p','--proc',dest="numproc",default=1,type="int")
     parser.add_option('-o','--output',dest="outputDir",default=phosimDir+'/output')
     parser.add_option('-w','--work',dest="workDir",default=phosimDir+'/work')
     parser.add_option('-b','--bin',dest="binDir",default=phosimDir+'/bin')
     parser.add_option('-d','--data',dest="dataDir",default=phosimDir+'/data')
     parser.add_option('--sed',dest="sedDir",default=phosimDir+'/data/SEDs')
     parser.add_option('--image',dest="imageDir",default=phosimDir+'/data/images')
     parser.add_option('-s','--sensor',dest="sensor",default="all")
     parser.add_option('-i','--instrument',dest="instrument",default="lsst")
     parser.add_option('-g','--grid',dest="grid",default="no")
     parser.add_option('-u','--universe',dest="universe",default="standard")
     parser.add_option('-e','--e2adc',dest="e2adc",default=1,type="int")
     parser.add_option('--keepscreens',dest="keepscreens",default=0,type="int")
     parser.add_option('--checkpoint',dest="checkpoint",default=12,type="int")

     if len(sys.argv)<2:
        usage()
        sys.exit()

     if sys.argv[1] in ('-h', '--help'):
        usage()
        sys.exit()

     opt, remainder = parser.parse_args(sys.argv[1:]) #parse_args returns a pair of values
     instanceCatalog=remainder[0]

     checkPaths(opt, phosimDir)

     grid_opts = {'numproc': opt.numproc}
     if opt.grid == 'condor':
          grid_opts = {'universe': opt.universe, 'checkpoint': opt.checkpoint}
     elif opt.grid == 'diagrid': 
 	  grid_opts = {'checkpoint': opt.checkpoint} 
     elif opt.grid == 'cluster':
          grid_opts = {'script_writer': jobChip}

     #the entire phosim workflow follows:
     focalplane = PhosimFocalplane(phosimDir, opt, grid_opts)
     focalplane.loadInstanceCatalog(instanceCatalog, opt.extraCommands)
     os.chdir(opt.workDir)
     focalplane.writeInputParamsAndCatalogs()
     focalplane.generateAtmosphere()
     focalplane.generateInstrumentConfig()
     focalplane.trimObjects(opt.sensor)
     focalplane.scheduleRaytrace(opt.instrument, bool(opt.e2adc),bool(opt.keepscreens))
     focalplane.cleanup(bool(opt.keepscreens))

if __name__ == "__main__":
    main()
