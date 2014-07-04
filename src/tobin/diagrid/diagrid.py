import sys
import os
import subprocess
import tarfile
from Pegasus.DAX3 import *

def initEnvironment(self):
    # Create a abstract dag
    self.dax=ADAG('phosim')

    observationID=self.observationID
    # Add input file to the DAX-level replica catalog
    fp=File("version")
    fp.addPFN(PFN("file://" + os.path.join(self.binDir,'version'), "local"))
    self.dax.addFile(fp)
    for f in ['m1_protAl_Ideal','m2_protAl_Ideal','m3_protAl_Ideal','silica_dispersion','lenses',
              'detectorar','focalplanelayout','silicon','location','central_wavelengths','spider',
              'tracking','segmentation']:
        fp=File(f+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.instrDir,f+'.txt'), "local"))
        self.dax.addFile(fp)
    for filt in range(6):
        fp=File('optics_'+str(filt)+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.instrDir,'optics_'+str(filt)+'.txt'), "local"))
        self.dax.addFile(fp)
        fp=File('filter_'+str(filt)+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.instrDir,'filter_'+str(filt)+'.txt'), "local"))
        self.dax.addFile(fp)
    for f in ['rayprofile','nprofile','o3cs','o3profile','o2cs','h2ocs','h2oprofile']:
        fp=File(f+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.dataDir, 'atmosphere', f+'.txt'), "local"))
        self.dax.addFile(fp)
    for f in ['darksky_sed','lunar_sed','sed_dome','sersic_const']:
        fp=File(f+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.dataDir, 'sky', f+'.txt'), "local"))
        self.dax.addFile(fp)
    for n in range(1,131):
        fp=File('iray'+str(n)+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.dataDir, 'cosmic_rays', 'iray'+str(n)+'.txt'), "local"))
        self.dax.addFile(fp)
    for f in ['cloudscreen_'+observationID+'_1.fits', 'cloudscreen_'+observationID+'_2.fits', 'airglowscreen_'+observationID+'.fits']:
        fp=File(f)
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),f), "local"))
        self.dax.addFile(fp)
    for layer in range(7):
        for f in ['coarsep','coarsex','coarsey','fineh','finep',
                 'largep','largex','largey','mediumh','mediump','mediumx','mediumy']:
            fp=File('atmospherescreen_'+observationID+'_'+str(layer)+'_'+f+'.fits')
            fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),'atmospherescreen_'+observationID+'_'+str(layer)+'_'+f+'.fits'), "local"))
            self.dax.addFile(fp)
    fp=File('tracking_'+observationID+'.pars')
    fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),'tracking_'+observationID+'.pars'), "local"))
    self.dax.addFile(fp)

    e_trim = Executable(namespace="phosim", name="trim", os="linux", arch="x86_64", installed=False)
    e_trim.addPFN(PFN("file://" + os.path.join(self.binDir,'trim'), "condorpool"))
    self.dax.addExecutable(e_trim)

    e_raytrace = Executable(namespace="phosim", name="raytrace", os="linux", arch="x86_64", installed=False)
    e_raytrace.addPFN(PFN("file://" + os.path.join(self.binDir,'raytrace'), "condorpool"))
    self.dax.addExecutable(e_raytrace)

    e_e2adc = Executable(namespace="phosim", name="e2adc", os="linux", arch="x86_64", installed=False)
    e_e2adc.addPFN(PFN("file://" + os.path.join(self.binDir,'e2adc'), "condorpool"))
    self.dax.addExecutable(e_e2adc)

def writeTrimDag(self,jobName,tc,nexp):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    jobID=self.dax.nextJobID()
    trim = Job(namespace="phosim", name="trim", id=jobID)
    trim.setStdin(jobName+'.pars')
    trim.uses( File(jobName+'.pars'), link=Link.INPUT)
    trim.uses( File("focalplanelayout.txt"), link=Link.INPUT)
    trim.uses( File("central_wavelengths.txt"), link=Link.INPUT)
    fp=File(jobName+'.pars')
    fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),jobName+'.pars'), "local"))
    self.dax.addFile(fp)
    for line in open('catlist_'+self.observationID+'.pars'):
        objectCatalog=line.split()[2]
        trim.uses( File(objectCatalog), link=Link.INPUT)
        fp=File(objectCatalog)
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),objectCatalog), "local"))
        self.dax.addFile(fp)

    for line in open(jobName+'.pars'):
        if "chipid" in line:
            cid=line.split()[2]
            trim.uses( File('trimcatalog_'+self.observationID+'_'+cid+'.pars'), link=Link.OUTPUT) 
    
    trim.addProfile(Profile(namespace="dagman", key="POST", value="posttrim"))
    trim.addProfile(Profile(namespace="dagman", key="POST.PATH.posttrim", value=os.path.join(self.binDir,"diagrid","chip")))
    arg='posttrim %s %d %d %s %d' % (self.observationID,tc,nexp,self.sedDir,checkpoint)
    trim.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))
    self.dax.addJob(trim)
    return jobID

def writeRaytraceDag(self,cid,eid,tc,run_e2adc):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    observationID=self.observationID
    fid=observationID + '_' + cid + '_' + eid
    instrument=self.instrDir.split("/")[-1]
    for ckpt in range(checkpoint+1):
        fidckpt=fid+'_'+str(ckpt)
        jobName='raytrace_'+fidckpt
        exec('raytrace'+str(ckpt)+'=Job(namespace="phosim", name="raytrace")')
        eval('raytrace'+str(ckpt)+'.setStdin("'+jobName+'.pars")')
        eval('raytrace'+str(ckpt)+'.uses(File("'+jobName+'.pars"), link=Link.INPUT)')
        fp=File(jobName+'.pars')
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),jobName+'.pars'), "local"))
        self.dax.addFile(fp)
        for f in ['cloudscreen_'+observationID+'_1.fits', 'cloudscreen_'+observationID+'_2.fits', 'airglowscreen_'+observationID+'.fits']:
            eval('raytrace'+str(ckpt)+'.uses(File("'+f+'"), link=Link.INPUT)')
        for layer in range(7):
            for f in ['coarsep','coarsex','coarsey','fineh','finep',
                     'largep','largex','largey','mediumh','mediump','mediumx','mediumy']:
                eval('raytrace'+str(ckpt)+'.uses(File("atmospherescreen_'+observationID+'_'+str(layer)+'_'+f+'.fits"), link=Link.INPUT)')
        eval('raytrace'+str(ckpt)+'.uses(File("version"), link=Link.INPUT)')
        for f in ['m1_protAl_Ideal','m2_protAl_Ideal','m3_protAl_Ideal','silica_dispersion','lenses',
                  'detectorar','focalplanelayout','silicon','location','central_wavelengths','spider','tracking']:
            eval('raytrace'+str(ckpt)+'.uses(File("'+f+'.txt"), link=Link.INPUT)')
        for filt in range(6):
            eval('raytrace'+str(ckpt)+'.uses(File("optics_'+str(filt)+'.txt"), link=Link.INPUT)')
            eval('raytrace'+str(ckpt)+'.uses(File("filter_'+str(filt)+'.txt"), link=Link.INPUT)')
        for f in ['rayprofile','nprofile','o3cs','o3profile','o2cs','h2ocs','h2oprofile']:
            eval('raytrace'+str(ckpt)+'.uses(File("'+f+'.txt"), link=Link.INPUT)')
        for f in ['darksky_sed','lunar_sed','sed_dome','sersic_const']:
            eval('raytrace'+str(ckpt)+'.uses(File("'+f+'.txt"), link=Link.INPUT)')
        for n in range(1,131):
            eval('raytrace'+str(ckpt)+'.uses(File("iray'+str(n)+'.txt"), link=Link.INPUT)')
        eval('raytrace'+str(ckpt)+'.uses(File("tracking_'+observationID+'.pars"), link=Link.INPUT)')
        eval('raytrace'+str(ckpt)+'.uses(File("SEDs_"+observationID+"_"+cid+".tar"), link=Link.INPUT)')
        if ckpt==0:
            if eid=='E000':
                tarName='SEDs_'+observationID+'_'+cid+'.tar'
                tar = tarfile.open(tarName, "w")
                fileName=os.path.join(self.dataDir, 'sky','sed_flat.txt')
                tar.add(fileName,"tmp.txt")
                tar.close()
                fp=File(tarName)
                fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),tarName), "local"))
                self.dax.addFile(fp)
            if ckpt!=checkpoint:
                eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fid+'_ckptdt.fits"), link=Link.OUTPUT, transfer=False, register=False)')
                eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fid+'_ckptfp.fits"), link=Link.OUTPUT, transfer=False, register=False)')
        else:
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fid+'_ckptdt.fits"), link=Link.INPUT)')
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fid+'_ckptfp.fits"), link=Link.INPUT)')
        if ckpt==checkpoint:
            fileName=instrument+'_e_'+fid+'.fits.gz'
            eval('raytrace'+str(ckpt)+'.uses(File(fileName), link=Link.OUTPUT)')
            if not run_e2adc:
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST", value="postraytrace"))')
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST.PATH.postraytrace", value=os.path.join(self.binDir,"diagrid","chip")))')
                arg='postraytrace %s %s %s %s %s' % (observationID,self.filt,cid,eid,self.outputDir)
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))')

        eval('self.dax.addJob(raytrace'+str(ckpt)+')')

        if ckpt==0:
            trim=self.dax.getJob(self.trimJobID[tc])
            eval('self.dax.depends(parent=trim,child=raytrace'+str(ckpt)+')') 
        else:
            eval('self.dax.depends(parent=raytrace'+str(ckpt-1)+',child=raytrace'+str(ckpt)+')')
            if ckpt>1:
                eval('self.dax.depends(parent=raytrace0,child=raytrace'+str(ckpt)+')')

        pfile=open('raytrace_'+fidckpt+'.pars','w')
        pfile.write(open('raytrace_'+fid+'.pars').read())
        pfile.write("checkpointcount %d\n" % ckpt)
        pfile.write("checkpointtotal %d\n" % checkpoint)
        pfile.close()
    if run_e2adc:
        jobName='e2adc_'+fid
        e2adc = Job(namespace="phosim", name="e2adc")
        e2adc.setStdin(jobName+'.pars')
        e2adc.uses( File(jobName+'.pars'), link=Link.INPUT)
        fp=File(jobName+'.pars')
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),jobName+'.pars'), "local"))
        self.dax.addFile(fp)
        e2adc.uses( File('segmentation.txt'), link=Link.INPUT)
        e2adc.uses( File('focalplanelayout.txt'), link=Link.INPUT)
        e2adc.uses( File('%s_e_%s.fits.gz' % (instrument,fid)), link=Link.INPUT)
        segfile=os.path.join(self.instrDir,'segmentation.txt')
        for line in open(segfile):
            aid=line.split()[0]
            if cid in line and aid != cid:
                e2adc.uses(File(instrument+'_a_'+observationID+'_'+aid+'_'+eid+'.fits.gz'), link=Link.OUTPUT)
        fileName=instrument+'_'+observationID+'_f'+self.filt+'_'+cid+'_'+eid+'.tar'
        e2adc.addProfile(Profile(namespace="dagman", key="POST", value="poste2adc"))
        e2adc.addProfile(Profile(namespace="dagman", key="POST.PATH.poste2adc", value=os.path.join(self.binDir,"diagrid","chip")))
        arg='poste2adc %s %s %s %s %s %s' % (observationID,self.filt,cid,eid,self.outputDir,self.instrDir)
        e2adc.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))
        self.dax.addJob(e2adc)
        eval('self.dax.depends(parent=raytrace'+str(checkpoint)+',child=e2adc)')


def submitDax(self):
    fp=open('phosim_'+self.observationID+'.dax','w')
    self.dax.writeXML(fp)
    fp.close()

    command='submit pegasus-plan --dax phosim_'+self.observationID+'.dax'
    #print command
    if subprocess.call(command, shell=True) != 0:
        raise RuntimeError("Error running %s" % command)

