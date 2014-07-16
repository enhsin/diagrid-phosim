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
    for f in ['focalplanelayout','central_wavelengths','segmentation']:
        fp=File(f+'.txt')
        fp.addPFN(PFN("file://" + os.path.join(self.instrDir,f+'.txt'), "local"))
        self.dax.addFile(fp)
    tarName='raytrace_'+observationID+'.tar'
    tar = tarfile.open(tarName, "w")
    tar.add(os.path.join(self.binDir,'version'),"version")
    for f in ['m1_protAl_Ideal','m2_protAl_Ideal','m3_protAl_Ideal','silica_dispersion','lenses',
              'detectorar','focalplanelayout','silicon','location','central_wavelengths','spider','tracking']:
        tar.add(os.path.join(self.instrDir,f+'.txt'),f+'.txt')
    for filt in range(6):
        tar.add(os.path.join(self.instrDir,'optics_'+str(filt)+'.txt'),'optics_'+str(filt)+'.txt')
        tar.add(os.path.join(self.instrDir,'filter_'+str(filt)+'.txt'),'filter_'+str(filt)+'.txt')
    for f in ['rayprofile','nprofile','o3cs','o3profile','o2cs','h2ocs','h2oprofile']:
        tar.add(os.path.join(self.dataDir, 'atmosphere', f+'.txt'),f+'.txt')
    for f in ['darksky_sed','lunar_sed','sed_dome','sersic_const']:
        tar.add(os.path.join(self.dataDir, 'sky', f+'.txt'),f+'.txt')
    for n in range(1,131):
        tar.add(os.path.join(self.dataDir, 'cosmic_rays', 'iray'+str(n)+'.txt'),'iray'+str(n)+'.txt')
    for f in ['cloudscreen_'+observationID+'_1.fits.gz', 'cloudscreen_'+observationID+'_2.fits.gz', 'airglowscreen_'+observationID+'.fits.gz']:
        tar.add(f)
        os.remove(f)
    for layer in range(7):
        for f in ['coarsep','coarsex','coarsey','fineh','finep',
                 'largep','largex','largey','mediumh','mediump','mediumx','mediumy']:
            tar.add('atmospherescreen_'+observationID+'_'+str(layer)+'_'+f+'.fits.gz')
            os.remove('atmospherescreen_'+observationID+'_'+str(layer)+'_'+f+'.fits.gz')
    tar.close()
    fp=File(tarName)
    fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),tarName), "local"))
    self.dax.addFile(fp)
    fp=File('tracking_'+observationID+'.pars')
    fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),'tracking_'+observationID+'.pars'), "local"))
    self.dax.addFile(fp)
    for line in open('catlist_'+observationID+'.pars'):
        objectCatalog=line.split()[2]
        fp=File(objectCatalog)
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),objectCatalog), "local"))
        self.dax.addFile(fp)

    e_trim = Executable(namespace="phosim", name="trim", os="linux", arch="x86_64", installed=False)
    e_trim.addPFN(PFN("file://" + os.path.join(self.binDir,'trim'), "condorpool"))
    self.dax.addExecutable(e_trim)

    e_raytrace = Executable(namespace="phosim", name="raytrace", os="linux", arch="x86_64", installed=False)
    e_raytrace.addPFN(PFN("file://" + os.path.join(self.binDir,'raytrace'), "condorpool"))
    e_raytrace.addProfile(Profile(namespace="condor", key="requirements", value='TARGET.ClusterName == "Hansen"'))
    self.dax.addExecutable(e_raytrace)

    e_e2adc = Executable(namespace="phosim", name="e2adc", os="linux", arch="x86_64", installed=False)
    e_e2adc.addPFN(PFN("file://" + os.path.join(self.binDir,'e2adc'), "condorpool"))
    e_e2adc.addProfile(Profile(namespace="condor", key="requirements", value='TARGET.ClusterName == "Hansen"'))
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

    for line in open(jobName+'.pars'):
        if "chipid" in line:
            cid=line.split()[2]
            trim.uses( File('trimcatalog_'+self.observationID+'_'+cid+'.pars'), link=Link.OUTPUT)

    trim.addProfile(Profile(namespace="dagman", key="POST", value="posttrim"))
    trim.addProfile(Profile(namespace="dagman", key="POST.PATH.posttrim", value=os.path.join(self.binDir,"diagrid","chip")))
    arg='posttrim %s %d %d %s %d %s' % (self.observationID,tc,nexp,self.sedDir,checkpoint,self.workDir)
    trim.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))
    self.dax.addJob(trim)
    return jobID

def writeRaytraceDag(self,cid,eid,tc,run_e2adc):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    observationID=self.observationID
    fid=observationID + '_' + cid + '_' + eid
    fidfilt=observationID + '_f'+self.filt+'_' + cid + '_' + eid
    instrument=self.instrDir.split("/")[-1]
    for ckpt in range(checkpoint+1):
        fidckpt=fid+'_'+str(ckpt)
        jobName='raytrace_'+fidckpt
        jobID=self.dax.nextJobID()
        exec('raytrace'+str(ckpt)+'=Job(namespace="phosim", name="raytrace", id=jobID)')
        eval('raytrace'+str(ckpt)+'.setStdin("'+jobName+'.pars")')
        eval('raytrace'+str(ckpt)+'.uses(File("'+jobName+'.pars"), link=Link.INPUT)')
        fp=File(jobName+'.pars')
        fp.addPFN(PFN("file://" + os.path.join(os.getcwd(),jobName+'.pars'), "local"))
        self.dax.addFile(fp)
        eval('raytrace'+str(ckpt)+'.uses(File("tracking_'+observationID+'.pars"), link=Link.INPUT)')
        eval('raytrace'+str(ckpt)+'.uses(File("raytrace_'+observationID+'.tar"), link=Link.INPUT)')
        eval('raytrace'+str(ckpt)+'.uses(File("segmentation.txt"), link=Link.INPUT)') # will replace with SEDs
        if ckpt>0:
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fidfilt+'_ckptdt_'+str(ckpt-1)+'.fits.gz"), link=Link.INPUT)')
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fidfilt+'_ckptfp_'+str(ckpt-1)+'.fits.gz"), link=Link.INPUT)')
        if ckpt<checkpoint:
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fidfilt+'_ckptdt_'+str(ckpt)+'.fits.gz"), link=Link.OUTPUT, transfer=False, register=False)')
            eval('raytrace'+str(ckpt)+'.uses(File("'+instrument+'_e_'+fidfilt+'_ckptfp_'+str(ckpt)+'.fits.gz"), link=Link.OUTPUT, transfer=False, register=False)')
        if ckpt==checkpoint:
            fileName=instrument+'_e_'+fidfilt+'.fits.gz'
            if run_e2adc:
                eval('raytrace'+str(ckpt)+'.uses(File(fileName), link=Link.OUTPUT, transfer=False, register=False)')
            else:
                eval('raytrace'+str(ckpt)+'.uses(File(fileName), link=Link.OUTPUT)')
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST", value="postraytrace"))')
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST.PATH.postraytrace", value=os.path.join(self.binDir,"diagrid","chip")))')
                arg='postraytrace %s %s %s %s %s %s' % (observationID,self.filt,cid,eid,self.workDir,instrument)
                eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))')

        eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="PRE", value=os.path.join(self.binDir,"diagrid","chip")))')
        arg='preraytrace %s %s %s' % (observationID,jobID,cid)
        eval('raytrace'+str(ckpt)+'.addProfile(Profile(namespace="dagman", key="PRE.ARGUMENTS", value=arg))')
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
        e2adc.uses( File('%s_e_%s.fits.gz' % (instrument,fidfilt)), link=Link.INPUT)
        e2adc.uses(File(instrument+'_'+fidfilt+'.tar'), link=Link.OUTPUT)
        e2adc.addProfile(Profile(namespace="dagman", key="POST", value="poste2adc"))
        e2adc.addProfile(Profile(namespace="dagman", key="POST.PATH.poste2adc", value=os.path.join(self.binDir,"diagrid","chip")))
        arg='poste2adc %s %s %s %s %s %s' % (observationID,self.filt,cid,eid,self.workDir,instrument)
        e2adc.addProfile(Profile(namespace="dagman", key="POST.ARGUMENTS", value=arg))
        self.dax.addJob(e2adc)
        eval('self.dax.depends(parent=raytrace'+str(checkpoint)+',child=e2adc)')

    os.remove('raytrace_'+fid+'.pars')

def submitDax(self):
    fp=open('phosim_'+self.observationID+'.dax','w')
    self.dax.writeXML(fp)
    fp.close()

    command='submit pegasus-plan --dax phosim_'+self.observationID+'.dax'
    if subprocess.call(command, shell=True) != 0:
        raise RuntimeError("Error running %s" % command)

