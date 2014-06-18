import os
import subprocess

## Condor method to setup directories
def initEnvironment(self):
    self.dagfile=open(self.workDir+'/dag_'+self.observationID+'.dag','w')
    if not os.path.exists(self.workDir+'/logs'):
        os.makedirs(self.workDir+'/logs')
    if not os.path.exists(self.workDir+'/errors'):
        os.makedirs(self.workDir+'/errors')
    if not os.path.exists(self.workDir+'/output'):
        os.makedirs(self.workDir+'/output')

## Condor method to submit job
def writeSubmit(self, job, jobName, fid='none', ckpt=0):
    if self.grid != 'condor':
        return
    assert job in ['raytrace', 'trim', 'e2adc']
    assert 'universe' in self.grid_opts
    universe = self.grid_opts['universe']
    submitfile=open(jobName+'.submit','w')
    submitfile.write('executable = %s/%s\n' % (self.binDir,job))
    submitfile.write('initialdir = %s\n' % self.workDir)
    submitfile.write('Universe = %s\n' % universe)
    submitfile.write('Input = %s.pars\n' % jobName)
    submitfile.write('Log = logs/log_%s.log\n' %jobName)
    submitfile.write('Output = output/out_%s.out\n' %jobName)
    submitfile.write('Error = errors/error_%s.error\n' %jobName)
    if universe == 'vanilla':
        submitfile.write('periodic_release=true\n')
        submitfile.write('should_transfer_files = YES\n')
    else:
        submitfile.write('should_transfer_files = IF_NEEDED\n')
    submitfile.write('when_to_transfer_output = ON_EXIT\n')
    submitfile.write('notification = NEVER\n')

    instrument=self.instrDir.split("/")[-1]
    if job == 'raytrace':
        submitfile.write('transfer_input_files = \\\n')
        submitfile.write('cloudscreen_%s_1.fits, \\\n' % self.observationID)
        submitfile.write('cloudscreen_%s_2.fits, \\\n' % self.observationID)
        submitfile.write('airglowscreen_%s.fits, \\\n' % self.observationID)
        for layer in range(7):
             for f in ['coarsep','coarsex','coarsey','fineh','finep',
                      'largep','largex','largey','mediumh','mediump','mediumx','mediumy']:
                 submitfile.write('atmospherescreen_%s_%d_%s.fits, \\\n' % (self.observationID,layer,f))
        submitfile.write('%s/version, \\\n' % self.binDir)
        for f in ['m1_protAl_Ideal','m2_protAl_Ideal','m3_protAl_Ideal','silica_dispersion','lenses',
                  'detectorar','focalplanelayout','silicon','location','central_wavelengths','spider','tracking']:
            submitfile.write('%s/%s.txt,    \\\n' % (self.instrDir,f))
        for filt in range(6):
            submitfile.write('%s/optics_%d.txt, \\\n' % (self.instrDir,filt))
            submitfile.write('%s/filter_%d.txt, \\\n' % (self.instrDir,filt))
        for f in ['rayprofile','nprofile','o3cs','o3profile','o2cs','h2ocs','h2oprofile']:
            submitfile.write('%s/atmosphere/%s.txt, \\\n' % (self.dataDir,f))
        for f in ['darksky_sed','lunar_sed','sed_dome','sersic_const']:
            submitfile.write('%s/sky/%s.txt, \\\n' % (self.dataDir,f))
        for n in range(1,131):
            submitfile.write('%s/cosmic_rays/iray%d.txt, \\\n' % (self.dataDir,n))
        submitfile.write('tracking_%s.pars' % self.observationID)

        if ckpt>0:
            submitfile.write(', %s_e_%s_ckptdt.fits' % (instrument,fid))
            submitfile.write(', %s_e_%s_ckptfp.fits' % (instrument,fid))

    else:
        if job == 'trim':
            submitfile.write('transfer_input_files =  %s/focalplanelayout.txt' % self.instrDir)
            submitfile.write(', %s/central_wavelengths.txt' % self.instrDir)
            for line in open('catlist_'+self.observationID+'.pars'):
                submitfile.write(', %s' % line.split()[2])
            submitfile.write('\n')
        elif job == 'e2adc':
            submitfile.write('transfer_input_files = %s/segmentation.txt, ' % self.instrDir)
            submitfile.write('%s/focalplanelayout.txt, ' % self.instrDir)
            submitfile.write('%s_e_%s.fits.gz\n' % (instrument,fid))
        submitfile.write('Queue 1\n')

    submitfile.close()

def writeTrimDag(self,jobName,tc,nexp):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    self.dagfile.write('JOB %s %s/%s.submit\n' % (jobName,self.workDir,jobName))
    self.dagfile.write('SCRIPT POST %s %s/condor/chip posttrim %s %d %d %s %s %d\n' %
                       (jobName,self.phosimDir,self.observationID,tc,nexp,self.dataDir,self.workDir,checkpoint))
    self.dagfile.write('RETRY %s 3\n' % (jobName))
    writeSubmit(self,'trim',jobName)

def writeRaytraceDag(self,cid,eid,tc,run_e2adc):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    observationID=self.observationID
    fid=observationID + '_' + cid + '_' + eid
    for ckpt in range(checkpoint+1):
        fidckpt=fid+'_'+str(ckpt)
        self.dagfile.write('JOB raytrace_%s %s/raytrace_%s.submit\n' % (fidckpt,self.workDir,fidckpt))
        self.dagfile.write('RETRY raytrace_%s 3\n' % (fidckpt))
        if ckpt==0:
            self.dagfile.write('SCRIPT PRE raytrace_%s_0 %s/condor/chip preraytrace %s %s\n' %
                          (fid,self.phosimDir,observationID+'_'+str(tc),self.workDir))
            self.dagfile.write('PARENT trim_%s_%d CHILD raytrace_%s_0\n' % (observationID,tc,fid))
        else:
            self.dagfile.write('PARENT raytrace_%s_%d CHILD raytrace_%s\n' % (fid,ckpt-1,fidckpt))
        writeSubmit(self,'raytrace','raytrace_'+fidckpt,fid,ckpt)
        pfile=open('raytrace_'+fidckpt+'.pars','w')
        pfile.write(open('raytrace_'+fid+'.pars').read())
        pfile.write("checkpointcount %d\n" % ckpt)
        pfile.write("checkpointtotal %d\n" % checkpoint)
        pfile.close()
    if run_e2adc:
        self.dagfile.write('JOB e2adc_%s %s/e2adc_%s.submit\n' % (fid,self.workDir,fid))
        self.dagfile.write('RETRY e2adc_%s 3\n' % fid)
        self.dagfile.write('SCRIPT PRE e2adc_%s %s/condor/chip pree2adc %s %s\n' % (fid,self.phosimDir,fid,self.workDir))
        self.dagfile.write('SCRIPT POST e2adc_%s %s/condor/chip poste2adc %s %s %s %s %s %s %s\n' %
                      (fid,self.phosimDir,observationID,self.filt,cid,eid,self.outputDir,self.instrDir,self.workDir))
        self.dagfile.write('PARENT raytrace_%s_%d CHILD e2adc_%s\n' % (fid,checkpoint,fid))
        writeSubmit(self,'e2adc','e2adc_'+fid,fid)
    else:
        self.dagfile.write('SCRIPT POST raytrace_%s_%d %s/condor/chip postraytrace %s %s %s %s %s %s\n' %
                      (fid,checkpoint,self.phosimDir,observationID,self.filt,cid,eid,self.outputDir,self.workDir))

def submitDag(self):
    self.dagfile.close()
    command='condor_submit_dag dag_'+self.observationID+'.dag'
    if subprocess.call(command, shell=True) != 0:
        raise RuntimeError("Error running %s" % command)

