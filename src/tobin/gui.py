#!/usr/bin/env python

from Tkinter import *
import ttk

import filebrowser
import time
import functools
import subprocess
import module_locator
import os
from collections import OrderedDict

workspacedir = os.path.join(os.environ['HOME'], "phosimrun")
if not os.path.exists(workspacedir):
	os.makedirs(workspacedir)
os.chdir(workspacedir)
currentdir = os.getcwd()
installdir = module_locator.module_path()
if os.path.basename(installdir) == "tobin":
	installdir = os.path.dirname(installdir)
phosimdir = os.path.dirname(installdir)
exampledir = os.path.join(phosimdir, "examples")
bindir = os.path.join(phosimdir, "bin")
datadir = os.path.join(phosimdir, "data")
seddir = os.path.join(datadir, "SEDs")
imagedir = os.path.join(datadir, "images")
phosimpy = os.path.join(bindir, "phosim.py")
validationdir = os.path.join(phosimdir, "validation")
validationscriptsdir = os.path.join(validationdir, "scripts")
#outputdir = os.path.join(currentdir, "output")
#workdir = os.path.join(currentdir, "work")
outputdir = currentdir
workdir = currentdir
calibrationfolders = [os.path.join(exampledir, "biases"), os.path.join(exampledir, "darks"), os.path.join(exampledir, "flats")]
if not os.path.exists(outputdir):
	os.mkdir(outputdir)
if not os.path.exists(workdir):
	os.mkdir(workdir)

browserapp = filebrowser.browserapp()

INSTANCECATALOG = "instancecatalog"
EXTRACOMMANDS = "extracommands"
E2ADC = "e2adc"
INSTRUMENT = "instrument"
GRID = "grid"
options = {INSTANCECATALOG:"", EXTRACOMMANDS:"none", E2ADC:"0", INSTRUMENT:"lsst", GRID:"no"}


def openmethod(paths, optionsindex=None, label=None):
	if len(paths) != 0:
		path = paths[0]
		name = os.path.basename(paths[0])
		label.config(text=name)
		options[optionsindex] = path

def launchworkspace():
	#browserapp.threadwindow(basefolder=[workspacedir, outputdir], basefoldername=["Workspace", "Output"], basefoldermakedirs=True, openmode=False)
	browserapp.threadwindow(basefolder=[workspacedir], basefoldername=["Workspace"], basefoldermakedirs=True, openmode=False)

def launchopen(optionsindex, label):
	browserapp.threadwindow(basefolder=[workspacedir, exampledir, validationdir, phosimdir], basefoldername=["Workspace", "Examples", "Validation", "Phosim"], basefoldermakedirs=True, openmode=True, openmethod=functools.partial(openmethod, optionsindex=optionsindex, label=label))

def launchinstrumentopen(optionsindex, label):
	browserapp.threadwindow(basefolder=[datadir, workspacedir, exampledir, validationdir, phosimdir], basefoldername=["Data", "Workspace", "Examples", "Validation", "Phosim"], basefoldermakedirs=True, openmode=True, openmethod=functools.partial(openmethod, optionsindex=optionsindex, label=label), foldersonly=True)

def updatee2adc():
	options[E2ADC] = str(e2adcvar.get())

def updategrid():
	options[GRID] = gridVar.get()

def phosimpyothercommands():
	return ["-o", outputdir, "-w", workdir, "-b", bindir, 
		"-d", datadir, "--sed="+seddir, "--image="+imagedir]

def runphosimpy(instancecatalog, extracommands, e2adcflag, grid="no", ckpt=0, instrument="lsst", sensor="all"):
	
	command = ["xterm", "-hold", "-title", "phosimrun", "-e", "python", phosimpy, instancecatalog, "-c", extracommands, 
		"-e", e2adcflag, "-g", grid, "--checkpoint="+str(ckpt), "-i", instrument,
		"-s", sensor] + phosimpyothercommands()
	
	print command
	
	result = subprocess.Popen(command)

def runuser():
	ckptnum=checkpointVar.get()
	sensorstr=sensorVar.get()
	runphosimpy(instancecatalog=options[INSTANCECATALOG], extracommands=options[EXTRACOMMANDS], e2adcflag=options[E2ADC],
		instrument=options[INSTRUMENT],sensor=sensorstr,grid=options[GRID],ckpt=ckptnum)
	
def runcalibration():
	cats = calibrationinstancecatalogs[calibrationinstancecatalog.get()]
	sets = numsets.get()
	
	print sets
	
	commands = []
	for cat in cats:
		commands += cat[:sets]
	
	for cat in commands:
		if "flat" in cat:
			runphosimpy(instancecatalog=cat, extracommands="none", e2adcflag="1", grid="diagrid", ckpt=54)
		else:
			runphosimpy(instancecatalog=cat, extracommands="none", e2adcflag="1", grid="diagrid", ckpt=0)

def runvalidation():
	script = vscripts[vscript.get()]
	
	command = command = ["xterm", "-hold", "-title", "phosimrun", "-e", script] + phosimpyothercommands()
		
	result = subprocess.Popen(command, cwd=phosimdir)

root = Tk()
root.wm_title("Phosim")

workspacebutton = Button(root, text="Open Folders", command=launchworkspace)
workspacebutton.grid(row=0, column=0)

notebook = ttk.Notebook(root)

## User
userframe = Frame(notebook)

instancecataloglabel = Label(userframe, text="Instance Catalog:")
instancecataloglabel.grid(row=0, column=0, sticky=W)

instancecatalogdisplay = Label(userframe, text="[no file selected]")
instancecatalogdisplay.grid(row=0, column=1, sticky=W)

instancecatalogbutton = Button(userframe, text="[...]", command=functools.partial(launchopen, INSTANCECATALOG, instancecatalogdisplay))
instancecatalogbutton.grid(row=0, column=2)

settinglabel = Label(userframe, text="Optional Setting:", font="Verdana 10 bold")
settinglabel.grid(row=1, column=0, sticky=W)

extracommandslabel = Label(userframe, text="Physics Override Commands:")
extracommandslabel.grid(row=2, column=0, sticky=W)

extracommandsdisplay = Label(userframe, text="none")
extracommandsdisplay.grid(row=2, column=1, sticky=W)

extracommandsbutton = Button(userframe, text="[...]", command=functools.partial(launchopen, EXTRACOMMANDS, extracommandsdisplay))
extracommandsbutton.grid(row=2, column=2)

instrumentlabel = Label(userframe, text="Instrument Site Directory:")
instrumentlabel.grid(row=3, column=0, sticky=W)

instrumentdisplay = Label(userframe, text="lsst")
instrumentdisplay.grid(row=3, column=1, sticky=W)

instrumentbutton = Button(userframe, text="[...]", command=functools.partial(launchinstrumentopen, INSTRUMENT, instrumentdisplay))
instrumentbutton.grid(row=3, column=2)

sensorlabel = Label(userframe, text="Sensor:")
sensorlabel.grid(row=4, column=0, sticky=W)
sensorVar = StringVar()
sensorVar.set("all")
sensorentry = Entry(userframe, text="all", textvariable=sensorVar)
sensorentry.grid(row=4, column=1)

e2adcvar = IntVar()

e2adccheckbox = Checkbutton(userframe, text="Run e2adc", command=updatee2adc, variable=e2adcvar)
e2adccheckbox.grid(row=5, column=0, columnspan=2, sticky=W)


gridlabel = Label(userframe, text="Execution Site:")
gridlabel.grid(row=6, column=0, sticky=W)
gridVar = StringVar()
gridVar.set("no")
gridbutton1 = Radiobutton(userframe, text="Local", variable=gridVar, value="no", command=updategrid)
gridbutton2 = Radiobutton(userframe, text="DiaGrid", variable=gridVar, value="diagrid", command=updategrid)
gridbutton1.grid(row=7, column=0, columnspan=1, sticky=W)
gridbutton2.grid(row=8, column=0, columnspan=1, sticky=W)
checkpointlabel = Label(userframe, text="Number of Checkpoints: ")
checkpointlabel.grid(row=8, column=1, sticky=W)
checkpointVar = IntVar()
checkpointVar.set(0)
checkpointspinner = Spinbox(userframe, from_=0, to=100, textvariable=checkpointVar)
checkpointspinner.grid(row=8, column=2, sticky=W)


outputbutton = Button(userframe, text="Run Simulation", command=runuser)
outputbutton.grid(row=9, column=0)

notebook.add(userframe, text="User")

## Calibration
calibrationframe = Frame(notebook)

calibrationinstancecatalogs = OrderedDict()
for calibrationfolder in calibrationfolders:
	cats = sorted(os.listdir(calibrationfolder))
	for instancecatalog in cats:
		if instancecatalog.endswith("_instcat_0"):
			instancecatalogname = instancecatalog.replace("_instcat_0", "")
			calibrationinstancecatalogs[instancecatalogname] = [[os.path.join(calibrationfolder, cat) for cat in cats if cat.startswith(instancecatalogname)]]
allcats = []
for cat in calibrationinstancecatalogs.values():
	allcats += cat
calibrationinstancecatalogs["all"] = allcats

calibrationinstancecataloglabel = Label(calibrationframe, text="data: ")
calibrationinstancecataloglabel.grid(row=0, column=0, sticky=W)

calibrationinstancecatalog = StringVar()
calibrationinstancecatalog.set(calibrationinstancecatalogs.keys()[0])
calibrationinstancecatalogmenu = OptionMenu(calibrationframe, calibrationinstancecatalog, *calibrationinstancecatalogs)
calibrationinstancecatalogmenu.grid(row=0, column=1, sticky=W)

numsetslabel = Label(calibrationframe, text="Number of sets: ")
numsetslabel.grid(row=1, column=0, sticky=W)

numsets = IntVar()
numsets.set(1)
numsetsspinner = Spinbox(calibrationframe, from_=1, to=10, textvariable=numsets)
numsetsspinner.grid(row=1, column=1, sticky=W)

calibrationbutton = Button(calibrationframe, text="Run Calibration", command=runcalibration)
calibrationbutton.grid(row=2, column=0)

notebook.add(calibrationframe, text="Calibration")

## Validation
validationframe = Frame(notebook)

vscriptlabel = Label(validationframe, text="Task:")
vscriptlabel.grid(row=0, column=0, sticky=W)

vscriptfiles = sorted(os.listdir(validationscriptsdir))
vscripts = OrderedDict()
for f in vscriptfiles:
	if not f.startswith("."):
		vscriptname = f.replace(".sh", "")
		vscriptname = ' & '.join(a+b for a,b in zip(vscriptname[::2], vscriptname[1::2]))
		vscripts[vscriptname] = os.path.join(validationscriptsdir, f)
vscript = StringVar()
vscript.set(vscripts.keys()[0])
vscriptmenu = OptionMenu(validationframe, vscript, *vscripts)
vscriptmenu.grid(row=0, column=1, sticky=W)

validationbutton = Button(validationframe, text="Run Validation", command=runvalidation)
validationbutton.grid(row=1, column=0)

notebook.add(validationframe, text="Validation")

notebook.grid(row=1, column=0)

mainloop()

browserapp.stop()
