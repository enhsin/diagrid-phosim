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

workspacedir = os.path.join(os.environ['HOME'], "phosimrun")#"/home/yyon/Workspace/phosimrun"
if not os.path.exists(workspacedir):
	os.makedirs(workspacedir)
os.chdir(workspacedir)
currentdir = os.getcwd()
installdir = module_locator.module_path()#os.path.realpath(__file__))
phosimdir = os.path.dirname(installdir)
exampledir = os.path.join(phosimdir, "examples")
bindir = os.path.join(phosimdir, "bin")
datadir = os.path.join(phosimdir, "data")
seddir = os.path.join(datadir, "SEDs")
imagedir = os.path.join(datadir, "images")
phosimpy = os.path.join(bindir, "phosim.py")
outputdir = os.path.join(currentdir, "output")
workdir = os.path.join(currentdir, "work")
calibrationfolders = [os.path.join(exampledir, "biases"), os.path.join(exampledir, "darks"), os.path.join(exampledir, "flats")]
if not os.path.exists(outputdir):
	os.mkdir(outputdir)
if not os.path.exists(workdir):
	os.mkdir(workdir)

browserapp = filebrowser.browserapp()

INSTANCECATALOG = "instancecatalog"
EXTRACOMMANDS = "extracommands"
E2ADC = "e2adc"
options = {INSTANCECATALOG:"", EXTRACOMMANDS:"", E2ADC:"0"}

def openmethod(paths, optionsindex=None, label=None):
	if len(paths) != 0:
		path = paths[0]
		name = os.path.basename(paths[0])
		label.config(text=name)
		options[optionsindex] = path

def launchworkspace():
	browserapp.threadwindow(basefolder=[workspacedir, outputdir], basefoldername=["Workspace", "Output"], basefoldermakedirs=True, openmode=False)

def launchopen(optionsindex, label):
	browserapp.threadwindow(basefolder=[workspacedir, exampledir], basefoldername=["Workspace", "Examples"], basefoldermakedirs=True, openmode=True, openmethod=functools.partial(openmethod, optionsindex=optionsindex, label=label))

def updatee2adc():
	options[E2ADC] = str(e2adcvar.get())

#launchworkspace()

def runphosimpy(instancecatalog, extracommands, e2adcflag):
	
	command = ["xterm", "-hold", "-title", "phosimrun", "-e", "python", phosimpy, instancecatalog, "-c", extracommands, 
		"-e", e2adcflag, "-o", outputdir, "-w", workdir, "-b", bindir, 
		"-d", datadir, "--sed="+seddir, "--image="+imagedir]
	
	print command
	
#	command = ["ls", "-l"]
#	
	log = open("output.log", "w")
	
	result = subprocess.Popen(command, stdout=log, stderr=log)

def runuser():
	runphosimpy(instancecatalog=options[INSTANCECATALOG], extracommands=options[EXTRACOMMANDS], e2adcflag=options[E2ADC])

#runphosimpy("","",1)

root = Tk()
root.wm_title("Phosim")

workspacebutton = Button(root, text="Open Folders", command=launchworkspace)
workspacebutton.grid(row=0, column=0)

notebook = ttk.Notebook(root)

userframe = Frame(notebook)

instancecataloglabel = Label(userframe, text="Instance Catalog:")
instancecataloglabel.grid(row=0, column=0, sticky=W)

instancecatalogdisplay = Label(userframe, text="[no file selected]")
instancecatalogdisplay.grid(row=0, column=1, sticky=W)

instancecatalogbutton = Button(userframe, text="[...]", command=functools.partial(launchopen, INSTANCECATALOG, instancecatalogdisplay))
instancecatalogbutton.grid(row=0, column=2)

extracommandslabel = Label(userframe, text="Physics Override Commands:")
extracommandslabel.grid(row=1, column=0, sticky=W)

extracommandsdisplay = Label(userframe, text="[no file selected]")
extracommandsdisplay.grid(row=1, column=1, sticky=W)

extracommandsbutton = Button(userframe, text="[...]", command=functools.partial(launchopen, EXTRACOMMANDS, extracommandsdisplay))
extracommandsbutton.grid(row=1, column=2)

e2adcvar = IntVar()

e2adccheckbox = Checkbutton(userframe, text="Run e2adc", command=updatee2adc, variable=e2adcvar)
e2adccheckbox.grid(row=2, column=0, columnspan=2, sticky=W)

outputbutton = Button(userframe, text="Run Simulation", command=runuser)
outputbutton.grid(row=3, column=0)

notebook.add(userframe, text="user")

calibrationframe = Frame(notebook)

calibrationinstancecatalogs = OrderedDict()
for calibrationfolder in calibrationfolders:
	for instancecatalog in os.listdir(calibrationfolder):
		if instancecatalog.endswith("_instcat_0"):
			calibrationinstancecatalogs[instancecatalog.replace("_instcat_0", "")] = os.path.join(calibrationfolder, instancecatalog)
#calibrationinstancecatalogs["all"]

calibrationinstancecataloglabel = Label(calibrationframe, text="data: ")
calibrationinstancecataloglabel.grid(row=0, column=0)

calibrationinstancecatalog = StringVar()
calibrationinstancecatalog.set(calibrationinstancecatalogs.keys()[0])
calibrationinstancecatalogmenu = OptionMenu(calibrationframe, calibrationinstancecatalog, *calibrationinstancecatalogs)
calibrationinstancecatalogmenu.grid(row=0, column=1)



notebook.add(calibrationframe, text="calibration")

notebook.grid(row=1, column=0)

mainloop()

#time.sleep(100)
#if browserapp.thread != None:
#	browserapp.thread.join()
browserapp.stop()
