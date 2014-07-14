#!/usr/bin/python
# made from tutorial http://wiki.wxpython.org/AnotherTutorial

"""
Copyright 2013 Ian Campbell, Kevin "Feng" Chen, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import wx
import subprocess
import functools
import shutil
import paramiko
import stat
import errno
import itertools
import copy
import  wx.lib.newevent
import threading
from fnmatch import fnmatch
from sys import argv

basefolder = None
basefoldername = None
basefoldermakedirs = False
if len(argv) > 1:
	basefolder = argv[1]
	if len(argv) > 2:
		basefoldername = argv[2]
		if len(argv) > 3:
			if argv[3].lower() == "true":
				basefoldermakedirs = True

sshconnections = []
clipboard = []
cut = False

class PopupDialog(wx.Dialog):
    """ A popup dialog for temporary user messages """
	
    def __init__(self, parent, title, msg):
        # Create a dialog
        wx.Dialog.__init__(self, parent, -1, title, size=(350, 150), style=wx.CAPTION | wx.STAY_ON_TOP)
        # Add sizers
        box = wx.BoxSizer(wx.VERTICAL)
        box2 = wx.BoxSizer(wx.HORIZONTAL)
        # Add an Info graphic
        bitmap = wx.EmptyBitmap(32, 32)
        bitmap = wx.ArtProvider_GetBitmap(wx.ART_INFORMATION, wx.ART_MESSAGE_BOX, (32, 32))
        graphic = wx.StaticBitmap(self, -1, bitmap)
        box2.Add(graphic, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, 10)
        # Add the message
        message = wx.StaticText(self, -1, msg)
        box2.Add(message, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 10)
        box.Add(box2, 0, wx.EXPAND)
        # Handle layout
        self.SetAutoLayout(True)
        self.SetSizer(box)
        self.Fit()
        self.Layout()
        self.CentreOnScreen()
        # Display the Dialog
        self.Show()
        # Make sure the screen gets fully drawn before continuing.
        wx.Yield()

class LocalPath(object):
	def __init__(self, path):
		self.path = path
		self.remote = False
		self.updatepath()
		self.parents = None
	
	def name(self):
		return os.path.split(os.path.abspath(self.path))[1]
	
	def listdir(self):
		files = os.listdir(self.path)
		return [self.join(f) for f in files]
		
	def join(self, relpath):
		return LocalPath(os.path.join(self.path, relpath))
		
	def isfile(self):
		return os.path.isfile(self.path)
		
	def isdir(self):
		return os.path.isdir(self.path)
	
	def dirname(self):
		return os.path.dirname(self.path)
		
	def basename(self):
		return os.path.basename(self.path)
		
	def mkdir(self):
		return os.mkdir(self.path)
		
	def remove(self):
		if self.file:
			return os.remove(self.path)
		elif self.dir:
			return shutil.rmtree(self.path)
			
	def copy(self, newpath):
		if newpath.remote:
			newpath.copyfrom(self)
		else:
#			newpath = newpath.join(self.basename())
			if self.file:
				shutil.copy(self.path, newpath.path)
			elif self.dir:
				shutil.copytree(self.path, newpath.path)
			
	def move(self, newpath):
		self.copy(newpath)
		self.remove()
		
	def invisible(self):
		return self.basename()[0] == "."
		
	def generatesimilarpath(self, path):
		return LocalPath(path)
		
	def getparents(self):
		if self.parents == None:
			parents = []
			currentparent = self.path
			while True:
				currentparent = os.path.dirname(currentparent)
				if currentparent in parents:
					break
				parents.append(currentparent)
			
			self.parents = [self.generatesimilarpath(p) for p in parents]
			
		return self.parents
		
	def ischild(self, parent, include=False):
		return (parent in self.getparents() or (parent == self and include))
		
	def __eq__(self, obj):
		if type(obj) == LocalPath:
			if self.normpath == obj.normpath:
				return True
		return False
		
	def __repr__(self):
		return "LocalPath<" + self.path + ">"
		
	def exists(self):
		return os.path.exists(self.path)
		
	def create(self):
		f = open(self.path,'w')
		f.close()
	
	def canextract(self, extractors):
		for ending, program in extractors:
			if self.path.endswith(ending):
				return True
		return False
		
	def cancompress(self):
		return True
	
	def candownload(self):
		return self.file

	def canimport(self):
		return True
		
	def canedit(self):
		return os.access(self.path, os.W_OK)
		
	def updatepath(self):
		self.updatenormpath()
		self.updateisdir()
		
	def updatenormpath(self):
		self.normpath = os.path.normpath(self.path)
		
	def updateisdir(self):
		self.dir = self.isdir()
		self.file = self.isfile()
		
	def makedirs(self):
		os.makedirs(self.path)
		
	def getsize(self):
		num = os.path.getsize(self.path)
		for x in ['bytes','KB','MB','GB','TB']:
			if num < 1024.0:
				return "%3.1f %s" % (num, x)
			num /= 1024.0

class RemotePath(LocalPath):
	def __init__(self, sftp, path, servername):
		self.sftp = sftp
		self.servername = servername
		LocalPath.__init__(self, path)
		self.remote = True
	
	def listdir(self):
		files = self.sftp.listdir(self.path)
		return [self.join(f) for f in files]
		
	def join(self, path):
		return RemotePath(self.sftp, os.path.join(self.path, path), self.servername)
		
	def __eq__(self, obj):
		if type(obj) == RemotePath:
			if self.normpath == obj.normpath:
				return True
		return False
		
	def __repr__(self):
		return "RemotePath<" + self.servername + ":" + self.path + ">"
		
	def isfile(self):
		return stat.S_ISREG(self.sftp.stat(self.path).st_mode)
		
	def isdir(self):
		return stat.S_ISDIR(self.sftp.stat(self.path).st_mode)
		
	def generatesimilarpath(self, path):
		return RemotePath(self.sftp, path, self.servername)
		
	def copy(self, newpath):
		if newpath.remote:
			pass
#			fo = self.sftp.open(self.path)
#			newpath.sftp.putfo(fo, newpath.path)
#			fo.close()
		else:
#			newpath = newpath.join(self.basename())
			if self.dir:
				newpath.mkdir()
				for content in self.listdir():
					content.copy(newpath)
			elif self.file:
				self.sftp.get(self.path, newpath.path)

	def copyfrom(self, path):
#		self.sftp.put(self.path, path.path)
		pass
		
	def exists(self):
		try:
			self.sftp.stat(self.path)
		except IOError, e:
			if e.errno == errno.ENOENT:
				return False
			raise
		else:
			return True
			
	def create(self):
		fo = self.sftp.open(self.path, "w")
		fo.close()

	def remove(self):
		if self.file:
			return self.sftp.remove(self.path)
		elif self.dir:
			for content in self.listdir():
				content.remove()
			return self.sftp.rmdir(self.path)
	
	def canextract(self, extractors):
		return False
		
	def cancompress(self):
		return False
	
	def candownload(self):
		return False

	def canimport(self):
		return False
	
	def mkdir(self):
		self.sftp.mkdir(self.path)
		
	def canedit(self):
		return False
		
	def getsize(self):
		return "Unknown"

class FileTree(wx.TreeCtrl):
	def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
				 size=wx.DefaultSize, style=wx.TR_DEFAULT_STYLE| wx.TR_HIDE_ROOT,
				 validator=wx.DefaultValidator, name="",
				 file_filter=("*.*")):
		wx.TreeCtrl.__init__( self, parent, id, pos, size, style, validator, name)
		self.file_filter = file_filter
		il = wx.ImageList(16,16)
		self.fldridx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_FOLDER,
													   wx.ART_OTHER, (16,16)))
		self.fldropenidx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_FOLDER,
														  wx.ART_OTHER, (16,16)))
		self.fileidx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE,
													   wx.ART_OTHER, (16,16)))
		self.rootidx = il.Add(wx.ArtProvider.GetBitmap(wx.ART_HARDDISK,
														wx.ART_OTHER, (16, 16)))

		self.allpaths = {}
		self.roots = {}
		
		self.AssignImageList(il)
		self.root = self.AddRoot("Root")#os.path.split(rootfolder)[1])
		self.SetItemImage(self.root,self.fldridx,wx.TreeItemIcon_Normal)
		self.SetItemImage(self.root, self.fldropenidx,wx.TreeItemIcon_Expanded)
		
#		self.AddTreeNodes(root, rootfolder)
		
		try:
			self.Expand(self.root)
		except:
			# not if we have a hidden root
			pass
			
	def addRootFolder(self, folder, expand=True, select=True, name=None):
		if name == None:
			name = folder.path
			if folder.remote:
				name += " on " + folder.servername
		item = self.addItem(self.root, folder, name, rootitem=True)
		self.roots[folder] = item
		if select:
			self.SelectItem(item)
		if expand:
			self.Expand(item)
			
	def expanding(self, event):
		item = event.GetItem()
		self.AddTreeNodes(item, self.GetPyData(item))
	
	def collapsing(self, event):
		item = event.GetItem()
		folder = self.GetPyData(item)
		pathstodelete = []
		for path in self.allpaths.keys():
			if path.ischild(folder):
				pathstodelete.append(path)
		for path in pathstodelete:
			self.rempath(path)
		self.DeleteChildren(item)
		
	def rempath(self, path):
		return self.allpaths.pop(path)
		
	def delitem(self, item):
		self.Delete(item)
		
	def addItem(self, parentItem, itempath, name=None, rootitem=False):
		if name == None:
			name = itempath.name()#os.path.split(os.path.abspath(itempath))[1]
		newItem = self.AppendItem(parentItem, name)#folder)
		self.SetItemHasChildren(newItem, True)
		if rootitem:
			self.SetItemImage(newItem, self.rootidx,wx.TreeItemIcon_Normal)
			self.SetItemImage(newItem, self.rootidx, wx.TreeItemIcon_Expanded)
		else:
			self.SetItemImage(newItem, self.fldridx,wx.TreeItemIcon_Normal)
			self.SetItemImage(newItem, self.fldropenidx, wx.TreeItemIcon_Expanded)
#		self.AddTreeNodes(newItem, itempath)
		self.SetPyData(newItem, itempath)
		self.Bind(wx.EVT_TREE_ITEM_EXPANDING, self.expanding)
		self.Bind(wx.EVT_TREE_ITEM_COLLAPSING, self.collapsing)
		self.allpaths[itempath] = newItem
		
		return newItem
	
	def getchildnodes(self, rootfolder):
		items = rootfolder.listdir()#os.listdir(rootfolder)
		items = sorted(items,key=lambda f : f.basename().lower())
		folders = []
		for item in items:
#			if item[0]==".":
#				continue
#			itempath = item#rootfolder.join(item)#os.path.join(rootfolder, item)
			if item.dir:#os.path.isfile(itempath):
				folders.append(item)
		return folders
	
	def AddTreeNodes(self, parentItem, rootfolder):
		folders = self.getchildnodes(rootfolder)
		for item in folders:
			self.addItem(parentItem, item)
		
	def GetPath(self):
		return self.GetPyData(self.GetSelection())
		
	def getitem(self, path):
		for item in self.allpaths:
			if path == item:
				return self.allpaths[item]
		print "couldnt find item", path.path
		
	def SetPath(self, path):
		pathroot = None
		for root in self.roots:
			if path.ischild(root, include=True):
				pathroot = root
		if pathroot == None:
			return
		parents = path.getparents()
		parents.reverse()
		if pathroot in parents:
			self.Expand(self.getitem(pathroot))
			
			rootindex = parents.index(pathroot)
			parents = parents[rootindex+1:]
			
			for sub in parents:
				self.Expand(self.getitem(sub))
		
		self.SelectItem(self.getitem(path))
		
	def refresh(self):
		toremove = []
		for path in self.allpaths.keys():
			item = self.allpaths[path]
			if not path.exists():
				item = self.rempath(path)
				toremove.append(item)
			else:
				if self.IsExpanded(item):
					children = []
					for newpath in self.allpaths.keys():
						if newpath.ischild(path):
							children.append(newpath)
					newchildren = self.getchildnodes(path)
					for newchild in newchildren:
						if not newchild in children:
							self.addItem(item, newchild)
		for item in toremove:
			self.delitem(item)

class MyFrame(wx.Frame):
	showingsubmit = False
	lochistory = []
	
	extractors = [(".tar.bz2", "tar xjvf"), 
		(".tar.gz", "tar xzvf"), 
		(".bz2", "bunzip2"), 
		(".rar", "rar x"), 
		(".gz", "gunzip"), 
		(".tar", "tar xf"), 
		(".tbz2", "tar xjvf"), 
		(".tgz", "tar xzvf"), 
		(".zip", "unzip"), 
		(".Z", "uncompress"), 
		(".7z", "7z x")]
	
	def __init__(self, parent, id, title, basefolder=None, basefoldername=None, basefoldermakedirs=False, openmode=False, openmethod=None, showotherfolders=True, copyinto=False, foldersonly=False):
		wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(800, 500))
		
		# menu bar
		menubar = wx.MenuBar()
		fileMenu = wx.Menu()
		importMenu = wx.Menu()
		
		newwinitem = fileMenu.Append(-1, 'New Window')
		self.Bind(wx.EVT_MENU, self.newwindow, newwinitem)
		
		refreshitem = fileMenu.Append(-1, 'Refresh')
		self.Bind(wx.EVT_MENU, self.refreshfiles, refreshitem)
		
		sshitem = importMenu.Append(-1, 'Connect to SSH Server')
		self.Bind(wx.EVT_MENU, self.sshconnect, sshitem)
		
		importitem = importMenu.Append(-1, 'Import from Local Computer')
		self.Bind(wx.EVT_MENU, self.importfile, importitem)
		
#		importitem = fileMenu.Append(-1, 'Import Files')
#		self.Bind(wx.EVT_MENU, self.importfile, importitem)
		
#		exportitem = fileMenu.Append(-1, 'Download Files')
#		self.Bind(wx.EVT_MENU, self.exportfile, exportitem)
		
#		submititem = toolsMenu.Append(-1, 'Submit')
#		self.Bind(wx.EVT_MENU, self.showsubmitwindow, submititem)
		
		menubar.Append(fileMenu, '&File')
		menubar.Append(importMenu, '&Import')
		self.SetMenuBar(menubar)
		
		# main splitters
#		self.mainsplitter = wx.SplitterWindow(self, -1, style=wx.SP_3D)
#		self.mainsplitter.SetSashGravity(1)
		self.filesplitter = wx.SplitterWindow(self, -1, style=wx.SP_3D)
		
		# file browser
#		self.dir = wx.GenericDirCtrl(self.filesplitter, -1, dir=os.environ['HOME'], style=wx.DIRCTRL_DIR_ONLY)
		
		tb = self.CreateToolBar( wx.TB_HORIZONTAL
                                 | wx.NO_BORDER
                                 | wx.TB_FLAT 
                                 | wx.TB_TEXT
                                 )
		backtoolbaritem = tb.AddSimpleTool(-1, wx.ArtProvider.GetBitmap(wx.ART_GO_BACK))
		uptoolbaritem = tb.AddSimpleTool(-1, wx.ArtProvider.GetBitmap(wx.ART_GO_UP))
		refreshtoolbaritem = tb.AddSimpleTool(-1, wx.ArtProvider.GetBitmap("gtk-refresh"))
		tb.Realize()
		
		self.Bind(wx.EVT_TOOL, self.goback, backtoolbaritem)
		self.Bind(wx.EVT_TOOL, self.goup, uptoolbaritem)
		self.Bind(wx.EVT_TOOL, self.refreshfiles, refreshtoolbaritem)
		
		dirpanel = wx.Panel(self.filesplitter)
		dirsizer = wx.BoxSizer(wx.VERTICAL)
		dirpanel.SetSizer(dirsizer)
		
		self.foldersonly = foldersonly
		
		self.dir = FileTree(dirpanel)#, rootfolder = os.environ['HOME'])
		dirsizer.Add(self.dir, 1, wx.EXPAND)
		
		lc1panel = wx.Panel(self.filesplitter)
		lc1sizer = wx.BoxSizer(wx.VERTICAL)
		lc1panel.SetSizer(lc1sizer)
		
		self.lc1 = wx.ListCtrl(lc1panel, -1, style=wx.LC_LIST)
		lc1sizer.Add(self.lc1, 1, wx.EXPAND)
		self.lc1.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.clickitem)
		self.lc1.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.rightclick)
		
		if openmode:
			openbutton = wx.Button(lc1panel, -1, "Open")
			openbutton.Bind(wx.EVT_BUTTON, self.openfiles)
			lc1sizer.Add(openbutton, 0, wx.ALIGN_RIGHT)
			self.openmethod = openmethod
		
		hasbasefolder = not (basefolder == None) and not (basefolder == [])
		select = True
		self.copyintofolder = None
		if hasbasefolder:
			if not isinstance(basefolder, list):
				basefolder = [basefolder]
				basefoldername = [basefoldername]
			for index, folder in enumerate(basefolder):
				name = basefoldername[index]
				folder = LocalPath(folder)
				if basefoldermakedirs:
					if not folder.exists():
						folder.makedirs()
				if folder.exists():
					self.dir.addRootFolder(folder, name=name, select=select)
					select = False
				if copyinto:
					if self.copyintofolder == None:
						self.copyintofolder = folder
		if showotherfolders:
			self.dir.addRootFolder(LocalPath(os.environ['HOME']), expand=not hasbasefolder, select=not hasbasefolder)
			sdatadir = LocalPath("/home/sdata/" + os.environ['USER'])
			if sdatadir.exists():
				self.dir.addRootFolder(sdatadir, select=False, expand=False)
		
		imgsize = 16
		imglist = wx.ImageList(imgsize,imgsize)
		imglist.Add(wx.ArtProvider.GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, (imgsize,imgsize)))
		imglist.Add(wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, (imgsize,imgsize)))
		self.lc1.AssignImageList(imglist, wx.IMAGE_LIST_SMALL)
		
		tree = self.dir#.GetTreeCtrl()
		
		self.filesplitter.SplitVertically(dirpanel, lc1panel, sashPosition=300)
		
		# command panel
		self.Layout()
		
		wx.EVT_TREE_SEL_CHANGED(self, tree.GetId(), self.treeselectionchanged)
		self.dir.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.rightclick)
		self.refreshfiles()
		
		self.Centre() 
		
#		self.togglesubmit()
#		self.togglesubmit()
	
	def openfiles(self, event=None):
		itemsselected = GetSelectedItems(self.lc1)
		paths = [self.getitempath(item) for item in itemsselected]
		
		if self.foldersonly and paths == []:
			paths = [self.getpath()]
		
		if self.copyintofolder != None:
			if any([not path.ischild(self.copyintofolder) for path in paths]):
				dlg = wx.MessageDialog(self, "Do you want to copy these files into the workspace directory?", "Copy files?", wx.YES_NO | wx.ICON_QUESTION)
				result = dlg.ShowModal() == wx.ID_YES
				if result:
					global clipboard, cut
					cut = False
					clipboard = [path for path in paths if not path.ischild(self.copyintofolder)]
					newpaths = self.paste(where=self.copyintofolder)
					
					paths = newpaths + [path for path in paths if path.ischild(self.copyintofolder)]
				else:
					return
		
		if any([path.remote for path in paths]):
			dlg = wx.MessageDialog(self, "Items are not in local folder.\nPlease move these items.", "Error", wx.OK | wx.ICON_WARNING)
			dlg.ShowModal()
			dlg.Destroy()
		else:
			paths = [path.path for path in paths]
			if paths == []:
				dlg = wx.MessageDialog(self, "No items selected", "Error", wx.OK | wx.ICON_WARNING)
				dlg.ShowModal()
				dlg.Destroy()
			else:
				if self.openmethod != None:
					self.openmethod(paths)
					self.Destroy()
		
	def treeselectionchanged(self, event=None):
		self.lochistory.append(self.dir.GetPath())
		self.refreshfilelist()
		
	def sshconnect(self, event):
		class sshwizard(wx.Frame):
			def __init__(self, parent, sshmethod):
				wx.Frame.__init__(self, parent, -1, 'SSH Connection Wizard')
				self.panel = wx.Panel(self, -1)
				self.sshmethod = sshmethod
				
				hbox = wx.BoxSizer(wx.VERTICAL)
				
				hostnamebox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.hostnamelabel = wx.StaticText(self.panel, label="Hostname:")
				hostnamebox.Add(self.hostnamelabel)
				
				self.hostnameentry = wx.TextCtrl(self.panel, size=(300, -1))
				self.hostnameentry.Bind(wx.EVT_KEY_DOWN, self.otherentry)
				self.hostnameentry.SetFocus()
				hostnamebox.Add(self.hostnameentry, 1)
				
				hbox.Add(hostnamebox, 0, wx.EXPAND)

				locbox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.loclabel = wx.StaticText(self.panel, label="Location:")
				locbox.Add(self.loclabel)
				
				self.locentry = wx.TextCtrl(self.panel, size=(300, -1))
				self.locentry.Bind(wx.EVT_KEY_DOWN, self.otherentry)
				locbox.Add(self.locentry, 1)
				
				hbox.Add(locbox, 0, wx.EXPAND)
				
				usernamebox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.usernamelabel = wx.StaticText(self.panel, label="Username:")
				usernamebox.Add(self.usernamelabel)
				
				self.usernameentry = wx.TextCtrl(self.panel, size=(300, -1))
				self.usernameentry.Bind(wx.EVT_KEY_DOWN, self.otherentry)
				usernamebox.Add(self.usernameentry, 1)
				
				hbox.Add(usernamebox, 0, wx.EXPAND)
				
				passwordbox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.passwordlabel = wx.StaticText(self.panel, label="Password:")
				passwordbox.Add(self.passwordlabel)
				
				self.passwordentry = wx.TextCtrl(self.panel, size=(300, -1), style=wx.TE_PASSWORD)
				self.passwordentry.Bind(wx.EVT_KEY_DOWN, self.pwentry)
				passwordbox.Add(self.passwordentry, 1)
				
				hbox.Add(passwordbox, 0, wx.EXPAND)
				
				connectbutton = wx.Button(self.panel, label="Connect")
				connectbutton.Bind(wx.EVT_BUTTON, self.connect)
				hbox.Add(connectbutton)
				
				self.panel.SetSizer(hbox)
				
				self.Centre()
				self.Show()
				
			def otherentry(self, event):
				if event.GetKeyCode() == wx.WXK_RETURN:
					win = event.GetEventObject()
					win.Navigate()
				else:
					event.Skip()
				
			def pwentry(self, event):
				if event.GetKeyCode() == wx.WXK_RETURN:
					self.connect()
				else:
					event.Skip()
				
			def connect(self, event=None):
				hostname = self.hostnameentry.GetValue()
				location = self.locentry.GetValue()
				if location == "":
					location = "/"
				username = self.usernameentry.GetValue()
				password = self.passwordentry.GetValue()
				self.sshmethod(hostname, username, password, location)
				self.Destroy()
		
		sshwizard(self, self.ssh)
	
	def showsubmitwindow(self, event=None):
		class submitwizard(wx.Frame):
			def __init__(self, parent, connectmethod):
				self.connectmethod = connectmethod
				
				wx.Frame.__init__(self, parent, -1, 'Submit Wizard')
				self.panel = wx.Panel(self, -1)
				
				hbox = wx.BoxSizer(wx.VERTICAL)
				
				self.procbox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.proclabel = wx.StaticText(self.panel, label="Number Processors: ")
				self.procbox.Add(self.proclabel)
				
				self.procentry = wx.TextCtrl(self.panel)
				self.procbox.Add(self.procentry, 1)
				
				hbox.Add(self.procbox, 0)
				
				
				self.timebox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.proclabel = wx.StaticText(self.panel, label="Time (minutes): ") 	
				self.timebox.Add(self.proclabel)
				
				self.timeentry = wx.TextCtrl(self.panel)
				self.timebox.Add(self.timeentry, 1)
				
				hbox.Add(self.timebox, 0)
				
				
				self.commandbox = wx.BoxSizer(wx.HORIZONTAL)
				
				self.commandlabel = wx.StaticText(self.panel, label="Command: ")
				self.commandbox.Add(self.commandlabel)
				
				self.commandentry = wx.TextCtrl(self.panel, size=(-1, -1))
				self.commandbox.Add(self.commandentry, 1)
				
				hbox.Add(self.commandbox, 0, wx.EXPAND)
				
				
				self.runbutton = wx.Button(self.panel, label="Run")
				self.runbutton.Bind(wx.EVT_BUTTON, self.run)
				hbox.Add(self.runbutton)
				
				self.panel.SetSizer(hbox)
				
				hbox.Layout()
				
				self.Centre()
				self.Show()
				
			def run(self, event=None):
				proc = self.procentry.GetValue()
				time = self.timeentry.GetValue()
				command = self.commandentry.GetValue()
				self.connectmethod(proc, time, command)
			
		wiz = submitwizard(self, self.runsubmit)
	
	def ssh(self, hostname, username, password, location="/"):
		ssh = paramiko.SSHClient()
		ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(hostname, username=username, password=password)
		sshconnections.append(ssh)
		sftp = ssh.open_sftp()
		sshconnections.append(sftp)
		
		path = RemotePath(sftp, location, username + "@" + hostname)
		if not path.exists():
			path = RemotePath(sftp, "/", username + "@" + hostname)
		
		self.dir.addRootFolder(path)
	
	def importfile(self, event=None, where=None):
		if where == None:
			where = self.getpath()
		subprocess.Popen(["importfile"], cwd=where.path)
	
	def exportfile(self, event=None, items=None):
		if items == None:
			items = GetSelectedItems(self.lc1)
			items = [self.getitempath(item) for item in items]
		files = items
		print files
		for thefile in files:
			subprocess.Popen(["exportfile", thefile.path], cwd=self.dir.GetPath().path)
	
	def runsubmit(self, proc, time, comm):
		command = ["submit", "-n", str(proc), "-w", str(time)] + str(comm).split(" ")
		print command
		ret = wx.MessageBox("Run the command:\n" + str(command) + "\nin the directory:\n" + self.dir.GetPath().path, 'Command', wx.YES_NO | wx.CENTRE |
		wx.NO_DEFAULT, self)
		if ret == wx.YES:
			print "running command"
			process = subprocess.Popen(command, cwd=self.dir.GetPath().path)
			process.wait()
			print "command ended"
			self.refreshfiles()
	
	def refreshfiles(self, event=None):
		self.dir.refresh()
		self.refreshfilelist()
		
	def refreshfilelist(self, event=None):
		path = self.dir.GetPath()
		if path != None:
			list = path.listdir()#os.listdir(path)
			self.lc1.ClearAll()
			folders = []
			files = []
			for path in list:
#				if not path.invisible():
#					path = name#self.getpath(name)
				if path.dir:#os.path.isdir(path):
					folders.append(path)
				else:
					if not self.foldersonly:
						files.append(path)
			folders.sort(reverse=True, key=lambda x : x.name().lower())
			files.sort(reverse=True, key=lambda x : x.name().lower())
			for path in files:
				self.lc1.InsertImageStringItem(0, path.name(), 1)
			for path in folders:
				self.lc1.InsertImageStringItem(0, path.name(), 0)
	
	def getitempath(self, item):
		return self.getpath(item.GetText())
		
	def getpath(self, name=None):
		if name == None:
			return self.dir.GetPath()
		else:
			return self.dir.GetPath().join(name)
#		return self.dir.GetPath() + "/" + name
	
	def clickitem(self, event=None, item=None):
		if item == None:
			item = event.GetItem()
		path = self.getitempath(item)
		if path.dir:#os.path.isdir(path):
			self.dir.SetPath(path)
#			self.refreshfiles()
#		else:
#			self.openitem(item)
	
	def openitem(self, item=None):
		subprocess.Popen(["xdg-open", item.path], cwd=self.getpath().path)
	
	def openitems(self, event=None, items=None):
		if items == None:
			items = [self.getitempath(event.GetItem())]
		for item in items:
			self.openitem(item)
		
	def newfolder(self, event=None, where=None):
		if where == None:
			where = self.getpath()
		test = wx.TextEntryDialog(None, "New folder name:", "New Folder")
		if test.ShowModal() == wx.ID_OK:
			foldername = test.GetValue()
			folderpath = where.join(foldername)
			folderpath.mkdir()#os.mkdir(folderpath)
			self.refreshfiles()
	
	def newfile(self, event=None, where=None):
		if where == None:
			where = self.getpath()
		test = wx.TextEntryDialog(None, "New file name:", "New File")
		if test.ShowModal() == wx.ID_OK:
			filename = test.GetValue()
			filepath = where.join(filename)
			filepath.create()
			self.refreshfiles()
			
	def deleteitems(self, event=None, items=None):
		ret = wx.MessageBox("Delete these items?", 'Delete', wx.YES_NO | wx.CENTRE |
		wx.NO_DEFAULT, self)
		if ret == wx.YES:
			for item in items:
				self.deletepath(item)
			self.refreshfiles()
		
	def deletepath(self, path):
		path.remove()
#		if path.isdir():#os.path.isdir(path):
#			shutil.rmtree(path)
#		else:
#			path.remove()#os.remove(path)
			
	def cutitems(self, event=None, items=None):
		global clipboard, cut
		cut = True
		clipboard = items
		
	def copyitems(self, event=None, items=None):
		global clipboard, cut
		cut = False
		clipboard = items
		
	def paste(self, event=None, where=None):
		if where == None:
			where = self.getpath()
		global clipboard, cut
		print "pasting", clipboard, "into", where
		loadDlg = PopupDialog(self, ("Working..."), ("Working.\nPlease wait...."))
		newdir = where
		newfiles = []
		for path in clipboard:
			if path.dirname() == newdir:
				basename = path.basename().split(".", 1)
				newdir = where.join(basename[0] + " (copy)." + basename[1])
			else:
				newdir = where.join(path.basename())
			path.copy(newdir)
			newfiles.append(newdir.join(path.basename()))
#			if path.isdir():#os.path.isdir(path):
#				shutil.copytree(path, newdir)
#			else:
#				shutil.copy(path, newdir)
			if cut:
				self.deletepath(path)
		self.refreshfiles()
		loadDlg.Destroy()
		return newfiles
		
	def renamefile(self, event=None, item=None):
		path = item
		name = path.basename()#os.path.basename(path)
		test = wx.TextEntryDialog(None, "New file name:", "New File", name)
		if test.ShowModal() == wx.ID_OK:
			newname = test.GetValue()
			newpath = self.getpath(newname)
			path.move(newpath)
#			shutil.copy(path, newpath)
#			self.deletepath(path)
			self.refreshfiles()
			
	def compress(self, event=None, items=None):
		test = wx.TextEntryDialog(None, "Archive file name:", "Compress Files")
		if test.ShowModal() == wx.ID_OK:
			newname = test.GetValue()
			process = subprocess.Popen(["zip", newname] + [item.basename() for item in items], cwd=self.getpath().path)
			process.wait()
			self.refreshfiles()
		
	def extract(self, event=None, item=None):
		itempath = item.path
		command = None
		for ending, prg in self.extractors:
			if itempath.endswith(ending):
				command = prg
		process = subprocess.Popen(command.split(" ") + [itempath], cwd=self.getpath().path)
		process.wait()
		self.refreshfiles()
	
#	def togglesubmit(self, event=None):
#		self.showingsubmit = not self.showingsubmit
#		if self.showingsubmit:
#			self.mainsplitter.SplitHorizontally(self.filesplitter, self.panel)
#		else:
#			self.mainsplitter.Unsplit()
	
	def fileinfo(self, event=None, item=None):
		size = item.getsize()
		dlg = wx.MessageDialog(self, "Path: " + item.path + "\nSize: "+size, "File Info", wx.OK)
		dlg.ShowModal()
		dlg.Destroy()
		
	def rightclick(self, event):
		eventsource = event.GetEventObject()
		if eventsource == self.lc1:
			hasitem = event.GetIndex() != -1
			itemclicked = event.GetItem()
			if hasitem:
				itemclickedpath = self.getitempath(itemclicked)
			else:
				itemclickedpath = self.getpath()
			itemsselected = [self.getitempath(item) for item in GetSelectedItems(self.lc1)]
		elif eventsource == self.dir:
			hasitem = True
			itemclicked = event.GetItem()
			itemclickedpath = self.getpath()
			itemsselected = [itemclickedpath]
		
		canedit = all([item.canedit() for item in itemsselected])
		
		items = []
		if itemclickedpath.dir:
			if canedit:
				items += [
					("New folder", functools.partial(self.newfolder, where=itemclickedpath)),
					("New file", functools.partial(self.newfile, where=itemclickedpath)),
					None,
					("Paste", functools.partial(self.paste, where=itemclickedpath))]
				if itemclickedpath:
					items += [None, ("Import File", functools.partial(self.importfile, where=itemclickedpath))]
				if hasitem:
					items += [None]
		if hasitem:
			if canedit:
				items += [
					("View File", functools.partial(self.openitems, items=[itemclickedpath]))]
				items += [
					None,
					("Cut", functools.partial(self.cutitems, items=itemsselected))]
			items += [
				("Copy", functools.partial(self.copyitems, items=itemsselected))]
			if all([item.candownload() for item in itemsselected]):
				items += [
					None, 
					("Download", functools.partial(self.exportfile, items=itemsselected))]
			if canedit:
				items += [
					None, 
					("Rename", functools.partial(self.renamefile, item=itemclickedpath)),
					None]
				if all([item.cancompress() for item in itemsselected]):
					items += [
						("Compress", functools.partial(self.compress, items=itemsselected))]
				if itemclickedpath.canextract(self.extractors):
					items += [
						("Extract Here", functools.partial(self.extract, item=itemclickedpath))]
				items += [None, 
					("Delete", functools.partial(self.deleteitems, items=itemsselected))]
			if itemclickedpath.file:
				items += [None, 
					("File Information", functools.partial(self.fileinfo, item=itemclickedpath))]
		
		if items != []:
			menu = wx.Menu()
			for item in items:
				if item == None:
					menu.AppendSeparator()
				else:
					label, funct = item
					menuitem = menu.Append(-1, label)
					self.Bind(wx.EVT_MENU, funct, menuitem)
			
			eventsource.PopupMenu(menu, event.GetPoint())
			menu.Destroy()
	
	def newwindow(self, event=None):
		launchWindow()
		
	def goup(self, event=None):
		self.dir.SetPath(self.dir.GetPath().getparents()[0])
		
	def goback(self, event=None):
		if len(self.lochistory) > 1:
			self.lochistory.pop()
			self.dir.SetPath(self.lochistory.pop())

def GetSelectedItems(listCtrl):
	selection = []
	index = listCtrl.GetFirstSelected()
	if index != -1:
		selection.append(index)
	while True:
		index = listCtrl.GetNextSelected(index)
		if index == -1:
			break
		selection.append(index)
	
	selection = [listCtrl.GetItem(index) for index in selection]
	
	return selection

class browserapp():
	def __init__(self):
		self.runningmainloop = False
		self.newwindowevent, self.EVT_NEW_WINDOW = wx.lib.newevent.NewEvent()
		self.thread = None
	
	def launchWindow(self, **kargs):
		frame = MyFrame(None, -1, "File Browser", **kargs)
		frame.Show(True)
		if not self.runningmainloop:
			self.mainloop()
			
	def threadwindow(self, **kargs):
		launchwindowmethod = functools.partial(self.launchWindow, **kargs)
		if self.runningmainloop:
			wx.CallAfter(launchwindowmethod)
		else:
			self.thread = threading.Thread(target=launchwindowmethod)
			self.thread.start()
			
	def stop(self):
		app.Exit()
		if self.thread != None:
			self.thread.join()

	def mainloop(self):
		global sshconnections, app
		self.runningmainloop=True
		try:
			app.MainLoop()
		finally:
			for connection in sshconnections:
				connection.close()
		self.runningmainloop=False

#class MyApp(wx.App):
#	def OnInit(self):
#		frame = MyFrame(None, -1, "File Browser")
#		frame.Show(True)
#		self.SetTopWindow(frame)
#		return True

app = wx.App(False)
#app = MyApp(0)
#	app.MainLoop()
if __name__ == "__main__":
	theapp = browserapp()
	theapp.launchWindow(basefolder=basefolder, basefoldername=basefoldername, basefoldermakedirs=basefoldermakedirs)
	#theapp.launchWindow(basefolder=basefolder, basefoldername=basefoldername, basefoldermakedirs=basefoldermakedirs)
#	mainloop()
