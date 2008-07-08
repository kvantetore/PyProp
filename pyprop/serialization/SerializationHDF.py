import numpy 
import os
import time

DEBUG = False

def SaveWavefunctionHDF(hdfFile, datasetPath, psi, conf=None):
	"""
	Saves the wavefunction data to a dataset in a HDF file.

	Parameters
	--------------------------------------------------------------------
	- hdfFile can be either a string or an instance of tables.File.
	1) type(hdf) == str: 
		hdfFile is a string containing a filename to a hdffile which 
		may or may not exist. If it exists, no other process my access 
		the file while this function is called
	2) type(hdf) == tables.File:
		hdfFile is a pytables file object. This object will be closed
		and reopened during the call to this function, so all references
		to objects within this file will be invalid when this function
		exists
	
	- datasetName is the path to the dataset where the wavefunction data
	should be stored.

	- psi is the wavefunction

	- other arguments: conf, a pyprop.Config object, which contains a cofig-
	                   parser object. The configparser object is stored as
	                   an attribute on the wavefunction-array.
	
	Example:
	SaveWavefunctionHDF("myfile.h5", "/mygroup/wavefunction", prop.psi):
	"""

	if type(hdfFile) == str:
		#if the first param is a filename
		filename = str(hdfFile)
	else:
		#assume the first param is a hdfFile
		filename = hdfFile.filename
		#we must first close the hdffile. this is required for
		#multiproc scenarios anyway, so we may as well do it here
		hdfFile.close()

	distr = psi.GetRepresentation().GetDistributedModel()
	if distr.IsSingleProc():
		t = - time.time()
		RemoveExistingDataset(filename, datasetPath)
		SaveLocalSlab(filename, datasetPath, psi)
		t += time.time()
		if DEBUG: print "Duration: %.10fs" % t
		SaveConfigObject(filename, datasetPath, conf)

	else:
		#let the processors save their part one by one
		procCount = distr.ProcCount
		procId = distr.ProcId
		localSize = psi.GetData().nbytes / 1024.**2
		distr.GlobalBarrier()
		if procId==0 and DEBUG: print "Saving %s..." % filename
		for i in range(procCount):
			if procId == i:
				if procId == 0:
					RemoveExistingDataset(filename, datasetPath)

				if DEBUG: print "    Process %i writing hyperslab of %iMB" % (procId, localSize)
				t = - time.time()
				SaveLocalSlab(filename, datasetPath, psi)
				t += time.time()
				if DEBUG: print "    Duration: %.10fs" % t

			distr.GlobalBarrier();
		if procId==0 and DEBUG: print "Done."

		#proc 0 save config object
		if procId == 0:
			SaveConfigObject(filename, datasetPath, conf)

def RemoveExistingDataset(filename, datasetPath):
	if not os.path.exists(filename):
		return

	f = tables.openFile(filename, "a")
	try:
		groupName, datasetName = GetDatasetName(datasetPath)
		try:
			f.removeNode(groupName, datasetName)
		except tables.NoSuchNodeError:
			pass
			
	finally:
		f.close()

def SaveLocalSlab(filename, datasetPath, psi):
	#open file
	f = tables.openFile(filename, "a")

	#get dataset
	fullShape = tuple(psi.GetRepresentation().GetFullShape())
	dataset = GetExistingDataset(f, datasetPath)
	if dataset == None:
		#create new dataset
		dataset = CreateDataset(f, datasetPath, fullShape)
	else:
		#check that is has the correct size
		if not dataset.shape == fullShape:
			raise "Invalid shape on existing dataset. Got %s, expected %s" % (dataset.shape, fullShape)

	#get hyperslab for this proc
	fileSlab = GetFileSlab(psi)
	
	#write data
	dataset[fileSlab] = psi.GetData()

	#close file
	dataset.close()
	f.close()


def LoadWavefunctionHDF(hdfFile, datasetPath, psi):
	"""
	Loads the wavefunction from a HDF5 dataset
	"""
	if type(hdfFile) == str:
		#if the first param is a filename
		filename = str(hdfFile)
	else:
		#assume the first param is a hdfFile
		filename = hdfFile.filename

	distr = psi.GetRepresentation().GetDistributedModel()
	if distr.IsSingleProc():
		LoadLocalSlab(filename, datasetPath, psi)

	else:
		#let the processors save their part one by one
		procCount = distr.ProcCount
		procId = distr.ProcId
		for i in range(procCount):
			if procId == i:
				if DEBUG: print "Loading slab in proc ", procId
				LoadLocalSlab(filename, datasetPath, psi)
			distr.GlobalBarrier();
	
def LoadLocalSlab(filename, datasetPath, psi):
	#open file
	f = tables.openFile(filename, "r")

	#get dataset
	fullShape = tuple(psi.GetRepresentation().GetFullShape())
	if not datasetPath.startswith("/"):
		datasetPath = "/" + datasetPath
	dataset = f.getNode(datasetPath)

	#check that the dataset has the correct shape
	if not dataset.shape == fullShape:
		raise "Invalid shape on existing dataset. Got %s, expected %s" % (dataset.shape, fullShape)

	#get hyperslab for this proc
	fileSlab = GetFileSlab(psi)
	
	#load data
	#print "loading data of shape ", dataset.shape, " into wavefunction of shape ", tuple(psi.GetData().shape)
	#print "current slab = ", fileSlab
	data = dataset[fileSlab]
	#print "datashape = ", data.shape
	psi.GetData()[:] = data

	#close file
	dataset.close()
	f.close()


def GetFileSlab(psi):
	"""
	Get the local slab for the current processor in the full wavefunction
	"""
	repr = psi.GetRepresentation()
	distr = repr.GetDistributedModel()

	#get shapes
	fullShape = tuple(repr.GetFullShape())
	localShape = tuple(psi.GetData().shape)

	#get StartIndex
	localStartIndex = numpy.zeros(len(fullShape), dtype=int)
	rank = len(fullShape)
	for i in range(rank):
		localStartIndex[i] = distr.GetLocalStartIndex(int(fullShape[i]), i)
	localStartIndex = tuple(localStartIndex)
		
	#set up hyperslabs
	fileSlab = []
	for i in range(rank):
		start = localStartIndex[i]
		end = start + localShape[i]
		fileSlab += [slice(start, end)]
	fileSlab = tuple(fileSlab)

	return fileSlab


def GetDatasetName(datasetPath):
	"""
	Translates a dataset path into a groupName and a datasetName returned
	as a tuple
	"""
	#get dataset
	if not datasetPath.startswith("/"):
		datasetPath = "/"+ datasetPath
	groupName = os.path.dirname(datasetPath)
	datasetName = os.path.basename(datasetPath)
	return groupName, datasetName

def GetExistingDataset(f, datasetPath):
	"""
	Get an existing dataset at the given path. None is returned if it does 
	not exist.
	"""
	if not datasetPath.startswith("/"):
		datasetPath = "/"+ datasetPath

	#check if dataset already exists
	dataset = None
	try:
		dataset = f.getNode(datasetPath)
	except tables.NoSuchNodeError:
		pass

	return dataset
	
def CreateDataset(f, datasetPath, fullShape):
	"""
	Creates a chunked array dataset of shape fullShape at the given path. 
	"""
	groupName, datasetName = GetDatasetName(datasetPath)
	group = "/"
	if len(groupName) > 1:
		group = f.createGroup("/", groupName[1:], createparents=True)
	atom = tables.ComplexAtom(itemsize=16)
	filters = tables.Filters(complevel=0)
	dataset = f.createCArray(group, datasetName, atom, fullShape, filters=filters)

	return dataset

def SaveConfigObject(filename, datasetPath, conf):
	if conf != None:
		h5file = tables.openFile(filename, "r+")
		try:
			h5file.setNodeAttr(datasetPath, "configObject", conf.cfgObj)
		finally:
			h5file.close()

def GetConfigFromHDF5(file, datasetPath = None, confObjName = "configObject"):
	"""
	Load a configparser object stored as an attribute on a wavefunction in a 
	HDF5 file. This can be used to recreate a pyprop.Config object, or the
	original config file.
	"""
	h5file = tables.openFile(file, "r")
	try:
		if datasetPath == None:
			for node in h5file.walkNodes():
				name = node._v_name
				if name == "wavefunction":
					cfgObj = node.getAttr(confObjName)
		else:
			cfgObj = h5file.getNodeAttr(datasetPath, confObjName)

	finally:
		h5file.close()
	
	return cfgObj

