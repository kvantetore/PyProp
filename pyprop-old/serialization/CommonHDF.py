import numpy 
import os
import time

from numpy import r_, s_, all, diff, array, asarray

DEBUG = False

#--------------------------------------------------------------------------------------
#                         Dataset Tools
#--------------------------------------------------------------------------------------

def RemoveExistingDataset(filename, datasetPath):
	"""
	Removes a node from a HDF5 file if it exists
	otherwise, do nothing
	"""
	if not os.path.exists(filename):
		return

	f = tables.openFile(filename, "a")
	try:
		groupName, datasetName = GetDatasetName(datasetPath)
		try:
			f.removeNode(groupName, datasetName, recursive=True)
		except tables.NoSuchNodeError:
			pass
			
	finally:
		f.close()


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
	#Separate path from node name 
	groupName, datasetName = GetDatasetName(datasetPath)

	#Iterative over path (for each /), determine if each node
	#in the path exists, if not, create it
	pathList = groupName.split("/")
	curPath = ""
	prevPath = ""
	for curGroupName in pathList:
		#Skip zero-length node names (from the split operation)
		if len(curGroupName) == 0:
			continue

		#Create node if it does not exist
		prevPath = curPath
		curPath += "/%s" % curGroupName
		if not f.__contains__(curPath):
			if prevPath == "":
				newGroup = f.createGroup("/", curGroupName)
			else:
				newGroup = f.createGroup(prevPath, curGroupName)

	#Finally, save the data set
	atom = tables.ComplexAtom(itemsize=16)
	filters = tables.Filters(complevel=0)
	dataset = f.createCArray(groupName, datasetName, atom, fullShape, filters=filters)

	return dataset

#--------------------------------------------------------------------------------------
#                         Load and Save slabs of large datasets
#--------------------------------------------------------------------------------------

def SaveLocalSlab(filename, datasetPath, localData, localSlab, fullShape):
	"""
	Saves the local slab of a global dataset to a file

	filename    - the hdf5-file to save to
	datasetPath - the path in <filename> where the dataset should be stored
	localData   - the n-dimensional array containing the local data to be stored
	localSlab   - the local slab (n-dimensional tuple containing the local range in
				  the global dataset
	fullShape   - shape of the global dataset

	"""
	#open file
	f = tables.openFile(filename, "a")
	try:
		#get dataset
		localShape = localData.shape
		dataset = GetExistingDataset(f, datasetPath)
		if dataset == None:
			#create new dataset
			dataset = CreateDataset(f, datasetPath, fullShape)
		else:
			#check that is has the correct size
			if not dataset.shape == fullShape:
				raise "Invalid shape on existing dataset. Got %s, expected %s" % (dataset.shape, fullShape)
		
		#write data
		isSlices = map(lambda s: isinstance(s, slice), localSlab)
		if all(isSlices):
			dataset[localSlab] = localData
		elif all(isSlices[1:]):
			for dataIdx, fileIdx in enumerate(localSlab[0]):
				curLocalSlab = (fileIdx,) + localSlab[1:]
				dataset[curLocalSlab] = localData[dataIdx, :]
		else:
			raise "Only the first rank may be fancy indexing"

	
	finally:
		#Make sure file is closed
		f.close()


def LoadLocalSlab(filename, datasetPath, localData, localSlab, fullShape):
	"""
	Loads the local slab of a global dataset from a file

	filename    - the hdf5-file to load from
	datasetPath - the path in <filename> where the dataset is stored
	localData   - the n-dimensional array where to store the local data
	localSlab   - the local slab (n-dimensional tuple containing the local range in
				  the global dataset
	fullShape   - shape of the global dataset

	"""
	#open file
	f = tables.openFile(filename, "r")
	try:
		#get dataset
		localShape = localData.shape
		dataset = GetExistingDataset(f, datasetPath)
		if dataset == None:
			raise Exception("Dataset '%s' not found in file '%s'" % (datasetPath, filename))
		else:
			#check that is has the correct size
			if not dataset.shape == fullShape:
				raise Exception("Invalid shape on existing dataset. Got %s, expected %s" % (dataset.shape, fullShape))
		
		#load data
		isSlices = map(lambda s: isinstance(s, slice), localSlab)
		if all(isSlices):
			localData[:] = dataset[localSlab]
		elif all(isSlices[1:]):
			for dataIdx, fileIdx in enumerate(localSlab[0]):
				curLocalSlab = (fileIdx,) + localSlab[1:]
				localData[dataIdx, :] = dataset[curLocalSlab]
		else:
			raise "Only the first rank may be fancy indexing"

	finally:
		#Make sure file is closed
		f.close()


#--------------------------------------------------------------------------------------
#                         Config Serialization
#--------------------------------------------------------------------------------------

def SaveConfigObject(filename, datasetPath, conf):
	if conf != None:
		h5file = tables.openFile(filename, "r+")
		try:
			dataset = GetExistingDataset(h5file, datasetPath)
			if datasetPath:
				h5file.setNodeAttr(dataset, "configObject", conf.cfgObj)
		except:
			print "WARNING: Could not store config object!"
		finally:
			h5file.close()


def GetConfigFromHDF5(file, datasetPath = None, confObjName = "configObject"):
	"""
	Load a configparser object stored as an attribute on a wavefunction in a 
	HDF5 file. This can be used to recreate a pyprop.Config object, or the
	original config file.
	"""
	h5file = tables.openFile(file, "r")
	cfgObj = None
	try:
		if datasetPath == None:
			for node in h5file.walkNodes():
				#name = node._v_name
				#if name == "wavefunction":
				if confObjName in node._v_attrs._v_attrnames:
					cfgObj = node._f_getAttr(confObjName)
		else:
			cfgObj = h5file.getNodeAttr(datasetPath, confObjName)

	finally:
		h5file.close()
	
	return cfgObj

