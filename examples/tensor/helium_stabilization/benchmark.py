

import tables
from pyprop.serialization import GetExistingDataset, CreateDataset, SaveConfigObject
	
def SaveTensorPotential(filename, pot):
	f = tables.openFile(filename, "w")
	try:
		#Standard parameters
		argList = [pot.PotentialData, 1.0, pot.psi.GetData(), pot.psi.GetData()]
		#Parameters for each geometry info
		for i, geom in enumerate(pot.GeometryList):
			argList += geom.GetMultiplyArguments(pot.psi)

		for i, arg in enumerate(argList):
			if asarray(arg).size != 0:
				f.createArray(f.root, "argument_%i" % i, arg)

	finally:
		f.close()

def GenerateBenchmarkPotential(outputFile=None, potentialList=["ElectronicCouplingPotential"], **args):
	args["useDefaultPotentials"] = False
	args["additionalPotentials"] = potentialList
	args["config"] = "config_benchmark.ini"

	if outputFile == None:
		outputFile = "benchmark/potential_%i.h5" % pyprop.ProcId

	prop = SetupProblem(**args)

	potList = prop.Propagator.BasePropagator.PotentialList
	if len(potList) != 1:
		raise Exception("More than one potential generated, only specify potentials that can be consolidated")

	SaveTensorPotential(outputFile, potList[0])

DEBUG = False

def SaveTensorPotential(filename, datasetPath, potential, distributedModel, conf=None):
	"""
	Saves the wavefunction data to a dataset in a HDF file, and stores a config object
	as an attribute on the wavefunction if applicable.

	Parameters
	--------------------------------------------------------------------
		filename is a string containing a filename to a hdffile which 
		may or may not exist. If it exists, no other process my access 
		the file while this function is called
	
	- datasetPath is the path to the dataset where the wavefunction data
	should be stored.

	- psi is the wavefunction

	- other arguments: conf, a pyprop.Config object, which contains a cofig-
	                   parser object. The configparser object is stored as
	                   an attribute on the potential-array.
	
	Example:
	SaveTensorPotential("myfile.h5", "/mygroup/potential", pot):
	"""

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = tuple(map(GetLocalIndices, potential.GeometryList))
	fullShape = tuple(distr.GetGlobalShape(array(localShape)))

	if distr.IsSingleProc():
		t = - time.time()
		RemoveExistingDataset(filename, datasetPath)
		SaveLocalSlab(filename, datasetPath, localData, localSlab, fullShape)
		t += time.time()
		if DEBUG: print "Duration: %.10fs" % t
		SaveConfigObject(filename, datasetPath, conf)
		SaveGeometryInfo(filename, datasetPath, potential.GeometryList)

	else:
		#let the processors save their part one by one
		procCount = distr.ProcCount
		procId = distr.ProcId
		localSize = localData.nbytes / 1024.**2
		distr.GlobalBarrier()
		if procId==0 and DEBUG: print "Saving %s..." % filename
		for i in range(procCount):
			if procId == i:
				if procId == 0:
					RemoveExistingDataset(filename, datasetPath)

				if DEBUG: print "    Process %i writing hyperslab of %iMB" % (procId, localSize)
				t = - time.time()
				SaveLocalSlab(filename, datasetPath, localData, localSlab, fullShape)
				t += time.time()
				if DEBUG: print "    Duration: %.10fs" % t

			distr.GlobalBarrier();
		if procId==0 and DEBUG: print "Done."

		#proc 0 save attribs
		if procId == 0:
			#config object
			SaveConfigObject(filename, datasetPath, conf)
			#basis pairs
			SaveGeometryInfo(filename, datasetPath, potential.GeometryList)


		#Make sure everyone is finished
		distr.GlobalBarrier();


def LoadTensorPotential(filename, datasetPath, potential, distributedModel):
	"""
	Loads the tensor potential data from a dataset in a HDF file
	"""

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = tuple(map(GetLocalIndices, potential.GeometryList))
	fullShape = tuple(distr.GetGlobalShape(array(localShape)))

	if not CheckLocalSlab(filename, datasetPath, potential.GeometryList, localSlab):
		raise Exception("Organization of basis pairs has changed from when this potential was generated. Please regenerate potential")

	if distr.IsSingleProc():
		t = - time.time()
		LoadLocalSlab(filename, datasetPath, localData, localSlab, fullShape)
		t += time.time()
		if DEBUG: print "Duration: %.10fs" % t

	else:
		#let the processors save their part one by one
		procCount = distr.ProcCount
		procId = distr.ProcId
		localSize = localData.nbytes / 1024.**2
		distr.GlobalBarrier()
		if procId==0 and DEBUG: print "Saving %s..." % filename
		for i in range(procCount):
			if procId == i:
				if DEBUG: print "    Process %i writing hyperslab of %iMB" % (procId, localSize)
				t = - time.time()
				LoadLocalSlab(filename, datasetPath, localData, localSlab, fullShape)
				t += time.time()
				if DEBUG: print "    Duration: %.10fs" % t

			distr.GlobalBarrier();
		if procId==0 and DEBUG: print "Done."

		#Make sure everyone is finished
		distr.GlobalBarrier();


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
		dataset[localSlab] = localData

	
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
	f = tables.openFile(filename, "a")
	try:
		#get dataset
		localShape = localData.shape
		dataset = GetExistingDataset(f, datasetPath)
		if dataset == None:
			#create new dataset
			raise Exception("Dataset '%s' not found in file '%s'" % (datasetPath, filename))
		else:
			#check that is has the correct size
			if not dataset.shape == fullShape:
				raise Exception("Invalid shape on existing dataset. Got %s, expected %s" % (dataset.shape, fullShape))
		
		#write data
		localData[:] = dataset[localSlab]

	finally:
		#Make sure file is closed
		f.close()


def GetLocalIndices(geomInfo):
	"""
	Gets the local slab of a 
	"""
	basisPairIndices = geomInfo.GetLocalBasisPairIndices()
	if (diff(basisPairIndices) == 1).all():
		basisPairIndices = s_[basisPairIndices[0]:basisPairIndices[0]+len(basisPairIndices)]
	return basisPairIndices


def CheckLocalSlab(filename, datasetPath, geometryList, localSlab):
	rank = len(localSlab)
	f = tables.openFile(filename, "a")
	try:
		node = GetExistingDataset(f, datasetPath)
		globalBasisPairList = map(lambda i: getattr(node._v_attrs, "basisPairs%i" % i), r_[:rank])
	finally:
		f.close()

	for origGlobalBasisPairs, localIndices, geomInfo in zip(globalBasisPairList, localSlab, geometryList):
		newBasisPairs = geomInfo.GetGlobalBasisPairs()[localIndices]
		origBasisPairs = origGlobalBasisPairs[localIndices]
		if not (origBasisPairs == newBasisPairs).all():
			return False
	return True

	
def SaveGeometryInfo(filename, datasetPath, geometryList):
	f = tables.openFile(filename, "a")
	try:
		node = GetExistingDataset(f, datasetPath)
		for i, geom in enumerate(geometryList):
			setattr(node._v_attrs, "basisPairs%i" % i, geom.GetGlobalBasisPairs())
	finally:
		f.close()



