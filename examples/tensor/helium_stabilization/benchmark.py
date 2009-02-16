

import tables
	
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
	
	- datasetName is the path to the dataset where the wavefunction data
	should be stored.

	- psi is the wavefunction

	- other arguments: conf, a pyprop.Config object, which contains a cofig-
	                   parser object. The configparser object is stored as
	                   an attribute on the wavefunction-array.
	
	Example:
	SaveWavefunctionHDF("myfile.h5", "/mygroup/wavefunction", prop.psi):
	"""

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = GetLocalSlab(localShape, distr)

	if distr.IsSingleProc():
		t = - time.time()
		RemoveExistingDataset(filename, datasetPath)
		SaveLocalSlab(filename, datasetPath, localData, distr)
		t += time.time()
		if DEBUG: print "Duration: %.10fs" % t
		SaveConfigObject(filename, datasetPath, conf)

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
				SaveLocalSlab(filename, datasetPath, localData, distr)
				t += time.time()
				if DEBUG: print "    Duration: %.10fs" % t

			distr.GlobalBarrier();
		if procId==0 and DEBUG: print "Done."

		#proc 0 save config object
		if procId == 0:
			SaveConfigObject(filename, datasetPath, conf)
	
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
		try:
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
			#Make sure dataset is closed
			dataset.close()
	
	finally:
		#Make sure file is closed
		f.close()



def GetLocalSlab(localData, distr):
	"""
	Get the local slab for the current processor in the global dataset
	"""

	#get shapes
	localShape = localData.shape
	fullShape = GetGlobalShape(localShape, distr)
	localStart = GetLocalStart(localShape, fullShape, distr)

	#set up hyperslabs
	fileSlab = []
	for i in range(rank):
		start = localStart[i]
		end = start + localShape[i]
		fileSlab += [slice(start, end)]
	fileSlab = tuple(fileSlab)

	return fileSlab

def GetGlobalShape(localShape, distr):
	raise Exception()

