
def SaveWavefunctionHDF(hdfFile, datasetPath, psi, conf=None):
	"""
	Saves the wavefunction data to a dataset in a HDF file, and stores a config object
	as an attribute on the wavefunction if applicable.

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
		SaveLocalWavefunctionSlab(filename, datasetPath, psi)
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
				SaveLocalWavefunctionSlab(filename, datasetPath, psi)
				t += time.time()
				if DEBUG: print "    Duration: %.10fs" % t

			distr.GlobalBarrier();
		if procId==0 and DEBUG: print "Done."

		#proc 0 save config object
		if procId == 0:
			SaveConfigObject(filename, datasetPath, conf)

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
		LoadLocalWavefunctionSlab(filename, datasetPath, psi)

	else:
		#let the processors save their part one by one
		procCount = distr.ProcCount
		procId = distr.ProcId
		for i in range(procCount):
			if procId == i:
				if DEBUG: print "Loading slab in proc ", procId
				LoadLocalWavefunctionSlab(filename, datasetPath, psi)
			distr.GlobalBarrier();
	

def SaveLocalWavefunctionSlab(filename, datasetPath, psi):
	#get hyperslab for this proc
	fullShape = tuple(psi.GetRepresentation().GetFullShape())
	fileSlab = GetFileSlab(psi)
	#save data
	SaveLocalSlab(filename, datasetPath, psi.GetData(), fileSlab, fullShape)


def LoadLocalWavefunctionSlab(filename, datasetPath, psi):
	#get hyperslab for this proc
	fullShape = tuple(psi.GetRepresentation().GetFullShape())
	fileSlab = GetFileSlab(psi)
	#load data	
	LoadLocalSlab(filename, datasetPath, psi.GetData(), fileSlab, fullShape)


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



