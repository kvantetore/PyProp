"""
(De)Serialization routines for Tensor Potentials

use the methods SaveTensorPotential and LoadTensorPotential
to load and save tensor potentials instead of regenerating
them every time
"""


def SaveTensorPotential(filename, datasetPath, potential, distributedModel, conf=None):
	"""
	Saves a tensor potential to a dataset in a HDF file, and stores a config object
	as an attribute.

	Parameters
	--------------------------------------------------------------------
	- filename is a string containing a filename to a hdffile which 
	may or may not exist. 
	
	- datasetPath is the path to the dataset where the wavefunction data
	should be stored.

	- potential is the TensorPotential object to be saved

	- distributedModel is a distributed model from a wavefunction decribing 
	the current distribution of things
	
	"""

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = tuple(map(GetLocalBasisPairSlice, potential.GeometryList))
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
	localSlab = tuple(map(GetLocalBasisPairSlice, potential.GeometryList))
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


def GetLocalBasisPairSlice(geomInfo):
	"""
	Gets the local slab of a 
	"""
	basisPairIndices = asarray(geomInfo.GetLocalBasisPairIndices(), dtype=int)
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

	for i, (origGlobalBasisPairs, localIndices, geomInfo) in enumerate(zip(globalBasisPairList, localSlab, geometryList)):
		newBasisPairs = geomInfo.GetGlobalBasisPairs()[localIndices]
		origBasisPairs = origGlobalBasisPairs[localIndices]
		if not numpy.all(origBasisPairs == newBasisPairs):
			newLen = len(geomInfo.GetGlobalBasisPairs())
			origLen = len(origGlobalBasisPairs)
			print "BASISPAIRDIFF(%i) (%i, %i) = %s != %s" % (i, newLen, origLen, newBasisPairs, origBasisPairs)
			return False
	return True

	
def SaveGeometryInfo(filename, datasetPath, geometryList):
	f = tables.openFile(filename, "a")
	try:
		node = GetExistingDataset(f, datasetPath)
		for i, geom in enumerate(geometryList):
			setattr(node._v_attrs, "basisPairs%i" % i, geom.GetGlobalBasisPairs())
			setattr(node._v_attrs, "storageId%i" % i, geom.GetStorageId())
	finally:
		f.close()



