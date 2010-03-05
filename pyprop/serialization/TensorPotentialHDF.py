"""
(De)Serialization routines for Tensor Potentials

use the methods SaveTensorPotential and LoadTensorPotential
to load and save tensor potentials instead of regenerating
them every time
"""


def SaveTensorPotential(filename, groupPath, potential, distributedModel, conf=None):
	"""
	Saves a tensor potential to a dataset in a HDF file, and stores a config object
	as an attribute.

	Parameters
	--------------------------------------------------------------------
	- filename is a string containing a filename to a hdffile which 
	may or may not exist. 
	
	- groupPath is the path to the group where the tensor potential
	should be stored.

	- potential is the TensorPotential object to be saved

	- distributedModel is a distributed model from a wavefunction decribing 
	the current distribution of things
	
	"""

	if groupPath == "/":
		raise Exception("Storing tensor potentials on '/' is a bad idea (trust us, it is)!")

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = tuple(map(GetLocalBasisPairSlice, potential.GeometryList))
	fullShape = tuple(distr.GetGlobalShape(array(localShape)))

	datasetPath = groupPath + "/potential"

	if distr.IsSingleProc():
		t = - time.time()
		RemoveExistingDataset(filename, groupPath)
		SaveLocalSlab(filename, datasetPath, localData, localSlab, fullShape)
		t += time.time()
		if DEBUG: print "Duration: %.10fs" % t
		SaveConfigObject(filename, groupPath, conf)
		SaveGeometryInfo(filename, groupPath, potential.GeometryList)

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
					RemoveExistingDataset(filename, groupPath)

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
			SaveConfigObject(filename, groupPath, conf)
			#basis pairs
			SaveGeometryInfo(filename, groupPath, potential.GeometryList)


		#Make sure everyone is finished
		distr.GlobalBarrier();


def LoadTensorPotential(filename, groupPath, potential, distributedModel):
	"""
	Loads the tensor potential data from a dataset in a HDF file
	"""

	if groupPath == "/":
		raise Exception("Storing tensor potentials on '/' is a bad idea (trust us, it is)!")

	datasetPath = groupPath + "/potential"

	distr = distributedModel
	localData = potential.PotentialData
	localShape = localData.shape
	localSlab = tuple(map(GetLocalBasisPairSlice, potential.GeometryList))
	fullShape = tuple(distr.GetGlobalShape(array(localShape)))

	if not CheckLocalSlab(filename, groupPath, potential.GeometryList, localSlab):
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


def CheckLocalSlab(filename, groupPath, geometryList, localSlab):
	rank = len(localSlab)
	f = tables.openFile(filename, "r")
	try:
		node = GetExistingDataset(f, groupPath)
		globalBasisPairList = map(lambda i: GetExistingDataset(f, groupPath + "/basisPairs%i" % i)[:], r_[:rank])
	finally:
		f.close()

	for i, (origGlobalBasisPairs, localIndices, geomInfo) in enumerate(zip(globalBasisPairList, localSlab, geometryList)):
		newBasisPairs = geomInfo.GetGlobalBasisPairs()[localIndices]
		origBasisPairs = origGlobalBasisPairs[localIndices]
		if not numpy.all(newBasisPairs == origBasisPairs):
			newLen = len(geomInfo.GetGlobalBasisPairs())
			origLen = len(origGlobalBasisPairs)
			print "BASISPAIRDIFF(%i) (%i, %i) = %s != %s" % (i, newLen, origLen, newBasisPairs, origBasisPairs)
			return False
	return True

	
def SaveGeometryInfo(filename, groupPath, geometryList):
	f = tables.openFile(filename, "a")
	try:
		node = GetExistingDataset(f, groupPath)
		for i, geom in enumerate(geometryList):
			curBasisPairNode = f.createArray(node, "basisPairs%i" % i, geom.GetGlobalBasisPairs())
			setattr(curBasisPairNode._v_attrs, "storageId%i" % i, geom.GetStorageId())
	finally:
		f.close()



