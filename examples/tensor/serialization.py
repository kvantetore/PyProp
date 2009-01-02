#
# SERIALIZATION.PY
#
# Various function which allows tensor potentials to be stored to and retrieved from
# HDF5 files.
#

def SetupTensorPotentialFromFile(psi, fileName, datasetPath, potentialName):
	"""
	Set up a tensorpotential using potential data stored in a HDF5 file.
	"""
	#Create a tensor potential generator 		
	generator = TensorPotentialGenerator(representation = psi.GetRepresentation())

	#Create tensor potential object
	tensorPotential = TensorPotential(psi)

	#Get potential data and configsection from file
	h5file = tables.openFile(fileName, "r")
	try:
		potentialData = h5file.getNode(datasetPath)[:]
		confSectList = h5file.getNodeAttr(datasetPath, "configSection")
	finally:
		h5file.close()
	
	configSection = pyprop.Section(potentialName)
	for opt, val in confSectList:
		configSection.Set(opt, eval(val))

	#Get geometry list
	geometryList = generator.GetGeometryList(configSection)

	#Create PotentialWrapper for TensorPotential
	configSection.Apply(tensorPotential)
	tensorPotential.GeometryList = geometryList
	tensorPotential.PotentialData = potentialData
	tensorPotential.Name = configSection.name.replace("__", "+")

	basisPairs = [geom.GetBasisPairs() for geom in geometryList]

	multiplyFuncName = "pyprop.core.TensorPotentialMultiply_" + "_".join([geom.GetStorageId() for geom in geometryList])
	tensorPotential.MultiplyFunction = eval(multiplyFuncName)

	return tensorPotential


def StoreTensorPotentials(prop, filename, datasetPath="/"):
	"""
	Store all tensor potentials on the propagator object to a HDF5 file
	"""
	potentialList = prop.Propagator.BasePropagator.PotentialList

	h5file = tables.openFile(filename, "w")
	try:
		for pot in potentialList:
			cfgSection = prop.Config.cfgObj.items(pot.Name.split("+")[0])
			potName = pot.Name.replace("+", "__")
			h5file.createArray(datasetPath, potName, pot.PotentialData)
			h5file.setNodeAttr(datasetPath + potName, "configSection", cfgSection)
	finally:
		h5file.close()
