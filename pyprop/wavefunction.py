import core
import distribution as distrib
import representation as repr
from createinstance import CreateInstanceRank

class WavefunctionFileFormat:
	Ascii  = 1
	Binary = 2
	HDF    = 3


def CreateWavefunction(config):
	"""
	Creates a Wavefunction from a config file. Use this function if
	you only need a wavefunction and not a complete proapagator.

	The wavefunction will have one data buffer allocated, and the content
	is unspecified.

	ex:
	conf = pyprop.Load("config.ini")
	psi = CreateWavefunction(config)
	x = psi.GetData().GetRepresentation().GetLocalGrid(0)
	psi.GetData()[:] = x * exp(- x**2)
	"""

	print "Creating DistributionModel..."
	distribution = distrib.CreateDistribution(config)

	print "Creating Representation..."
	representation = repr.CreateRepresentation(config, distribution)

	print "Creating Wavefunction..."
	psi = CreateWavefunctionInstance(representation)

	return psi


def CreateWavefunctionFromFile(filename, datasetPath="/wavefunction"):
	"""
	Loads a wavefunction directly from a HDF5-file. The config object
	is read from the attribute configObject on the node specified by
	datasetPath
	"""
	conf = config.LoadConfigFromFile(filename, datasetPath)
	psi = CreateWavefunction(conf)
	serialization.LoadWavefunctionHDF(filename, datasetPath, psi)

	return psi



def CreateWavefunctionInstance(representation, allocateData=True):
	#Create instance
	print "    Creating instance"
	rank = len(representation.GetFullShape())
	psi = CreateInstanceRank("core.Wavefunction", rank)
	
	#Set reresentation
	print "    Setting representation"
	psi.SetRepresentation(representation)
	
	#Allocate data
	if allocateData:
		print "    Allocating data"
		psi.AllocateData()
	
	return psi
	
	
