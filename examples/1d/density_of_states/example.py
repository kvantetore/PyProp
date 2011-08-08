import scipy, scipy.linalg
from numpy import sqrt, pi, zeros, double, array, abs, sum, argsort, real
from numpy import conj, dot, diff
import pyprop

def AddToProjectNamespace(obj):
	pyprop.ProjectNamespace[obj.__name__] = obj
	return obj

def CalculateDensityOfStates(configFile):

	#Set up config
	conf = pyprop.Load(configFile)

	#Set up pyprop problem
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Calculate eigenvalues
	E, V = SetupEigenstates(prop)

	#Calculate DOS
	dos = 1.0/diff(E)

	#Calculate exact dos
	xmax = conf.RadialRepresentation.xmax
	dos_exact = xmax / ( pi *sqrt(2 * E[1:]) )

	#Estimate highest reliable energy, about 2/3*E_max for B-splines
	maxIdx = int(2*len(E) / 3.0)
	maxReliableEnergy = E[maxIdx]
	print "Estimated highest reliable energy = %1.1e" % maxReliableEnergy

	return E[1:], dos_exact, dos


def SetupEigenstates(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	SAE Helium system.

	The eigenvalues are found by setting up a radial matrix for each l-value
	and using the generalized eigenvalue solver in scipy to find all
	eigenvalues and vectors. 
	
	eigenvalues is a list of 1-d eigenvalue arrays. 	
	
	"""

	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)

	pot = prop.Propagator.BasePropagator.PotentialList[0]
	M = SetupHamiltonMatrix(prop.psi, pot)

	E, V = scipy.linalg.eigh(a=M, b=S)

	idx = argsort(real(E))
	E = real(E[idx])
	eigenValues = E

	#Sort and normalize eigenvectors
	sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
	V = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()
	eigenVectors = V

	return eigenValues, eigenVectors


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupHamiltonMatrix(prop.psi, overlap)
	return matrix


def SetupHamiltonMatrix(psi, potential):
	matrixSize = psi.GetData().shape[0]
	matrix = zeros((matrixSize, matrixSize), dtype=double)

	print "Processing potential: %s" % (potential.Name, )

	basisPairs = potential.BasisPairs[0]
	for i, (x,xp) in enumerate(basisPairs):
		indexLeft = x
		indexRight = xp
		matrix[indexLeft, indexRight] += potential.PotentialData[i]

	return matrix


def GetRadialGridInfo(conf):
	"""Display radial grid information

	From a config file, some information of the radial grid
	is produced:

	  -minimum radial spacing
	  -maximum radial spacing
	  -number of radial points

	Input: pyprop config object
	Returns: None

	"""
	distr = pyprop.CreateDistribution(conf)
	repr = pyprop.CreateRepresentation(conf, distr)
	bsplObj = repr.GetRepresentation(0).GetBSplineObject()
	rGrid = bsplObj.GetBreakpointSequence()

	print "Min. r-spacing: %f" % diff(rGrid)[0]
	print "Max. r-spacing: %f" % diff(rGrid)[-1]
	print "Num. of r points: %i" % len(rGrid)



#---------------------------------------------------------------------------------------------------
# Potentials
#---------------------------------------------------------------------------------------------------
@AddToProjectNamespace
class KineticEnergyPotential:
	def __init__(self):
		self.TimeStep = None
		self.CurTime = None

	def ApplyConfigSection(self, conf):
		self.Mass = conf.Get("mass")

#	def GetPotentialValue(self, pos):
#		return -1.0 / (2.0 * self.Mass)

	def UpdatePotentialData(self, data, psi, dt, t):
		data[:] = -1.0 / (2.0 * self.Mass)

