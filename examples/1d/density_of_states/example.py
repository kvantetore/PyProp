from numpy.linalg import eigh
from numpy import sqrt, pi
import pyprop


def CalculateDensityOfStates(configFile):

	#Set up config
	conf = pyprop.Load(configFile)

	#Set up pyprop problem
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Calculate eigenvalues
	E, V = SetupEigenstates(prop)

	#Calculate DOS
	dos = 1.0/diff(E[0])

	#Calculate exact dos
	xmax = conf.RadialRepresentation.xmax
	dos_exact = xmax / ( pi *sqrt(2 * E[0][1:]) )

	return E[1:], dos_exact, dos


def SetupEigenstates(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	SAE Helium system.

	The eigenvalues are found by setting up a radial matrix for each l-value
	and using the generalized eigenvalue solver in scipy to find all
	eigenvalues and vectors. 
	
	eigenvalues is a list of 1-d eigenvalue arrays. Each array corresponding
	to a
	"""
	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)
	bspl = prop.psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=double)

	eigenValues = []
	eigenVectors = []

	lCount = prop.psi.GetData().shape[0]
	eigenvectorScaling = 1

	M = SetupRadialMatrix(prop, potentialIndices, l)

	E, V = scipy.linalg.eigh(a=M, b=S)

	idx = argsort(real(E))
	E = real(E[idx])
	eigenValues.append(E)

	#Sort and normalize eigenvectors
	sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
	V = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()
	eigenVectors.append(V)

	#assure correct phase convention (first oscillation should start out real positive)
	#for i, curE in enumerate(E):
	#	bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
	#	phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
	#	V[:,i] *= exp(-1.0j * phase)

	return eigenValues, eigenVectors


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupRadialMatrix(prop, [overlap])
	return matrix


def SetupHamiltonMatrix(prop, whichPotentials):
	matrixSize = prop.psi.GetData().shape[0]
	matrix = zeros((matrixSize, matrixSize), dtype=double)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		angularBasisPairs = potential.BasisPairs[0]
		idx = [idx for idx, (i,j) in enumerate(zip(angularBasisPairs[:,0], angularBasisPairs[:,1])) if i==j==angularIndex]
		if len(idx) != 1:
			raise "Invalid angular indices %s" % idx
		idx = idx[0]

		basisPairs = potential.BasisPairs[1]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[idx, i]

	return matrix
