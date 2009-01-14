import scipy
import scipy.linalg

#------------------------------------------------------------------------------------
#                       Stabilization functions
#------------------------------------------------------------------------------------

def FormatDuration(duration):
	duration = int(duration)
	h = (duration / 3600)
	m = (duration / 60) % 60
	s = (duration % 60)

	str = []
	if h>0: str += ["%ih" % h]
	if m>0: str += ["%im" % m]
	if s>0: str += ["%is" % s]

	return " ".join(str)

def RunStabilization(**args):
	groundstateFilename = args.get("groundstateFilename", "helium_groundstate.h5")
	groundstateDatasetPath = args.get("groundstateDatasetPath", "/wavefunction")

	#Set up propagation problem
	potList = ["LaserPotentialVelocity1", "LaserPotentialVelocity2", "LaserPotentialVelocity3", "Absorber"]
	prop = SetupProblem(additionalPotentials=potList, **args)

	#Find radial eigenstates
	E, V = SetupRadialEigenstates(prop)
	overlapMatrix = SetupOverlapMatrix(prop)
		
	#Setup initial state
	SetRadialEigenstate(prop.psi, V, n=1, l=0)
	print "Norm = ", prop.psi.Normalize()
	PrintOut("Initial State Energy = %s" % (prop.GetEnergy(), ))
	initPsi = prop.psi.Copy()

	timeList = []
	corrList = []
	normList = []
	radialDensityList = []
	angularDensityList = []
	boundTotalList = []

	tempPsi = prop.psi.Copy()

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()
	for t in prop.Advance(outputCount, yieldEnd=True):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2
		radialDensity = numpy.sum(abs(prop.psi.GetData()), axis=0)

		#Calculate l-distribution
		tempPsi.GetData()[:] = prop.psi.GetData()
		tempPsi.GetRepresentation().MultiplySqrtOverlap(False, tempPsi)
		angularDensity = numpy.sum(abs(tempPsi.GetData())**2, axis=1)

		#Calculate bound state distribution
		boundE, boundDistr, boundTotal = CalculateBoundDistribution(prop.psi, E, V, overlapMatrix)

		timeList.append(t)
		corrList.append(corr)
		normList.append(norm)
		radialDensityList.append(radialDensity)
		angularDensityList.append(angularDensity)
		boundTotalList.append(boundTotal)

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime
	
		#Print stats
		PrintOut("t = %.2f; N = %.10f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))

	#Save the time-valued variables
	prop.TimeList = timeList
	prop.CorrList = corrList
	prop.NormList = normList
	prop.RadialDensityList = radialDensityList
	prop.AngularDensityList = angularDensityList
	prop.BoundTotalList = boundTotalList

	PrintOut("")
	prop.Propagator.PampWrapper.PrintStatistics()

	#Saving final wavefunction
	outputFilename = args.get("outputFilename", "final.h5")
	outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
	prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)

	return prop


def SetupBigMatrix(prop, whichPotentials):
	"""
	Generates a huge dense matrix of a list of potentials in prop.

	whichPotentials is a list of integers specifying the indices of 
	the tensor-potentials to include in the dense matrix
	"""
	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	
	#Allocate the hamilton matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" % (matrixSize, matrixSize, matrixSize**2 * 16 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]

		for i, (x0,x0p) in enumerate(basisPairs0):
			for j, (x1,x1p) in enumerate(basisPairs1):
				indexLeft = (x1 * Count0) + x0
				indexRight = (x1p * Count0) + x0p
				BigMatrix[indexLeft, indexRight] += potential.PotentialData[i, j]

	return BigMatrix

def SetupRadialEigenstates(prop, potentialIndices=[0]):
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
	S = SetupOverlapMatrix(prop)

	eigenValues = []
	eigenVectors = []

	lCount = prop.psi.GetData().shape[0]
	eigenvectorScaling = 1

	for l in range(lCount):
		l = int(l)
		M = SetupRadialMatrix(prop, potentialIndices, l)

		E, V = scipy.linalg.eig(a=M, b=S)

		idx = argsort(real(E))
		eigenValues.append(real(E[idx]))

		#Sort an normalize eigenvectors
		sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
		V = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()
		eigenVectors.append(V)

	return eigenValues, eigenVectors

def SetRadialEigenstate(psi, eigenVectors, n, l):
	"""
	Sets psi to an eigenvector from a list of eigenvectors as calculated by
	SetupRadialEigenstates()

	n, l is the quantum number of the state to set in psi
	"""
	radialIndex = n - l - 1
	vec = eigenVectors[l][:, radialIndex]
	psi.Clear()
	psi.GetData()[l, :] = vec


def CalculateRadialCorrelation(psi, eigenVectors, n, l, overlap):
	radialIndex = n - l - 1
	psi.Clear()
	vec = eigenVectors[l][:, radialIndex]
	return sum( conj(vec) * dot(overlap, psi.GetData()[l, :]) )
	

def CalculateEnergyDistribution(psi, eigenValues, eigenVectors, overlap):
	dE = min(diff(eigenValues[0]))
	E = r_[eigenValues[0][0]:eigenValues[0][-1]:dE]
	energyDistr = []

	for l, (curE, curV) in enumerate(zip(eigenValues, eigenVectors)):
		#Get projection on eigenstates
		psiSlice = psi.GetData()[l, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = abs(dot(curV.transpose(), overlapPsi))**2

		#Calculate density of states
		interiorSpacing = list((diff(curE[:-1]) + diff(curE[1:])) / 2.)
		leftSpacing = (curE[1] - curE[0]) / 2.
		rightSpacing = (curE[-1] - curE[-2]) / 2.
		spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
		density = 1.0 / spacing

		#Interpolate to get equispaced dP/dE
		energyDistr.append( scipy.interp(E, curE, proj * density, left=0, right=0) )

	return E, energyDistr
	
def CalculateBoundDistribution(psi, eigenValues, eigenVectors, overlap):
	boundDistr = []
	boundE = []
	boundTotal = 0

	for l, (curE, curV) in enumerate(zip(eigenValues, eigenVectors)):
		boundIdx = where(curE < 0)[0]

		#Get projection on eigenstates
		psiSlice = psi.GetData()[l, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = abs(dot(curV[:,boundIdx].transpose(), overlapPsi))**2

		#Interpolate to get equispaced dP/dE
		boundDistr.append( proj )
		boundE.append( curE[boundIdx] )
		boundTotal += sum(proj)

	return boundE, boundDistr, boundTotal
	

def SetupRadialMatrix(prop, whichPotentials, angularIndex):
	matrixSize = prop.psi.GetData().shape[1]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

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


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupRadialMatrix(prop, [overlap], 0)
	return matrix
	

