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

	#Find Groundstate
	findGroundstate = args.get("findGroundstate", True)
	initPsi = None
	if findGroundstate:
		#Find Groundstate with piram
		initProp = SetupProblem(eigenvalueCount=2, **args)
		solver = pyprop.PiramSolver(initProp)
		solver.Solve()
	
		#Get groundstate wavefunction
		initPsi = initProp.psi
		solver.SetEigenvector(initPsi, 0)

		initProp.SaveWavefunctionHDF(groundstateFilename, groundstateDatasetPath)

		#free memory
		del solver
		del initProp

	
		
	#Set up propagation problem
	potList = ["LaserPotentialVelocity1", "LaserPotentialVelocity2", "LaserPotentialVelocity3", "Absorber"]
	prop = SetupProblem(additionalPotentials=potList, **args)

	#Find radial eigenstates
	E, V = SetupRadialEigenstates(prop)
		
	#Setup initial state
	if initPsi == None:
		prop.LoadWavefunctionHDF(groundstateFilename, groundstateDatasetPath)
		initPsi = prop.psi.Copy()
	else:
		prop.psi.GetData()[:] = initPsi.GetData()

	timeList = []
	corrList = []
	normList = []
	radialDensityList = []
	angularDensityList = []

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

		timeList.append(t)
		corrList.append(corr)
		normList.append(norm)
		radialDensityList.append(radialDensity)
		angularDensityList.append(angularDensity)

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

	PrintOut("")
	prop.Propagator.PampWrapper.PrintStatistics()

	#Saving final wavefunction
	outputFilename = args.get("outputFilename", "final.h5")
	outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
	prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)

	return prop


def SetupBigMatrix(prop, whichPotentials):
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

def SetupRadialEigenstates(prop):
	S = SetupOverlapMatrix(prop)

	eigenValues = []
	eigenVectors = []

	for l in range(prop.psi.GetData().shape[0]):
		l = int(l)
		M = SetupRadialMatrix(prop, [0], l)

		E, V = scipy.linalg.eig(a=M, b=S)
		idx = argsort(real(E))
		eigenValues.append(real(E[idx]))
		eigenVectors.append(V[:,idx])

	return eigenValues, eigenVectors



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
	

