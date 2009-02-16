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

def RunHasbani(**args):
	args["silent"] = True
	if not "config" in args:
		args["config"] = "config_hasbani.ini"

	#Set up init problem
	initProp = SetupProblem(**args)
	tempPsi = initProp.psi

	#Find radial eigenstates
	E, V = SetupRadialEigenstates(initProp)
	overlapMatrix = SetupOverlapMatrix(initProp)

	frequencyList = r_[0.7:0.9:0.02]
	frequencyData = []
	for frequencyIndex, frequency in enumerate(frequencyList):
		#setup propagation problem
		potList = ["LaserPotentialVelocity1", "LaserPotentialVelocity2", "LaserPotentialVelocity3", "Absorber"]
		prop = SetupProblem(frequency=frequency, additionalPotentials=potList, **args)

		#Setup initial state
		SetRadialEigenstate(prop.psi, V, n=1, l=0)
		prop.psi.Normalize()
		initPsi = prop.psi.Copy()

		#Propagate
		for t in prop.Advance(False):
			pass

		data = {}
		data["frequency"] = frequency

		#calculate values
		data["norm"] = prop.psi.GetNorm()
		data["corr"] = abs(initPsi.InnerProduct(prop.psi))**2
		data["radialDensity"] = numpy.sum(abs(prop.psi.GetData()), axis=0)

		#Calculate l-distribution
		tempPsi.GetData()[:] = prop.psi.GetData()
		tempPsi.GetRepresentation().MultiplySqrtOverlap(False, tempPsi)
		data["angularDensity"] = numpy.sum(abs(tempPsi.GetData())**2, axis=1)

		#Calculate bound state distribution
		data["boundE"], data["boundV"], data["boundDistr"], data["boundTotal"] = CalculateBoundDistribution(prop.psi, E, V, overlapMatrix, boundThreshold=-0.00)

		data["ionizedE"], data["ionizedDistr"] = CalculateEnergyDistribution(prop.psi, E, V, overlapMatrix)

		frequencyData.append(data)
		del prop

	return frequencyData

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
	PrintOut("Initial State Energy = %s" % (prop.GetEnergy(), ))
	prop.psi.Normalize()
	initPsi = prop.psi.Copy()

	timeList = []
	corrList = []
	normList = []
	radialDensityList = []
	angularDensityList = []
	boundTotalList = []
	outsideAbsorberList = []

	tempPsi = prop.psi.Copy()

	#print "Potentials:\n    %s" % "\n    ".join([pot.Name for pot in prop.Propagator.BasePropagator.PotentialList])

	#Setup absorber tests
	absorberStart = prop.Config.Absorber.absorber_start
	absorberEnd = absorberStart + prop.Config.Absorber.absorber_length
	outsideAbsorberBox = SetupRadialMaskPotential(prop, absorberEnd, 100)
	insideAbsorberBox = SetupRadialMaskPotential(prop, 0, absorberStart)

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()

	#Propagate until end of pulse
	duration = prop.Duration
	#prop.Duration  = 2*pi
	#for t in prop.Advance(False):
	#	pass
	#Remove bound states
	#RemoveBoundDistribution(prop.psi, E, V, overlapMatrix)
	#prop.psi.GetData()[0,:] = 0

	prop.Duration = duration
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
		boundE, boundV, boundDistr, boundTotal = CalculateBoundDistribution(prop.psi, E, V, overlapMatrix)

		#Calculate amount of stuff outside absorber
		outsideAbsorber = abs(outsideAbsorberBox.GetExpectationValue(prop.psi, prop.GetTempPsi(), 0, 0))

		timeList.append(t)
		corrList.append(corr)
		normList.append(norm)
		radialDensityList.append(radialDensity)
		angularDensityList.append(angularDensity)
		boundTotalList.append(boundTotal)
		outsideAbsorberList.append(outsideAbsorber)

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime
	
		#Print stats
		PrintOut("t = %.2f; N = %.10f; Corr = %.10f, outside = %.10f, ETA = %s" % (t, norm, corr, outsideAbsorber, FormatDuration(eta)))

	#Save the time-valued variables
	prop.TimeList = timeList
	prop.CorrList = corrList
	prop.NormList = normList
	prop.RadialDensityList = radialDensityList
	prop.AngularDensityList = angularDensityList
	prop.BoundTotalList = boundTotalList
	prop.OutsideAbsorberList = outsideAbsorberList

	PrintOut("")
	#prop.Propagator.PampWrapper.PrintStatistics()

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

#---------------------------------------------------------------------------------------
#            Box-norm analyisis
#---------------------------------------------------------------------------------------

def SetupRadialMaskPotential(prop, maskStart, maskEnd):
	conf = prop.Config.RadialMaskPotential
	conf.mask_start = maskStart
	conf.mask_end = maskEnd

	pot = prop.Propagator.BasePropagator.GeneratePotential(conf)
	pot.SetupStep(0)

	return pot


