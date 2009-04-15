import tables
import sys

from pyprop.utilities import ElectricFieldAtomicFromIntensitySI as field_from_intensity
from pyprop.utilities import AngularFrequencyAtomicFromWavelengthSI as freq_from_wavelength
#pyprop.serialization.DEBUG = True

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
	if "configFile" not in args: args["configFile"] = "config.ini"

	gsFilenamePostfix = args.get("gsFilenamePostfix", "")
	gsFilenameGenerated = "output/initialstates/helium_groundstate_%s_%s.h5" % ("_".join(GetRadialGridPostfix(**args)), gsFilenamePostfix)
	groundstateFilename = args.get("groundstateFilename", gsFilenameGenerated)
	groundstateDatasetPath = args.get("groundstateDatasetPath", "/wavefunction")

	saveWavefunctionDuringPropagation = args.get("saveWavefunctionDuringPropagation", False)
	storeInitialState = args.get("storeInitialState", False)

	#Find Groundstate
	findGroundstate = args.get("findGroundstate", True)
	initPsi = None
	if findGroundstate:
		#Find Groundstate with piram
		initProp = SetupProblem(eigenvalueCount=1, **args)
		#solver = pyprop.PiramSolver(initProp)
		#solver.Solve()

		#Setup inverse iterator	
		invIt = pyprop.GMRESShiftInvertSolver(initProp)
		initProp.Config.Arpack.matrix_vector_func = invIt.InverseIterations

		#Setup solver
		solver = pyprop.PiramSolver(initProp)
		solver.Solve()
	
		#Get groundstate wavefunction
		#initPsi = initProp.psi
		solver.SetEigenvector(initProp.psi, 2)
		initProp.psi.Normalize()
		initPsi = initProp.psi.Copy()

		#Store wavefunction
		if storeInitialState:
			initProp.SaveWavefunctionHDF(groundstateFilename, groundstateDatasetPath)

		#Print ground state energy
		energyExpt = initProp.GetEnergyExpectationValue()
		groundstateEnergy = solver.GetEigenvalues()[0].real
		PrintOut("Ground state energy = %s (%s)" % (1.0 / groundstateEnergy + initProp.Config.GMRES.shift, energyExpt))

		#free memory
		#initPsi = None
		del solver
		del initProp
		del invIt
		
		
	#Set up propagation problem
	if args.get("laserOff", False):
		potList = ["Absorber"]
		PrintOut("Setting up new problem WITHOUT laser potentials...")
	else:
		potList = ["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity", "Absorber"]
		PrintOut("Setting up new problem with laser potentials...")

	sys.stdout.flush()
	prop = SetupProblem(additionalPotentials=potList, **args)
	PrintOut("Done setting up problem! (initPsi = %s)" % initPsi)
	sys.stdout.flush()
	
	#Setup initial state
	if initPsi == None:
		PrintOut("Loading wavefunction... (%s,%s,%s)" % prop.psi.GetData().shape)
		sys.stdout.flush()
		prop.LoadWavefunctionHDF(groundstateFilename, groundstateDatasetPath)
		initPsi = prop.psi.Copy()
	else:
		#PrintOut("Uh-oh, trying to get data from deleted object!")
		sys.stdout.flush()
		prop.psi.GetData()[:] = initPsi.GetData()
		prop.psi.Normalize()
		#initPsi = prop.psi.Copy()

	#Set up box norm potentials
	singleIonizationBox = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.SingleIonizationBox)
	singleIonizationBox.SetupStep(0)
	doubleIonizationBox = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.DoubleIonizationBox)
	doubleIonizationBox.SetupStep(0)

	#Get electron coupling potential
	propPotList = prop.Propagator.BasePropagator.PotentialList
	electronicPotential = propPotList[where([pot.Name == "ElectronicCouplingPotential" for pot in propPotList])[0]]

	tmpPsi = prop.psi.Copy()

	distr = prop.psi.GetRepresentation().GetDistributedModel()

	timeList = []
	normList = []
	corrList = []
	singleIonization = []
	doubleIonization = []
	r12ExpVal = []

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 300)
	startTime = time.time()
	for step, t in enumerate(prop.Advance(outputCount)):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2

		#keep these values
		timeList.append(t)
		normList.append(norm)
		corrList.append(corr)

		#calculate ionization
		tmpPsi.Clear()
		singleIonization += [singleIonizationBox.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]
		tmpPsi.Clear()
		doubleIonization += [doubleIonizationBox.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]

		#calculate expectation value of 1/r12
		tmpPsi.Clear()
		r12ExpVal += [electronicPotential.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]

		#Save wavefunction
		if saveWavefunctionDuringPropagation:
			outputFilename = args.get("outputFilename", "wavefunction") + "_%04i.h5" % step
			outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
			prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime

		#Print stats
		PrintOut("t = %.2f; N = %.15f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))
		PrintOut(prop.Propagator.Solver.GetErrorEstimateList())

	#Final output
	norm = prop.psi.GetNorm()
	corr = abs(initPsi.InnerProduct(prop.psi))**2
	timeList.append(t)
	normList.append(norm)
	corrList.append(corr)

	#calculate ionization
	tmpPsi.Clear()
	singleIonization += [singleIonizationBox.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]
	tmpPsi.Clear()
	doubleIonization += [doubleIonizationBox.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]

	#calculate expectation value of 1/r12
	tmpPsi.Clear()
	r12ExpVal += [electronicPotential.GetExpectationValue(prop.psi, tmpPsi, 0, 0)]

	PrintOut("t = %.2f; N = %.10f; Corr = %.10f" % (t, norm, corr))

	PrintOut("")
	#prop.Propagator.PampWrapper.PrintStatistics()
	prop.Propagator.Solver.PrintStatistics()

	#Saving final wavefunction
	outputFilename = args.get("outputFilename", "final") + ".h5"
	outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
	prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)

	#Save sample times, norm and initial correlation
	if pyprop.ProcId == 0:
		h5file = tables.openFile(outputFilename, "r+")
		try:
			h5file.createArray("/", "SampleTimes", timeList)
			h5file.createArray("/", "Norm", normList)
			h5file.createArray("/", "InitialCorrelation", corrList)
			h5file.createArray("/", "SingleIonization", singleIonization)
			h5file.createArray("/", "DoubleIonization", doubleIonization)
			h5file.createArray("/", "ElectronicCouplingExpectationValue", r12ExpVal)
		finally:
			h5file.close()


def SubmitStabilizationRun():
	"""
	Calculate total ionization for a range of intensities to determine stabilization
	"""
	configFile = "config_stabilization_freq_5.ini"
	timestep = 0.01
	frequency = 3.0
	cycles = 6
	pulse_duration = cycles * 2 * pi / frequency
	#duration = 1 * cycles * 2 * pi / frequency
	duration = 1.0 * pulse_duration

	radialGrid = {\
	'xsize' : 80, \
	'xmax' : 80, \
	'order' : 5, \
	'bpstype' : 'exponentiallinear', \
	'xpartition' : 20, \
	'gamma' : 2.0}

	amplitudeList = arange(1.0, 20.0)
	#amplitudeList = arange(25.0, 31.0)
	#amplitudeList = arange(41.0, 61.0)
	#amplitudeList = [1]
	
	#outputDir = "stabilization_freq_%s_scan_1s2p_%s/" % (frequency, "_".join(GetRadialGridPostfix(config=configFile, radialGrid=radialGrid)))
	outputDir = "stabilization_freq_%s_scan_M1_%s/" % (frequency, "_".join(GetRadialGridPostfix(config=configFile, radialGrid=radialGrid)))
	if not os.path.exists(outputDir):
		print "Created output dir: %s" % outputDir
		os.makedirs(outputDir)
	
	#for I in amplitudeList:
	for I in [10]:
		name = outputDir + "stabilization_I_%i_kb20_dt_%1.e_T_%2.1f" % (I, timestep, duration)
		RunSubmitFullProcCount(RunStabilization, \
		#RunSubmitFullProcCount(SetupAllStoredPotentials, \
		#RunSubmitFullProcCount(SetupIonizationBoxStoredPotentials, \
			account='fysisk', \
			procPerNode=1, \
			procMemory="4000M", \
			walltime=timedelta(hours=3), \
			config=configFile, \
			timestep=timestep, \
			duration=duration,\
			frequency=frequency, \
			amplitude=I/frequency, \
			pulse_duration=pulse_duration, \
			outputCount=300, \
			outputFilename=name, \
			findGroundstate=True, \
			storeInitialState=True, \
			saveWavefunctionDuringPropagation=False, \
			writeScript=False, \
			useStoredPotentials = False, \
			radialGrid=radialGrid, \
			#multipoleCutoff = 5, \
			laserOff = False, \
			gsFilenamePostfix = "_M1", \
			shift = -2.1, \
			)
			

def SubmitHasbaniExampleRun(workingDir):
	"""
	Calculate total ionization for a range of frequencies
	"""
	outputDir = "example_hasbani/"
	#frequencyList = [.2, .3, .45, .5, .6, .65, .7, .75, .8, .9, 1., 1.1, 1.2]
	frequencyList = [.34, .36, .38, .41, .43, 0.71, .72, .74, .74]
	
	for w in frequencyList:
		name = outputDir + "example_hasbani_omega_%i.h5" % (w * 100)
		Submit(executable="run_stabilization.py", \
			runHours=5, \
			jobname="helium", \
			numProcs=111, \
			frequency=w, \
			config="config_hasbani.ini", \
			outputCount=300, \
			workingDir=workingDir, \
			outputFilename=name, \
			findGroundstate = False, \
			writeScript=False)


def SubmitNikolopoulosExampleRun():
	"""
	Calculate single and double ionization for a range of intensities
	"""
	outputDir = "example_nikolopoulos_intensity_scan/"
	
	intensityList = [2e14, 4e14, 8e14, 1e15, 2e15, 4e15, 8e15, 1e16]
	
	for I in intensityList:
		name = outputDir + "example_nikolopoulos_I_%s%s" % tuple(str(I).split("+"))
		RunSubmitFullProcCount(RunStabilization, \
			procPerNode=2, \
			procMemory="2000mb", \
			walltime=timedelta(hours=6), \
			frequency=1.65, \
			amplitude=field_from_intensity(I)/1.65, \
			config="config_largegrid.ini", \
			outputCount=300, \
			outputFilename=name, \
			findGroundstate = True, \
			writeScript=False)


def SubmitStabilizationEigenvaluesJob():
	"""
	Calculate a number of the lowest eigenvalues of Helium for a range of L's
	"""
	lmax = 7
	radialGrid = {\
		'xpartition': 20,\
		'xmax': 200,\
		'bpstype': 'exponential',\
		'xsize': 40,\
		'order': 5,\
		'gamma': 6.0}
	
	#for L in range(lmax+1):
	for L in [0,1,2,3]:
		RunSubmitFullProcCount(RunFindBoundstates, \
			account='fysisk', \
			procPerNode=4, \
			procMemory="1000M", \
			walltime=timedelta(hours=1, minutes=30), \
			config="config_stabilization_freq_5.ini", \
			writeScript=False, \
			radialGrid=radialGrid, \
			L = L, \
			lmax = lmax, \
			shift = L == 0 and -2.9 or -2.1)


def FindStabilization(runFilePath):
	#runFilePath = "stabilization_freq_5_cycle4"
	lmax = 5
	boundstateFiles = ["eigenstates/eigenvalues_stabilization_L%i_20stk.h5" % l for l in range(lmax+1)]

	runFiles = os.listdir(runFilePath)

	for file in runFiles:
		if file[-2:] == "h5":
			dataFile = runFilePath + "/" + file
			print "Processing file: %s" % file
			curIonization = FindIonizationProbability(dataFile, boundstateFiles[:], ionizationThreshhold=-2.0)
			h5file = tables.openFile(dataFile, "r+")
			try:
				h5file.createArray("/", "TotalIonization", [curIonization])
			except:
				"Could not create array! Moving on..."
				pass
			finally:
				h5file.close()
