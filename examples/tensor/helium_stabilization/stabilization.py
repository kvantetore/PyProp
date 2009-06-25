import tables
import sys

from pyprop.utilities import ElectricFieldAtomicFromIntensitySI as field_from_intensity
from pyprop.utilities import AngularFrequencyAtomicFromWavelengthSI as freq_from_wavelength

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
	pyprop.PrintMemoryUsage("Before RunStabilization")
	if "configFile" not in args: args["configFile"] = "config.ini"

	initialStateL = args["initialStateL"]
	initialStateIndex = args["initialStateIndex"]

	saveWavefunctionDuringPropagation = args.get("saveWavefunctionDuringPropagation", False)
		
	#Set up propagation problem
	potList = []
	if not args.get("laserOff", False):
		PrintOut("Setting up new problem with laser potentials...")
		potList += ["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity"]
	else:
		PrintOut("Setting up new problem WITHOUT laser potentials...")

	if not args.get("absorberOff", False):
		potList += ["Absorber"]
		PrintOut("Setting up new problem with absorber...")
	else:
		PrintOut("Setting up new problem WITHOUT absorber...")

	sys.stdout.flush()
	pyprop.PrintMemoryUsage("Before SetupProblem")
	prop = SetupProblem(additionalPotentials=potList, **args)
	pyprop.PrintMemoryUsage("After SetupProblem")

	pyprop.PrintMemoryUsage("Before Loading InitialState")
	conf = SetupConfig(**args)
	initialEnergyLoaded = LoadBoundstate(prop.psi, conf, initialStateL, initialStateIndex)
	initPsi = prop.psi.Copy()
	initialEnergyCalculated = prop.GetEnergyExpectationValue()
	PrintOut("Initial State Energy = %s (loaded) %s (calculated)" % (initialEnergyLoaded, initialEnergyCalculated))
	pyprop.PrintMemoryUsage("After Loading InitialState")

	PrintOut("Done setting up problem! (initPsi = %s)" % initPsi)
	sys.stdout.flush()
	
	#Set up box norm potentials
	pyprop.PrintMemoryUsage("Before Box Potential 1")
	singleIonizationBox = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.SingleIonizationBox)
	singleIonizationBox.SetupStep(0)
	pyprop.PrintMemoryUsage("Before Box Potential 2")
	doubleIonizationBox = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.DoubleIonizationBox)
	doubleIonizationBox.SetupStep(0)
	pyprop.PrintMemoryUsage("After Box Potential")

	#Get electron coupling potential
	propPotList = prop.Propagator.BasePropagator.PotentialList
	electronicPotential = propPotList[where([pot.Name == "ElectronicCouplingPotential" for pot in propPotList])[0]]
	
	pyprop.PrintMemoryUsage("Before Wavefunction Copy")
	tmpPsi = prop.psi.Copy()
	pyprop.PrintMemoryUsage("After Wavefunction Copy")

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
		pyprop.PrintMemoryUsage("At t = %s" % t)

	#Final output
	pyprop.PrintMemoryUsage("After Propagation")
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

	pyprop.PrintMemoryUsage("After RunStabilization")

def SubmitStabilizationRun():
	"""
	Calculate total ionization for a range of intensities to determine stabilization
	"""
	configFile = "config_stabilization_freq_5.ini"
	timestep = 0.01
	frequency = 5.0
	cycles = 6
	phase = "zero"
	pulse_duration = cycles * 2 * pi / frequency
	duration = 1.5 * pulse_duration
	#duration = 1.0
	lmax = 5
	L = range(lmax+1)

	radialGrid = {\
		'xsize' : 80, \
		'xmax' : 80, \
		'order' : 5, \
		'bpstype' : 'exponentiallinear', \
		'xpartition' : 20, \
		'gamma' : 2.0, \
	}

	#amplitudeList = arange(1.0, 31.0)
	#amplitudeList = arange(25.0, 31.0)
	#amplitudeList = arange(41.0, 61.0)

	#amplitudeList = [1,5,10,15,20]
	amplitudeList = [10]

	statePostfix = "1s2p_1"
	args = {
		'account':'matematisk', 
		'preconditionType': 'ifpack', 
		'procPerNode':1, 
		'procMemory':"4000M", 
		'config':configFile, 
		'storeInitialState':True, 
		'saveWavefunctionDuringPropagation':True, 
		'writeScript':False, 
		'useStoredPotentials': False, 
		'radialGrid':radialGrid, 
		'lmax':lmax, 
		'L': L, 
		}

	radialPostfix = "_".join(GetRadialGridPostfix(**args))
	angularPostFix = "_".join(GetAngularGridPostfix(**args))

	#Set initial state to ground state or excited states
	if statePostfix == "1s1s_1":
		args["initialStateL"] = 0
		args["initialStateIndex"] = 0
	elif statePostfix == "1s2s_3":	
		args["initialStateL"] = 0
		args["initialStateIndex"] = 1
	elif statePostfix == "1s2p_3":	
		args["initialStateL"] = 1
		args["initialStateIndex"] = 0
	elif statePostfix == "1s2p_1":	
		args["initialStateL"] = 1
		args["initialStateIndex"] = 1
	else:
		raise Exception("Unknown statePostfix %s" % (statePostfix))

	outputDir = "output/freq_%s_%s_%s_%s/" % (frequency, radialPostfix, angularPostFix, statePostfix)
	if not os.path.exists(outputDir):
		print "Created output dir: %s" % outputDir
		os.makedirs(outputDir)


	#Check that boundstate files exists
	boundstateArgs = dict(args)
	boundstateArgs["L"] = args["initialStateL"]
	boundstateFile = GetBoundstatesFilename(**boundstateArgs)
	if not os.path.exists(boundstateFile):
		raise Exception("Required boundstate File %s does not exist. Run SetupBoundstates first." % boundstateFile)

	depend = None
	#if not os.path.exists(GetInitialStateFilename(**args)):
	#	dependJob = RunSubmitFullProcCount(RunStabilizationInputFile, walltime=timedelta(hours=0.5), **args)
	#	depend = "afterany:%s" % (dependJob,)

	for I in amplitudeList:
		name = outputDir + "stabilization_I_%i_kb20_dt_%1.e_T_%2.1f_phase_%s" % (I, timestep, duration, phase)
		#args['walltime'] = timedelta(hours=6)
		args['walltime'] = timedelta(hours=10)
		args['timestep'] = timestep
		args['duration'] = duration
		args['frequency'] = frequency
		args['phase'] = phase
		args['amplitude'] = I/frequency
		args['pulse_duration'] = pulse_duration
		args['outputCount'] = 30
		args['outputFilename'] = name
		args['findGroundstate'] = False
		args['absorberOff'] = False
		args['laserOff'] = False
		args['depend'] = depend
		RunSubmitFullProcCount(RunStabilization, **args)

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
	lmax = 5
	radialGrid = {\
		'xsize' : 22, \
		'xmax' : 200, \
		'order' : 5, \
		'bpstype' : 'exponentiallinear', \
		'xpartition' : 8.63130809996, \
		'gamma' : 2.19549141339}

	
	#for L in range(lmax+1):
	#for L in [0,1,2,3]:
	for L in [0]:
		RunSubmitFullProcCount(RunFindBoundstates, \
			account='fysisk', \
			procPerNode=1, \
			procMemory="4000M", \
			walltime=timedelta(hours=24, minutes=00), \
			config="config_stabilization_freq_5.ini", \
			writeScript=False, \
			radialGrid=radialGrid, \
			L = L, \
			lmax = lmax, \
			shift = L == 0 and -2.9 or -2.1, \
			eigenvalueCount = 30, \
			eigenvalueBasisSize = 90)


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
