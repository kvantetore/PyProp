import tables

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

	groundstateFilename = args.get("groundstateFilename", "helium_groundstate.h5")
	groundstateDatasetPath = args.get("groundstateDatasetPath", "/wavefunction")

	#Find Groundstate
	findGroundstate = args.get("findGroundstate", True)
	initPsi = None
	if findGroundstate:
		#Find Groundstate with piram
		initProp = SetupProblem(eigenvalueCount=1, **args)
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
	potList = ["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity", "Absorber"]
	prop = SetupProblem(additionalPotentials=potList, **args)
	
	#Setup initial state
	#if initPsi == None:
	#	prop.LoadWavefunctionHDF(groundstateFilename, groundstateDatasetPath)
	#	initPsi = prop.psi.Copy()
	#else:
	prop.psi.GetData()[:] = initPsi.GetData()
	prop.psi.Normalize()
	initPsi = prop.psi.Copy()

	for pot in prop.Propagator.BasePropagator.PotentialList:
		PrintOut( "Potential %s: \n %s" % (pot.Name,  pot.MultiplyFunction.__doc__) )

	timeList = []
	normList = []
	corrList = []

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()
	for t in prop.Advance(outputCount):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2

		#keep these values
		timeList.append(t)
		normList.append(norm)
		corrList.append(corr)

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime
	
		#Print stats
		PrintOut("t = %.2f; N = %.15f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))

	#Final output
	norm = prop.psi.GetNorm()
	corr = abs(initPsi.InnerProduct(prop.psi))**2
	PrintOut("t = %.2f; N = %.10f; Corr = %.10f" % (t, norm, corr))

	PrintOut("")
	prop.Propagator.PampWrapper.PrintStatistics()

	#Saving final wavefunction
	outputFilename = args.get("outputFilename", "final.h5")
	outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
	prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)


	#Save sample times, norm and initial correlation
	if pyprop.ProcId == 0:
		h5file = tables.openFile(outputFilename, "r+")
		try:
			h5file.createArray("/", "SampleTimes", timeList)
			h5file.createArray("/", "Norm", normList)
			h5file.createArray("/", "InitialCorrelation", corrList)
		finally:
			h5file.close()


def SubmitStabilizationRun():
	"""
	Calculate total ionization for a range of intensities to determine stabilization
	"""
	outputDir = "stabilization_freq_5/photoelectron_spectrum_scan"
	frequency = 5.0
	amplitudeList = arange(1,31)
	
	for I in amplitudeList:
		name = "stabilization_I_%i" % I
		Submit(executable="run_stabilization.py", \
			runHours=1, \
			jobname="stabilization", \
			nodes=17, \
			config="config.ini", \
			amplitude=I/frequency, \
			outputCount=300, \
			outputFilename=name)

