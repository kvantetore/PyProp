execfile("example.py")

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
	if "configFile" not in args: args["configFile"] = "config_helium.ini"

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
	potList = ["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity"]
	prop = SetupProblem(additionalPotentials=potList, **args)
	
	#Setup initial state
	if initPsi == None:
		prop.LoadWavefunctionHDF(groundstateFilename, groundstateDatasetPath)
		initPsi = prop.psi.Copy()
	else:
		prop.psi.GetData()[:] = initPsi.GetData()

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()
	for t in prop.Advance(outputCount):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime
	
		#Print stats
		PrintOut("t = %.2f; N = %.10f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))

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

