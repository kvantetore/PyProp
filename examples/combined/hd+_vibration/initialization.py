
def SetupConfig(**args):
	configSilent = args.get("configSilent", False)
	pyprop.Redirect.Enable(configSilent)
	try:
		#Decide which config file to use
		configFile = "config.ini"
		if "config" in args:
			configFile = args["config"]
		
		#Load the config file
		conf = pyprop.Load(configFile)
		
		if "dt" in args:
			dt = args["dt"]
			conf.Propagation.timestep = dt
			print "Using timestep %.4f" % dt
		
		#Modify the config
		if "imtime" in args:
			imtime = args["imtime"]
			if imtime: print "Using imaginary time"
			propSection = conf.Propagation
			dt = abs(propSection.timestep)
			renormalize = False
			if imtime:
				dt = -1.0j * dt
				renormalize = True
		
			propSection.timestep = dt
			propSection.renormalization = renormalize
		
		if 'species' in args:	
			species = args["species"]
			conf.VibrationalPotential.species = species
			print "Using species %s" % species
		
		if "molecule" in args:
			molecule = args["molecule"].lower()
			print "Using molecule %s" % molecule
			if molecule == "h2+":
				conf.RadialPropagator.mass = ReducedMass(mass_proton_au, mass_proton_au)
				conf.ProbePulsePotential.static_dipole_moment = 0
			elif molecule == "d2+":
				conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_deuteron_au)
				conf.ProbePulsePotential.static_dipole_moment = 0
			elif molecule == "hd+":
				conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_proton_au)
				conf.ProbePulsePotential.static_dipole_moment = 1. / 3.
			else:
				raise Exception("Unknown molecule '%s'" % molecule)
		
		if "staticDipoleMoment" in args:
			staticDipoleMoment = args["staticDipoleMoment"]
			print "Using static dipole moment %.4f" % staticDipoleMoment
			conf.ProbePulsePotential.static_dipole_moment = staticDipoleMoment
		
		if 'radialScaling' in args:
			radialScaling = args["radialScaling"]
			print "Using radial scaling %i" % radialScaling
			rank0 = conf.RadialRepresentation.rank0
			radialRange = rank0[1] - rank0[0]
			radialSize = rank0[2]
		
			rank0[1] = radialRange / abs(radialScaling) + rank0[0]
			if radialScaling > 0:
				rank0[2] = radialSize / abs(radialScaling)
		
		if "radialGrid" in args:
			radialGrid = args["radialGrid"]
			print "Using radial grid %s" % radialGrid
			conf.RadialRepresentation.rank0 = list(radialGrid)
			
		#Probe Pulse
		if 'pulseDelay' in args:
			pulseDelay = args["pulseDelay"]
			print "Using pulse delay %.4f fs" % (pulseDelay / femtosec_to_au)
			conf.ProbePulsePotential.delay = pulseDelay
		
		if "pulsePhase" in args:
			pulsePhase = args["pulsePhase"]
			print "Using pulse phase %.4f" % pulsePhase
			conf.ProbePulsePotential.phase = pulsePhase
		
		if "pulseDuration" in args:
			pulseDuration = args["pulseDuration"]
			print "Using pulse duration %i fs" % (pulseDuration / femtosec_to_au)
			conf.ProbePulsePotential.duration = pulseDuration
		
		if "pulseIntensity" in args:
			pulseIntensity = args["pulseIntensity"]
			print "Using pulse intensity %s W/cm^2" % pulseIntensity
			conf.ProbePulsePotential.intensity = pulseIntensity
		
		#Control Pulse
		if 'controlDelay' in args:
			controlDelay = args["controlDelay"]
			print "Using control pulse delay %.4f fs" % (controlDelay / femtosec_to_au)
			conf.ControlPulsePotential.delay = controlDelay
		
		if "controlPhase" in args:
			controlPhase = args["controlPhase"]
			print "Using control pulse phase %.4f" % controlPhase
			conf.ControlPulsePotential.phase = controlPhase
		
		if "controlDuration" in args:
			controlDuration = args["controlDuration"]
			print "Using control pulse duration %i fs" % (controlDuration / femtosec_to_au)
			conf.ControlPulsePotential.duration = controlDuration
		
		if "controlIntensity" in args:
			controlIntensity = args["controlIntensity"]
			print "Using control pulse intensity %s W/cm^2" % controlIntensity
			conf.ControlPulsePotential.intensity = controlIntensity
		
		#Propagation
		if "duration" in args:
			duration = args["duration"]
			print "Propagating for %i fs" % (duration / femtosec_to_au)
			conf.Propagation.duration = duration
		
		if "silent" in args:
			silent = args["silent"]
			conf.Propagation.silent = silent
		
		#Pump Pulse
		if "pumpFrequency" in args:
			print "Using pump frequency %.4f" % args["pumpFrequency"]
			conf.InitialCondition.pump_frequency = args["pumpFrequency"]
		
		if "pumpProbability" in args:
			print "Using pump probability %.4f" % args["pumpProbability"]
			conf.InitialCondition.pump_probability = args["pumpProbability"]
		
		if "pumpCount" in args:
			print "Using pump count %.4f" % args["pumpCount"]
			conf.InitialCondition.pump_count = args["pumpCount"]
		
		if "potentialSlope" in args:
			conf.VibrationalPotential.potential_slope=args["potentialSlope"]
	
	finally:
		pyprop.Redirect.Disable()

	return conf

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def FindGroundstate(**args):
	args['config'] = "groundstate.ini"
	args['imtime'] = True
	prop = SetupProblem(**args)

	silent = False
	if "silent" in args:
		silent = args["silent"]

	outputCount = 10
	if silent == True:
		outputCount = 1

	for t in prop.Advance(outputCount):
		E = prop.GetEnergy()
		print "t = %f, E = %.17f" % (t/femtosec_to_au, E)

	E = prop.GetEnergy()
	print "Ground State Energy = %f" % E

	return prop

def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	
	for E in solver.GetEigenvalues():
		print "%.17f" % E

	return solver
	

def GetInputFile(**args):
	if not "molecule" in args:
		raise Exception("Please specify molecule")

	conf = SetupConfig(**args)
	molecule = args["molecule"]
	gridSize = conf.RadialRepresentation.rank0[1]
	boxSize = conf.RadialRepresentation.rank0[2]

	return "inputfiles/%s_%i_%i.h5" % (molecule, gridSize, boxSize)


def SetupEigenstates(**args):
	args["config"] = "groundstate.ini"
	args["species"] = "Ion"
	outputfile = args["outputfile"]

	solver = FindEigenvalues(**args)
	E = solver.GetEigenvalues()
		
	prop = solver.BaseProblem
	shape = (len(E), ) + prop.psi.GetData().shape

	pyprop.serialization.RemoveExistingDataset(outputfile, "/eigenstates")
	pyprop.serialization.RemoveExistingDataset(outputfile, "/eigenvalues")

	f = tables.openFile(outputfile, "a")
	try:
		f.createArray(f.root, "eigenvalues", E)
		atom = tables.ComplexAtom(itemsize=16)
		eigenstates = f.createCArray(f.root, "eigenstates", atom, shape)
		
		for i in range(len(E)):
			solver.SetEigenvector(prop.psi, i)
			eigenstates[i,:] = prop.psi.GetData()

	finally:
		f.close()

def SetupInitialState(**args):
	outputfile = args["outputfile"]

	#Setup initial state
	prop = FindGroundstate(**args)
	prop.SaveWavefunctionHDF(outputfile, "/initial_state")
	f = tables.openFile(outputfile, "a")
	try:
		E = prop.GetEnergy()
		SaveArray(f, "/", "initial_energy", [E])
	finally:
		f.close()

def GetInitialStateEnergy(outputfile):
	f = tables.openFile(outputfile, "r")
	try:
		E = f.root.initial_energy[:][0]
	finally:
		f.close()
	return E


def GetInputFile(**args):
	if not "molecule" in args:
		raise Exception("Please specify molecule")

	conf = SetupConfig(**args)
	molecule = args["molecule"]
	gridSize = conf.RadialRepresentation.rank0[1]
	boxSize = conf.RadialRepresentation.rank0[2]

	return "inputfiles/%s_%i_%i.h5" % (molecule, gridSize, boxSize)


def GetHamiltonMatrix(**args):
	molecule = args["molecule"]
	args["config"] = "groundstate.ini"
	args["species"] = "ion"
	prop = SetupProblem(**args)

	size = prop.psi.GetData().size
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()
	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()
		
	return matrix
	
def GetFullSpectrum(**args):
	print "Setting up hamilton matrix"
	H = GetHamiltonMatrix(**args)
	print "Diagonalizing matrix"
	E, V = eig(H)
	print "Done"

	prop = SetupProblem(**args)
	dr = prop.psi.GetRepresentation().GetRepresentation(0).GetRange(0).Dx

	E = real(E)
	index = argsort(E)
	E = E[index]
	V = V[:,index] / sqrt(dr)
	
	return E, V

def SetupFullSpectrum(**args):
	outputfile = args["outputfile"]

	#For the 1s_sigma_g
	args["potentialSlope"] = 1
	E1, V1 = GetFullSpectrum(**args)

	#For the 2p_sigma_u
	args["potentialSlope"] = 2
	E2, V2 = GetFullSpectrum(**args)
	f = tables.openFile(outputfile, "a")

	try:
		SaveArray(f, "/complete/binding", "eigenvalues", E1)
		SaveArray(f, "/complete/binding", "eigenstates", V1)
		SaveArray(f, "/complete/unbinding", "eigenvalues", E2)
		SaveArray(f, "/complete/unbinding", "eigenstates", V2)

	finally:
		f.close()


def SetupInputFile(**args):
	outputfile = GetInputFile(**args)
	args["outputfile"] = outputfile

	setupEigenstates=False
	if "setupEigenstates" in args:
		setupEigenstates = args["setupEigenstates"]
	setupFullSpectrum=True
	if "setupFullSpectrum" in args:
		setupFullSpectrum = args["setupFullSpectrum"]
	setupInitialState=True
	if "setupInitialState" in args:
		setupInitialState = args["setupInitialState"]

	if setupEigenstates: SetupEigenstates(**args)
	if setupFullSpectrum: SetupFullSpectrum(**args)
	if setupInitialState: SetupInitialState(**args)


def LoadEigenstates(**args):
	inputfile = GetInputFile(**args)

	conf = SetupConfig(**args)
	c = conf.RadialRepresentation
	dr = (c.rank0[1] - c.rank0[0]) / float(c.rank0[2])

	f = tables.openFile(inputfile, "r")
	try:
		E1 = f.root.complete.binding.eigenvalues[:]
		V1 = conj(f.root.complete.binding.eigenstates[:].transpose())
		E2 = f.root.complete.unbinding.eigenvalues[:]
		V2 = conj(f.root.complete.unbinding.eigenstates[:].transpose())
	finally:
		f.close()

	return E1, V1*dr, E2, V2*dr

def LoadBoundEigenstates(**args):	
	threshold = -0.5
	E1, V1, E2, V2 = LoadEigenstates(**args)

	thresholdIndex = max(where(E1<=threshold)[0])
	E1 = E1[:thresholdIndex]
	V1 = V1[:thresholdIndex,:]

	return E1, V1

def LoadRotatedBoundEigenstates(**args):
	"""
	Loads the eigenstates, and rotates them such that all
	the franck condon coefficients are real and positive
	"""
	args["silent"] = True
	args["configSilent"] = True

	#Get eigenstates
	E, V = LoadBoundEigenstates(**args)

	#Get initial state
	prop = SetupProblem(**args)
	LoadInitialState(prop, **args)

	#Expand initial state in eigenstaes
	initProj = dot(V, prop.psi.GetData()[:,0])
	initPhases = numpy.arctan2(initProj.imag, initProj.real)
	for i in range(len(E)):
		V[i,:] *= exp(-1.0j * initPhases[i])
	
	return E, V
	

def LoadContinuumEigenstates(**args):
	threshold = -0.5
	upperThreshold = 0.0

	E1, V1, E2, V2 = LoadEigenstates(**args)
	
	maxIndex1 = max(where(E1<=upperThreshold)[0]) + 1
	thresholdIndex = max(where(E1<=threshold)[0])
	E1 = E1[thresholdIndex:maxIndex1]
	V1 = V1[thresholdIndex:maxIndex1, :]

	maxIndex2 = max(where(E2<=upperThreshold)[0]) + 1
	E2 = E2[:maxIndex2]
	V2 = V2[:maxIndex2, :]

	return E1, V1, E2, V2


def LoadInitialState(prop, **args):
	"""
	Loads an initial state into prop.psi

	which initial condition is used is dependent on
	args["initType"] 
		if it is "fc", a FranckCondon distrib is used
		if it is "adk", the Niederhausen & Thumm distrib is used
	
	"""
	
	initType = "fc"
	if "initType" in args:
		initType = args["initType"]
	
	if initType == "fc":
		print "Using FranckCondon initial condition"
		LoadInitialStateFranckCondon(prop, GetInputFile(**args))

	elif initType == "adk":
		print "Using ADI initial condition"
		LoadInitialStateNiederhausen(prop, **args)

	else:
		raise Exception("Invalid initType (%s), specify 'fc' or 'adk'" % initType)


def LoadInitialStateFranckCondon(prop, inputfile):
	"""
	Loads a FranckCondon initial condition from a file
	(as created from SetupInitialCondition, and puts it
	into the wavefunction of prop
	"""
	f = tables.openFile(inputfile, "r")
	try:
		#Setup initial state
		prop.psi.GetData()[:,0] = f.root.initial_state[:]
		prop.psi.GetData()[:,1] = 0

		f.close()
	finally:
		f.close()
	del f


def LoadInitialStateNiederhausen(prop, **args):
	"""
	Loads the initial condition used by Niederhausen & Thumm, specified
	by the population of the vibrational levels.
	"""
	vibrationalPopulation = [0.17, 0.25, 0.23, 0.16, 0.09, 0.05, 0.03, 0.025, 0.006]
	E, V = LoadBoundEigenstates(**args)

	prop.psi.GetData()[:] = 0
	for i in range(len(vibrationalPopulation)):
		amplitude = sqrt(vibrationalPopulation[i])
		prop.psi.GetData()[:,0] += amplitude * V[i,:]

	prop.psi.Normalize()

