
def SetupConfig(**args):
	#Decide which config file to use
	configFile = "config.ini"
	if "config" in args:
		configFile = args["config"]

	#Load the config file
	conf = pyprop.Load(configFile)

	if "dt" in args:
		dt = args["dt"]
		conf.Propagation.timestep = dt

	#Modify the config
	if "imtime" in args:
		imtime = args["imtime"]
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

	if "molecule" in args:
		molecule = args["molecule"].lower()
		if molecule == "h2+":
			conf.RadialPropagator.mass = ReducedMass(mass_proton_au, mass_proton_au)
			conf.ElectronicCoupling.static_dipole_moment = 0
		elif molecule == "d2+":
			conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_deuteron_au)
			conf.ElectronicCoupling.static_dipole_moment = 0
		elif molecule == "hd+":
			conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_proton_au)
			conf.ElectronicCoupling.static_dipole_moment = 1. / 3.
		else:
			raise Exception("Unknown molecule '%s'" % molecule)
	
	if "staticDipoleMoment" in args:
		staticDipoleMoment = args["staticDipoleMoment"]
		conf.ElectronicCoupling.static_dipole_moment = staticDipoleMoment

	if 'radialScaling' in args:
		radialScaling = args["radialScaling"]
		rank0 = conf.RadialRepresentation.rank0
		radialRange = rank0[1] - rank0[0]
		radialSize = rank0[2]
	
		rank0[1] = radialRange / radialScaling + rank0[0]
		rank0[2] = radialSize / radialScaling

	if 'pulseDelay' in args:
		pulseDelay = args["pulseDelay"]
		conf.ElectronicCoupling.delay = pulseDelay

	if "pulsePhase" in args:
		pulsePhase = args["pulsePhase"]
		conf.ElectronicCoupling.phase = pulsePhase

	if "pulseDuration" in args:
		pulseDuration = args["pulseDuration"]
		conf.ElectronicCoupling.duration = pulseDuration
	
	if "pulseIntensity" in args:
		pulseIntensity = args["pulseIntensity"]
		conf.ElectronicCoupling.intensity = pulseIntensity

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration

	if "silent" in args:
		silent = args["silent"]
		conf.Propagation.silent = silent

	if "potentialSlope" in args:
		conf.VibrationalPotential.potential_slope=args["potentialSlope"]
		
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
	
	for t in prop.Advance(10):
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
	H = GetHamiltonMatrix(**args)
	E, V = eig(H)

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
	#SetupEigenstates(**args)
	SetupFullSpectrum(**args)
	SetupInitialState(**args)


