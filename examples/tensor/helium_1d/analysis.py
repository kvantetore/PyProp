#------------------------------------------------------------------------
#                Symmetrization related functions
#------------------------------------------------------------------------
def GetSymmetrizedWavefunction(psi):
	"""
	Symmetrizes and anti-symmetrizes the wavefunction with respect to
	particle exchange.

	Returns a tuple of the symmetrized and anti-symmetrized wavefunction 
	(symPsi, antiSymPsi)
	"""
	exchgPsi = psi.Copy()
	exchgPsi.GetData()[:] = exchgPsi.GetData().transpose()

	#create symmetrized wavefunction
	symPsi = psi.Copy()
	symPsi.GetData()[:] += exchgPsi.GetData()
	symPsi.GetData()[:] *= 0.5
	
	antiSymPsi = exchgPsi
	antiSymPsi.GetData()[:] -= psi.GetData()
	antiSymPsi.GetData()[:] *= 0.5

	return symPsi, antiSymPsi


def SymmetrizeWavefunction(psi, symmetrize):
	if symmetrize:
		symfactor = 1
	else:
		symfactor = -1

	exchgPsi = psi.GetData().transpose()

	#create symmetrized wavefunction
	exchgPsi.GetData()[:] *= symfactor
	psi.GetData()[:] += exchgPsi.GetData()
	psi.GetData()[:] *= 0.5


#------------------------------------------------------------------------
#                        Setting up 1D eigenstates
#------------------------------------------------------------------------
def SetupEigenstates1D(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	1D Hydrogen system.

	eigenvalues is a list of 1-d eigenvalue arrays. Each array corresponding
	to a
	"""
	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)
	bspl = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	eigenValues = []
	eigenVectors = []

	eigenvectorScaling = 1

	M = SetupPotentialMatrix(prop, potentialIndices)

	E, V = scipy.linalg.eig(a=M, b=S)

	idx = argsort(real(E))
	E = real(E[idx])
	eigenValues = E

	#Sort and normalize eigenvectors
	sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
	eigenVectors = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()

	#assure correct phase convention (first oscillation should start out real positive)
	for i, curE in enumerate(E):
		bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
		phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
		eigenVectors[:,i] *= exp(-1.0j * phase)

	return eigenValues, eigenVectors


def SaveEigenstates(**args):
	postfix = GetGridPostfix(**args)
	outputFile = "eigenstates/eigenstates_1D_%s.h5" % ("_".join(postfix))

	#Setup problem
	prop = SetupProblem(**args)

	#Setp eigenvalues and eigenvectors
	eigenValues, eigenVectors = SetupEigenstates1D(prop)

	#Save eigenvalues and eigenvectors to file
	if outputFile != None:
		f = tables.openFile(outputFile, "w")
		try:
			#Create main group
			eigGroup = f.createGroup(f.root, "Eig")

			#save config object
			eigGroup._v_attrs.configObject = prop.Config.cfgObj

			#save eigenvalues and eigenstates
			f.createArray(eigGroup, "eigenvalues", eigenValues)
			f.createArray(eigGroup, "eigenvectors", eigenVectors)

		finally:
			f.close()


def SetupPotentialMatrix(prop, whichPotentials):
	matrixSize = prop.psi.GetData().shape[0]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		basisPairs = potential.BasisPairs[0]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[i]

	return matrix


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupPotentialMatrix(prop, [overlap])
	return matrix


#------------------------------------------------------------------------
#                        Product State Analysis (examples)
#------------------------------------------------------------------------

def RunSingleIonizationStabilizationScan():
	#filenameTemplate = "stabilization_freq_3.0_scan_1s2p_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/stabilization_I_%i_kb20_dt_1e-02_T_12.6.h5"
	filenameTemplate = "output/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0/stabilization_I_%i_kb20_dt_1e-02_T_11.3.h5"
	outputPrefix = "stabilization_scan_freq_5.0_1s2s"
	intensity = r_[1:4]

	filenames = [filenameTemplate % i for i in intensity]
	absorb, totalIon, singleIon, doubleIon = RunSingleIonizationScan(filenames, outputPrefix)

	pylab.plot(intensity, totalIon, "-", label="Total Ion.")
	pylab.plot(intensity, singleIon, "--", label="Single Ion.")
	pylab.plot(intensity, doubleIon, ":", label="Double Ion.")
	xlabel("Intensity")
	title("Ionization Probability (1s2s) for w=5")
	ylim(0,1)
	pylab.legend(loc="lower right")
	pylab.savefig("%s.png" % outputPrefix)

	f = tables.openFile("%s.h5" % outputPrefix, "a")
	try:
		f.createArray(f.root, "intensity", intensity)
	finally:
		f.close()


def RunSingleIonizationEnergyDistributionScan():
	filenameTemplate = "raymond/stabilization_freq_5_scan_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/stabilization_I_%i_kb20_dt_1e-02.h5"
	outputPrefix = "dpde_scan"
	intensity = r_[1:37:1]

	filenames = [filenameTemplate % i for i in intensity]
	E, dpde = RunGetSingleIonizationEnergyDistribution(filenames)

	f = tables.openFile("%s.h5" % outputPrefix, "w")
	try:
		f.createArray(f.root, "intensity", intensity)
		f.createArray(f.root, "energy", E)
		f.createArray(f.root, "dpde", array(dpde))
	finally:
		f.close()

def RunDoubleIonizationAngularDistributionScan(dE=None, dTheta=None):
	filenameTemplate = "raymond/stabilization_freq_5_scan_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/extra_cycles_propagation/stabilization_I_%i_kb20_dt_1e-02.h5"
	outputPrefix = "stabilization_freq5_angular_distrib_scan"
	intensity = r_[1:37:1]

	filenames = [filenameTemplate % i for i in intensity]
	f = tables.openFile("output/%s.h5" % outputPrefix, "w")
	try:
		for i, filename in enumerate(filenames):
			energy, theta, (distribution,) = RunGetDoubleIonizationAngularDistribution([filename], dE, dTheta)
			f.createArray(f.root, "distrib_%i" % i, distribution)
		f.createArray(f.root, "intensity", intensity)
		f.createArray(f.root, "theta", theta)
		f.createArray(f.root, "energy", energy)
		
	finally:
		f.close()


#------------------------------------------------------------------------
#                        Product State Analysis (implementation)
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#         Run a series of calculations and save the result
#------------------------------------------------------------------------

def RunFolderUpdateAnalysis(folder, noOverwrite=False):
	"""
	Checks every folder and subfolder, and updates every h5 file containing a /wavefunction
	- single ionization dpde, 
	- double ionization dpde,
	- double ionization dpdedomega

	"""
	
	#walk all h5 files files 
	for dirpath, subfolders, files in os.walk(folder):
		for file in filter(lambda s: s.endswith(".h5"), files):
			filename = os.path.join(dirpath, file)
			RunFileUpdateAnalysis(filename, noOverwrite)
	


def RunFileUpdateAnalysis(filename, noOverwrite=False):
	"""
	Updates 'filename' (if it contains a /wavefunction):
	- single ionization dpde, 
	- double ionization dpde,
	- double ionization dpdedomega
	"""

	#check file contains wavefunction
	f = tables.openFile(filename)
	try:
		if not "/wavefunction" in f:
			print "Found no /wavefunction, so I skip this file: %s" % file
			return
		if "/dpdomega_double" in f and noOverwrite:
			print "Skipping this file: %s" % filename
			return
	finally:
		f.close()

	print "Updating analysis for file %s" % (filename, )
	print "    ionization prob"
	pyprop.Redirect.Enable(silent=True)
	(symProb,), (antiProb,) = RunGetSingleIonizationProbability([filename]) 
	pyprop.Redirect.Disable()
	print "    dP/dOmega (double)"
	pyprop.Redirect.Enable(silent=True)
	e1, th1, (dp1,) = RunGetDoubleIonizationAngularDistribution([filename])
	pyprop.Redirect.Disable()
	print "    dP/dE (double)"
	pyprop.Redirect.Enable(silent=True)
	e2, (dp2,) = RunGetDoubleIonizationEnergyDistribution([filename])
	pyprop.Redirect.Disable()
	print "    dP/dE (single)"
	pyprop.Redirect.Enable(silent=True)
	e3, (dp3,) = RunGetSingleIonizationEnergyDistribution([filename])
	pyprop.Redirect.Disable()
	print "    dP/dOmega (single)"
	pyprop.Redirect.Enable(silent=True)
	e4, th4, (dp4,) = RunGetSingleIonizationAngularDistribution([filename])
	pyprop.Redirect.Disable()

	print "    saving results"
	f = tables.openFile(filename, "a")
	try:
		if "dpdomega" in f.root:
			f.removeNode(f.root, "dpdomega", recursive=True)
		if "dpdomega_double" in f.root:
			f.removeNode(f.root, "dpdomega_double", recursive=True)
		f.createArray(f.root, "dpdomega_double", dp1)
		f.root.dpdomega_double._v_attrs.energy = e1
		f.root.dpdomega_double._v_attrs.theta = th1

		if "dpde_double" in f.root:
			f.removeNode(f.root, "dpde_double", recursive=True)
		f.createArray(f.root, "dpde_double", dp2)
		f.root.dpde_double._v_attrs.energy = e2

		if "dpde_single" in f.root:
			f.removeNode(f.root, "dpde_single", recursive=True)
		f.createArray(f.root, "dpde_single", dp3)
		f.root.dpde_single._v_attrs.energy = e3

		if "dpdomega_single" in f.root:
			f.removeNode(f.root, "dpdomega_single", recursive=True)
		f.createArray(f.root, "dpdomega_single", dp4)
		f.root.dpdomega_single._v_attrs.energy = e4
		f.root.dpdomega_single._v_attrs.theta = th4

		attrs = f.root.wavefunction._v_attrs
		attrs.Absorbed = symProb[0]
		attrs.Ionization = symProb[1]
		attrs.SingleIonization = symProb[2]
		attrs.DoubleIonization = symProb[3]
		attrs.AntiAbsorbed = antiProb[0]
		attrs.AntiIonization = antiProb[1]
		attrs.AntiSingleIonization = antiProb[2]
		attrs.AntiDoubleIonization = antiProb[3]


	finally:
		f.close()


def RunGetProductStatePopulations(fileList, scanParameter, outputFile, removeBoundStates=True):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	maxE = 30.

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: 0.0 < E
	isIonizedCutoff = lambda E: 0.0 < E <= maxE
	isBound = lambda E: E <= 0.
	doubleIonEnergies, doubleIonStates = GetFilteredSingleParticleStates("he+", isIonizedCutoff, config=conf)
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)


	f = tables.openFile(outputFile, "w")
	try:
		f.createArray(f.root, "scanParameter", scanParameter)
		f.createArray(f.root, "filenames", fileList)
		f.createVLArray(f.root, "singleBoundEnergies", atom=tables.ObjectAtom()).append(singleBoundEnergies)
		f.createVLArray(f.root, "singleIonEnergies", atom=tables.ObjectAtom()).append(singleIonEnergies)
		f.createVLArray(f.root, "doubleIonEnergies", atom=tables.ObjectAtom()).append(doubleIonEnergies)

		f.createVLArray(f.root, "doubleIonStates", atom=tables.ObjectAtom()).append(doubleIonStates)
		f.createVLArray(f.root, "singleIonStates", atom=tables.ObjectAtom()).append(singleIonStates)
		f.createVLArray(f.root, "singleBoundStates", atom=tables.ObjectAtom()).append(singleBoundStates)

		#Calculate Energy Distribution (dP/dE1 dE2)
		for i, filename in enumerate(fileList):
			psi = pyprop.CreateWavefunctionFromFile(filename)
			sym, anti = GetSymmetrizedWavefunction(psi)
		
			#get absorbed prob
			absorbedProbability = 1.0 - real(psi.InnerProduct(psi))
			
			#remove boundstate projection
			if boundStates != None:
				RemoveBoundStateProjection(psi, boundStates)
			ionizationProbability = real(psi.InnerProduct(psi))
		
			#Get single ionization populations
			singleIonPop = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)
			singleIonization = sum([sum([p for i1, i2, p in pop]) for l1, l2, pop in singleIonPop ])
			#Get double ionization populations
			doubleIonPop = GetPopulationProductStates(psi, doubleIonStates, doubleIonStates)

			#save
			grp = f.createGroup(f.root, "parameter_%i" % i)
			grp._v_attrs.AbsorbedProbability = absorbedProbability
			grp._v_attrs.TotalIonization = ionizationProbability
			grp._v_attrs.SingleIonization = singleIonization
			grp._v_attrs.SymmetrizedProbability = real(sym.InnerProduct(sym))
			grp._v_attrs.AntiSymmetrizedProbability = real(anti.InnerProduct(anti))
			f.createVLArray(grp, "singleIonPop", atom=tables.ObjectAtom()).append(singleIonPop)
			f.createVLArray(grp, "doubleIonPop", atom=tables.ObjectAtom()).append(doubleIonPop)
			del grp
		
	finally:
		f.close()


def FromFileCalculateSingleIonizationDPDE(filename, scanIndex=-1, maxE=15, dE=0.1):
	E = r_[0:maxE:dE]

	f = tables.openFile(filename, "r")
	try:
		singleBoundEnergies = f.root.singleBoundEnergies[0]
		singleIonEnergies = f.root.singleIonEnergies[0]
		params = f.root.scanParameter[:]

		def calculateSingleIndex(index):
			grp = f.getNode("/parameter_%i" % index)
			singleIonPop = grp.singleIonPop[0]

			#Create an array for dpde
			dpde = zeros(len(E), dtype=double)
			
			#for every l-pair (lBound, lIon) we have a set of (iBound, iIon) states
			#In order to create an approx to dp/de_ion, we interpolate for each iBound,
			#and add coherently to the dpde array
			for lBound, lIon, lPop in singleIonPop:
				#number of states in this l-shell
				nBound = len(singleBoundEnergies[lBound])
				nIon = len(singleIonEnergies[lIon])
			
				#iterate over all bound states in this l-shell
				pop = array([d[2] for d in lPop]).reshape(nBound, nIon)
				for iBound in range(nBound):
					#interpolate over ionized singleIonPop
					curPop = pop[iBound, :-1] / diff(singleIonEnergies[lIon])
					dpde += interp(E, singleIonEnergies[lIon][:-1], curPop)
			
			return dpde


		if scanIndex != -1:
			return params[scanIndex], E, calculateSingleIndex(scanIndex)
		else:
			return params, E, map(calculateSingleIndex, r_[:len(params)])

	finally:
		f.close()

def FromFileCalculateDoubleIonizationDPDE(filename, scanIndex=-1, maxE=15, dE=0.1):
	E = r_[0:maxE:dE]

	f = tables.openFile(filename, "r")
	try:
		doubleIonEnergies = f.root.doubleIonEnergies[0]
		params = f.root.scanParameter[:]

		def calculateSingleIndex(index):
			grp = f.getNode("/parameter_%i" % index)
			doubleIonPop = grp.doubleIonPop[0]

			dpde = zeros((len(E), len(E)), dtype=double)
			for l1, l2, lPop in doubleIonPop:
				#number of states in this l-shell
				n1 = len(doubleIonEnergies[l1])
				n2 = len(doubleIonEnergies[l2])
			
				#scale states with 1/dE_1 dE_2
				pop = array([d[2] for d in lPop]).reshape(n1, n2)
				E1 = doubleIonEnergies[l1]
				E2 = doubleIonEnergies[l2]
				pop[:-1,:-1] /= outer(diff(E1), diff(E2))
				
				#2d interpolation over all states in this shell
				interpolator = scipy.interpolate.RectBivariateSpline(E1[:-1], E2[:-1], pop[:-1, :-1], kx=1, ky=1)
				dpde += interpolator(E, E)
	
			return dpde


		if scanIndex != -1:
			return params[scanIndex], E, calculateSingleIndex(scanIndex)
		else:
			return params, E, map(calculateSingleIndex, r_[:len(params)])

	finally:
		f.close()



#------------------------------------------------------------------------
#        Load Single Particle States
#------------------------------------------------------------------------


def GetSingleStatesFile(**args):
	locations = ["eigenstates"]
	gridPostfix = "_".join(GetGridPostfix(**args))
	singleStatesFile = filter(os.path.exists, ["%s/eigenstates_1D_%s.h5" % (loc, gridPostfix) for loc in locations])
	if len(singleStatesFile) == 0:
		raise Exception("Could not find single-particle data files!")
	return singleStatesFile[0]


def LoadSingleParticleStates(singleStatesFile):
	f = tables.openFile(singleStatesFile, "r")

	#Create BSpline object
	cfgSection = pyprop.Section("BsplineRepresentation", f.root.Eig._v_attrs.configObject)
	bspl = pyprop.BSPLINE()
	bspl.ApplyConfigSection(cfgSection)

	#setup grid to check for correct phase convention
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	try:
		node = f.getNode("/Eig")
		E = node.eigenvalues[:]
		V = node.eigenvectors[:]

		#assure correct phase convention (first oscillation should start out real positive)
		#for i, curE in enumerate(E):
		#	bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
		#	phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
		#	V[:,i] *= exp(-1.0j * phase)

	finally:
		f.close()

	return E, V


def GetFilteredSingleParticleStates(stateFilter, **args):
	"""
	Returns the single particle states and energies corresponding to the 1D model and 
	an energy filter specified by stateFilter
	"""

	#load single particle states
	singleStatesFile = GetSingleStatesFile(**args)
	singleEnergies, singleStates = LoadSingleParticleStates(singleStatesFile)
	
	#filter states 
	filteredStates = [singleStates[:,i] for (i,E) in enumerate(singleEnergies) if stateFilter(E)]
	filteredEnergies = filter(stateFilter, singleEnergies)

	return filteredEnergies, filteredStates	


def SetupRadialCoulombStatesEnergyNormalized(psi, Z, Emax, dE, lmax):
	E = r_[dE:Emax:dE]
	k = sqrt(E*2)

	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)
	rcount = psi.GetRepresentation().GetRepresentation(1).GetFullShape()[0]
	
	#Setup Radial Waves
	states = []
	for l in range(lmax+1):
		V = zeros((rcount, len(k)), dtype=complex)
		for i, curk in enumerate(k):
			coeff = GetRadialCoulombWaveBSplines(Z, l, curk, bspline)
			V[:,i] = sqrt(2*dE/pi/curk) * coeff
		states.append(V)
	
	energies = [E]*(lmax+1)

	return energies, states


def GetSingleParticleCoulombStates(Z, dk, mink, maxk, lmax, radialRepr):
	"""
	Gets coulomb wave functions for every k between mink and maxk (in dk steps), 
	for every l up to (and including) lmax evaluated in bsplines.

	The structure returned is similar to that of LoadSingleParticleStates
	"""
	bspl = radialRepr.GetBSplineObject()
	k = r_[mink:maxk:dk]
	rcount = radialRepr.GetFullShape()[0]
	
	states = []
	for l in range(lmax+1):
		V = zeros((rcount, len(k)), dtype=complex)
		for i, curk in enumerate(k):
			coeff = GetRadialCoulombWaveBSplines(Z, l, curk, bspl)
			V[:,i] = coeff
		states.append(V)

	return [k]*(lmax+1), states


def GetRadialCoulombWaveBSplines(Z, l, k, bsplineObj):
	#Get the Coulomb function in grid space
	r = bsplineObj.GetQuadratureGridGlobal()
	wav = zeros(len(r), dtype=double)
	SetRadialCoulombWave(Z, l, k, r, wav)
	cplxWav = array(wav, dtype=complex)

	#get bspline coeffs
	coeff = zeros(bsplineObj.NumberOfBSplines, dtype=complex)
	bsplineObj.ExpandFunctionInBSplines(cplxWav, coeff)

	return coeff


#------------------------------------------------------------------------
#        Calculations on general product state combinations
#------------------------------------------------------------------------
def CalculatePopulationProductStates(V1, V2, psiData):
	numBsplines = psiData.shape[0]
	numStates1 = shape(V1)[0]
	numStates2 = shape(V2)[0]
	print "numStates1 = %s, numStates2 = %s" % (numStates1, numStates2)
	tempData = zeros((numBsplines, numStates2))
	populations = zeros((numStates1, numStates2))	

	#Project on V2 states
	for i, v2 in enumerate(V2):
		tempData[:,i] = dot(conj(v2), psiData[:])

	#Project on V1 states
	for j, v1 in enumerate(V1):
		populations[j,:] = dot(conj(v1), tempData[:])
	
	populations *= conj(populations)

	popList = []
	for i1 in range(numStates1):
		for i2 in range(numStates2):
			popList += [[i1, i2, 2 * real(populations[i1,i2])]]

	return popList
		

def GetPopulationProductStates(psi, singleStates1, singleStates2):
	"""
	Calculates the population of psi in a set of single electron product states

	P_i = 2|< SingleState1_i(1), SingleState2_j(2) | psi(1,2) >|^2

	singleStates 1 and 2 are lists of 1D eigenstates created by SetupEigenstates
	
	the projection is carried out for every combination of singlestate1 and singlestate2
	"""

	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	repr = psi.GetRepresentation()
	repr.MultiplyIntegrationWeights(tempPsi)
	data = tempPsi.GetData()

	#Get the population for every combination of v1 and v2
	popList = CalculatePopulationProductStates(singleStates1, singleStates2, data)

	return popList


def RemoveProductStatesProjection(psi, singleStates1, singleStates2):
	"""
	Makes psi orthogonal to a set of single electron product states

	psi =  (I - sum_{i,j} |j(2), i(1)> <i(1), j(2)|) | psi(1,2) >

	where i and j are single particles states.
	
	singleStates1 and 2 are lists of angular momentum states, containing an array 
	of radial states for the given angular momentum number such as generated
	by SetupRadialEigenstates in the Helium SAE example
	
	the projection is carried out for every combination of singlestate1 and singlestate2i
	"""

	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	repr = psi.GetRepresentation()
	repr.MultiplyIntegrationWeights(tempPsi)
	distr = psi.GetRepresentation().GetDistributedModel()

	data = tempPsi.GetData()
	population = []

	for l1, V1 in enumerate(singleStates1):
		print "%i/%i" % (l1, len(singleStates1))

		for l2, V2 in enumerate(singleStates2):

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(psi, lfilter), dtype=int32)
			if len(angularIndices) == 0:
				continue
		
			#Remove projection for every combination of v1 and v2
			projV = RemoveProjectionRadialProductStates(l1, V1, l2, V2, data, angularIndices, psi.GetData())


	return population


#------------------------------------------------------------------------
#                        Product state model ionization probabilities
#------------------------------------------------------------------------
def GetProductStateModelIonizationProbabilities(modelFirstElectron, modelSecondElectron):
	"""
	Calculate two-electron product state model ionization probabilities
	from SAE ionization data.
	"""
#	if modelFirstElectron == modelSecondElectron:
#		singleIon, 
