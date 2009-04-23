#------------------------------------------------------------------------
#                Spherical Harmonics (partial waves) analysis
#------------------------------------------------------------------------

def RunSphericalHarmonicDistribution(wavefunctionFile):
	if isinstance(wavefunctionFile, str):
		rightPsi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
	else:
		rightPsi = wavefunctionFile
	leftPsi = rightPsi.Copy()

	#data = <left | right> = sum conj(left) * S * right
	repr = rightPsi.GetRepresentation()
	angRepr = repr.GetRepresentation(0)
	repr.MultiplyIntegrationWeights(rightPsi)
	data = real(conj(leftPsi.GetData()) * rightPsi.GetData())

	data = sum(sum(data, axis=2), axis=1)
	return data

	for i in range(pyprop.ProcCount):
		if i == pyprop.ProcId:
			print "Proc %i" % pyprop.ProcId
			for dataIndex, angIndex in enumerate(repr.GetLocalGrid(0)):
				coupledIndex = angRepr.Range.GetCoupledIndex(int(angIndex))
				print "%s = %s" % (coupledIndex, data[dataIndex])
		pyprop.pypar.barrier()


def GetLocalCoupledSphericalHarmonicIndices(psi, coupledIndexFilter):
	"""
	Returns the processor local indices which corresponds to a filter on
	l1, l2, L, M
	"""
	angularRank = 0

	#Get info about angular representation
	repr = psi.GetRepresentation().GetRepresentation(angularRank)
	distr = psi.GetRepresentation().GetDistributedModel()
	nL = repr.GetFullShape()[0]

	#Filter global indices which are on this processor
	localStart = distr.GetLocalStartIndex(int(nL), int(angularRank))
	localEnd = localStart + nL
	isLocal = lambda i: localStart <= i < localEnd
	globalIndices = filter(isLocal, r_[:nL])

	#Filter indices where coupled index is corresponding to supplied filter
	curFilter = lambda i: coupledIndexFilter(repr.Range.GetCoupledIndex(i))
	globalFilteredIndices = filter(curFilter, globalIndices)

	#map global to local indices
	globalToLocal = lambda i: i - localStart
	localFilteredIndices = map(globalToLocal, globalFilteredIndices)

	return localFilteredIndices


def GetCoupledIndexList(psi):
	"""
	Returns a list of the coupled indices in psi
	"""
	angularRank = 0
	repr = psi.GetRepresentation().GetRepresentation(angularRank)
	distr = psi.GetRepresentation().GetDistributedModel()
	nL = repr.GetFullShape()[0]
	coupledIndexList = map(repr.Range.GetCoupledIndex, range(nL))
	return coupledIndexList


def GetSymmetrizationIndexPairs(psi):
	"""
	Returns index pairs corresponding to symmetrization of the basis:

	L = L' L = L' and l1 = l2' and l2 = ll'

	returns a list of indices (i1, i2) 
	such that psi[i1, r1, r2] -> psi[i2, r2, r1] is particle exchange
	psi(1, 2) -> psi(2, 1)
	"""
	angularRank = 0

	coupledIndexList = GetCoupledIndexList(psi)
	nL = len(coupledIndexList)
	
	#construct all possible pairs
	coupledPairs = [(i1, i2) for i1 in range(nL) for i2 in range(nL)]
	#filter out symmetry pairs
	symmetryFilter = lambda l, r: l.L==r.L and l.M==r.M and l.l1==r.l2 and l.l2==r.l1
	#symmetryFilter = lambda l, r: l.L==r.L and l.M==r.M and l.l1==r.l1 and l.l2==r.l2
	symmetryIndexFilter = lambda i: symmetryFilter(coupledIndexList[i[0]], coupledIndexList[i[1]])
	symmetryPairs = filter(symmetryIndexFilter, coupledPairs)

	return symmetryPairs

def GetSymmetrizedWavefunction(psi):
	"""
	Symmetrizes and anti-symmetrizes the wavefunction with respect to
	particle exchange.

	Returns a tuple of the symmetrized and anti-symmetrized wavefunction 
	(symPsi, antiSymPsi)
	"""
	AssertSingleProc()

	sym = GetSymmetrizationIndexPairs(psi)
	exchgPsi = GetWavefunctionParticleExchange(psi, sym)

	#create symmetrized wavefunction
	symPsi = psi.Copy()
	symPsi.GetData()[:] += exchgPsi.GetData()
	symPsi.GetData()[:] *= 0.5
	
	antiSymPsi = exchgPsi
	antiSymPsi.GetData()[:] -= psi.GetData()
	antiSymPsi.GetData()[:] *= 0.5

	return symPsi, antiSymPsi

def SymmetrizeWavefunction(psi, symmetrize):
	AssertSingleProc()

	if symmetrize:
		symfactor = 1
	else:
		symfactor = -1

	sym = GetSymmetrizationIndexPairs(psi)
	exchgPsi = GetWavefunctionParticleExchange(psi, sym)

	#create symmetrized wavefunction
	exchgPsi.GetData()[:] *= symfactor
	psi.GetData()[:] += exchgPsi.GetData()
	psi.GetData()[:] *= 0.5

#------------------------------------------------------------------------
#                        Product State Analysis (examples)
#------------------------------------------------------------------------

def RunSingleIonizationStabilizationScan():
	filenameTemplate = "stabilization_freq_3.0_scan_1s2p_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/stabilization_I_%i_kb20_dt_1e-02_T_12.6.h5"
	outputPrefix = "stabilization_scan_freq_3.0_1s2p"
	intensity = r_[1:16]

	filenames = [filenameTemplate % i for i in intensity]
	absorb, totalIon, singleIon, doubleIon = RunSingleIonizationScan(filenames, outputPrefix)

	pylab.plot(intensity, totalIon, "-", label="Total Ion.")
	pylab.plot(intensity, singleIon, "--", label="Single Ion.")
	pylab.plot(intensity, doubleIon, ":", label="Double Ion.")
	xlabel("Intensity")
	title("Ionization Probability (1s2p) for w=3")
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

def RunDoubleIonizationAngularDistributionScan():
	filenameTemplate = "raymond/stabilization_freq_5_scan_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/extra_cycles_propagation/stabilization_I_%i_kb20_dt_1e-02.h5"
	outputPrefix = "stabilization_freq5_angular_distrib_scan"
	intensity = r_[1:37:1]

	filenames = [filenameTemplate % i for i in intensity]
	f = tables.openFile("output/%s.h5" % outputPrefix, "w")
	try:
		for i, filename in enumerate(filenames):
			energy, theta, (distribution,) = RunGetDoubleIonizationAngularDistribution([filename])
			f.createArray(f.root, "distrib_%i" % i, distribution)
		f.createArray(f.root, "intensity", intensity)
		f.createArray(f.root, "theta", theta)
		f.createArray(f.root, "energy", energy)
		
	finally:
		f.close()


#------------------------------------------------------------------------
#                        Product State Analysis (implementation)
#------------------------------------------------------------------------

def GetFilteredSingleParticleStates(model, stateFilter, **args):
	"""
	Returns the single particle states and energies corresponding to a SAE model and 
	an energy filter specified by stateFilter
	"""

	#load single particle states
	singleStatesFile = GetSingleStatesFile(model=model, **args)
	lList, singleEnergies, singleStates = LoadSingleParticleStates(singleStatesFile)
	
	#filter states 
	getStates = lambda E, V: transpose(array([ V[:, i] for i in range(V.shape[1]) if stateFilter(E[i]) ]))
	getEnergies = lambda E: filter(stateFilter, E)
	filteredStates = map(getStates, singleEnergies, singleStates)
	filteredEnergies = map(getEnergies, singleEnergies)

	return filteredEnergies, filteredStates	


def GetAssociatedLegendrePoly(lmax, theta):
	leg = []
	for l in range(lmax+1):
		for m in range(-l, l+1):
			leg.append(scipy.special.sph_harm(m, l, 0, theta))
	return array(leg)


def GetSingleParticleCoulombStates(Z, dk, mink, maxk, lmax, radialRepr):
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


def RunGetSingleIonizationProbability(fileList, removeBoundStates=True):
	"""
	Calculates total and single ionization probability
	for a list of wavefunction files
	"""

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: E > 0.0
	isBound = lambda E: not isIonized(E)
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)

	def getIonProb(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		#SymmetrizeWavefunction(psi, True)
		#return GetSingleIonizationProbability(psi, boundStates, singleBoundStates, singleIonStates)
		sym, anti = GetSymmetrizedWavefunction(psi)
		symProb = GetSingleIonizationProbability(sym, boundStates, singleBoundStates, singleIonStates)
		#antiProb = GetSingleIonizationProbability(anti, boundStates, singleBoundStates, singleIonStates)
		#return (symProb, antiProb)
		return (symProb,)

	return zip(*map(getIonProb, fileList))


def RunGetSingleIonizationEnergyDistribution(fileList, removeBoundStates=True):
	"""
	Calculates the energy distribution (dP/dE) of the 
	single ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	maxE = 15.
	dE = 0.01

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None


	#Get single particle states
	isIonized = lambda E: 0.0 <= E <= maxE
	isBound = lambda E: E < 0.0
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)

	#Calculate Energy Distribution (dP/dE)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStates, singleIonStates, singleIonEnergies, dE, maxE)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


def RunGetDoubleIonizationEnergyDistribution(fileList, removeBoundStates=True):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	maxE = 15.
	dE = 0.1

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: 0.0 < E <= maxE
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)

	#Calculate Energy Distribution (dP/dE1 dE2)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetDoubleIonizationEnergyDistribution(psi, boundStates, singleIonStates, singleIonEnergies, dE, maxE)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


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



def RunGetDoubleIonizationAngularDistribution(fileList, removeBoundStates=True, removeSingleIonStates=False):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""

	Z = 2
	#dk = 0.01
	#mink = 0.1
	#maxk = 5
	lmax = 5
	#theta = array([0., pi/2, pi], dtype=double)
	theta = linspace(0, pi, 20)

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: 0.0 < E
	isFilteredIonized = lambda E: 0.0 < E < 15
	isBound = lambda E: E <= 0.
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)
	doubleIonEnergies, doubleIonStates = GetFilteredSingleParticleStates("he+", isFilteredIonized, config=conf)

	dE = maximum(min(diff(singleIonEnergies[0])), 0.05)
	interpE = r_[dE:15:dE]
	interpK = sqrt(interpE * 2)


	#Get Coulomb states
	#psi = pyprop.CreateWavefunctionFromFile(fileList[0])
	#repr = psi.GetRepresentation().GetRepresentation(1)
	#singleK, coulombStates = GetSingleParticleCoulombStates(Z=Z, dk=dk, mink=mink, maxk=maxk, lmax=lmax, radialRepr=repr)
	#k = singleK[0]
	#del psi

	#Get spherical harmonics
	assocLegendre = GetAssociatedLegendrePoly(lmax, theta)

	angDistr = []
	for i, filename in enumerate(fileList):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		sym, anti = GetSymmetrizedWavefunction(psi)
		psi = sym
		print "norm = %s" % psi.GetNorm()**2
	
		#get absorbed prob
		absorbedProbability = 1.0 - real(psi.InnerProduct(psi))
		
		#remove boundstate projection
		if boundStates != None:
			RemoveBoundStateProjection(psi, boundStates)
		ionizationProbability = real(psi.InnerProduct(psi))

		#remove single ionization
		if removeSingleIonStates:
			RemoveProductStatesProjection(psi, singleBoundStates, singleIonStates)
			RemoveProductStatesProjection(psi, singleIonStates, singleBoundStates)

		doubleIonProb = real(psi.InnerProduct(psi))
		print "Double Ionization Probability = %s" % (doubleIonProb,)

		#calculate angular distr for double ionized psi
		angDistr.append(GetDoubleAngularDistribution(psi, Z, interpE, doubleIonEnergies, doubleIonStates, assocLegendre, theta))

	return interpE, theta, angDistr


def	GetSingleIonizationProbability(psi, boundStates, singleBoundStates, singleIonStates):
	"""
	Calculates the single ionization probability by first projecting 
	away doubly bound states, then projecting on products of single 
	particle states, where one particle is bound, and the other is ionized

	returns 
		- absorbed probabilty: norm when psi enters routine
		- ionization probability: norm when all doubly bound states are projected away
		- single ionization probability: projection on single particle states
	"""

	#get absorbed prob
	absorbedProbability = 1.0 - real(psi.InnerProduct(psi))

	#remove boundstate projection
	if boundStates != None:
		RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	##calculate populations in product states containing bound he+ states
	populations = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)

	#remove single ion projection
	#RemoveProductStatesProjection(psi, singleBoundStates, singleIonStates)
	#RemoveProductStatesProjection(psi, singleIonStates, singleBoundStates)
	#singleIonizationProbability2 = ionizationProbability - real(psi.InnerProduct(psi))
	singleIonizationProbability2 = 0

	#Calculate single ionization probability
	lpop = [sum([p for i1, i2, p in pop]) for l1, l2, pop in populations]
	singleIonizationProbability = sum(lpop)
	
	print "Absorbed Probability     = %s" % (absorbedProbability)
	print "Ioniziation Probability  = %s" % (ionizationProbability)
	print "Single Ionization Prob.  = %s, %s" % (singleIonizationProbability, singleIonizationProbability2)
	print "Double Ionization Prob.  = %s" % (ionizationProbability - singleIonizationProbability)
	print "Single Ionization ratio  = %s" % (singleIonizationProbability/ionizationProbability)

	#return absorbedProbability, ionizationProbability, singleIonizationProbability
	return singleIonizationProbability2


def GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStates, singleIonStates, singleIonEnergies, dE, maxE):
	"""
	Calculates dP/dE for the single ionized part of the wavefunction, by first
	projecting away doubly bound states, and then projecting on single particle
	product states, and binning the probability by energy.
	"""

	#get absorbed prob
	absorbedProbability = 1.0 - real(psi.InnerProduct(psi))

	#remove boundstate projection
	if boundStates != None:
		RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	#calculate populations in product states containing bound he+ states
	populations = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)

	#Create an array for dpde
	E = r_[0:maxE:dE]
	dpde = zeros(len(E), dtype=double)

	#for every l-pair (lBound, lIon) we have a set of (iBound, iIon) states
	#In order to create an approx to dp/de_ion, we interpolate for each iBound,
	#and add coherently to the dpde array
	for lBound, lIon, lPop in populations:
		#number of states in this l-shell
		nBound = singleBoundStates[lBound].shape[1]
		nIon = singleIonStates[lIon].shape[1]

		#iterate over all bound states in this l-shell
		pop = array([d[2] for d in lPop]).reshape(nBound, nIon)
		for iBound in range(nBound):
			#interpolate over ionized populations
			curPop = pop[iBound, :-1] / diff(singleIonEnergies[lIon])
			dpde += interp(E, singleIonEnergies[lIon][:-1], curPop)
		
	return E, dpde


def GetDoubleIonizationEnergyDistribution(psi, boundStates, singleIonStates, singleEnergies, dE, maxE):
	"""
	Calculates double differential d^2P/(dE_1 dE_2) by 
	1) projecting on a set of product of single particle ionized states.
	2) binning the probability by energy

	"""
	#get absorbed prob
	absorbedProbability = 1.0 - real(psi.InnerProduct(psi))

	#remove boundstate projection
	if boundStates != None:
		RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	populations = GetPopulationProductStates(psi, singleIonStates, singleIonStates)


	E = r_[0:maxE:dE]
	dpde = zeros((len(E), len(E)), dtype=double)

	#for every l-pair (l1, l2) we have a set of (i1, i2) states
	#In order to create an approx to dp/de_1 de_2, we make a 2d 
	#interpolation for each l-shell
	#and add coherently to the dpde array
	for l1, l2, lPop in populations:
		#number of states in this l-shell
		n1 = singleIonStates[l1].shape[1]
		n2 = singleIonStates[l2].shape[1]

		#scale states with 1/dE_1 dE_2
		pop = array([d[2] for d in lPop]).reshape(n1, n2)
		E1 = singleEnergies[l1]
		E2 = singleEnergies[l2]
		meshE1, meshE2 = meshgrid(E1[:-1], E2[:-1])
		pop[:-1,:-1] /= outer(diff(E1), diff(E2))
		
		#2d interpolation over all states in this shell
		interpolator = scipy.interpolate.RectBivariateSpline(E1[:-1], E2[:-1], pop[:-1, :-1], kx=1, ky=1)
		dpde += interpolator(E, E)
	
	return E, dpde



def GetSingleStatesFile(**args):
	radialGridPostfix = "_".join(GetRadialGridPostfix(**args))
	model = args.get("model", "he+")
	singleStatesFile = "output/singleelectron/eigenstates_sae_model_%s_%s.h5" % (model, radialGridPostfix)
	return singleStatesFile


def LoadSingleParticleStates(singleStatesFile):
	f = tables.openFile(singleStatesFile, "r")

	#Create BSpline object
	cfgSection = pyprop.Section("RadialRepresentation", f.root.RadialEig._v_attrs.configObject)
	bspl = pyprop.BSPLINE()
	bspl.ApplyConfigSection(cfgSection)

	#setup grid to check for correct phase convention
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	eigenvalues = []
	eigenvectors = []
	try:
		lList = f.root.RadialEig.l[:]
		for l in lList:
			node = f.getNode("/RadialEig/L%03i" % l)
			E = node.eigenvalues[:]
			V = node.eigenvectors[:]

			#assure correct phase convention (first oscillation should start out real positive)
			for i, curE in enumerate(E):
				bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
				phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
				V[:,i] *= exp(-1.0j * phase)

			eigenvalues.append(E)
			eigenvectors.append(V)

	finally:
		f.close()

	return lList, eigenvalues, eigenvectors	


def GetPopulationProductStates(psi, singleStates1, singleStates2):
	"""
	Calculates the population of psi in a set of single electron product states

	P_i =  |< SingleState1_i(1), SingleState2_j(2) | psi(1,2) >|^2

	singleStates 1 and 2 are lists of angular momentum states, containing an array 
	of radial states for the given angular momentum number such as generated
	by SetupRadialEigenstates in the Helium SAE example
	
	the projection is carried out for every combination of singlestate1 and singlestate2i
	is returned in a similar structure
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
		if V1.size == 0:
			continue

		for l2, V2 in enumerate(singleStates2):
			if V2.size == 0:
				continue

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(psi, lfilter), dtype=int32)
			if len(angularIndices) == 0:
				continue
		
			#Get the population for every combination of v1 and v2
			projV = CalculatePopulationRadialProductStates(l1, V1, l2, V2, data, angularIndices)
			cursum = sum([p for i1, i2, p in projV])
			print l1, l2, len(projV), cursum
			population.append((l1, l2, projV))

	return population


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


def GetDoubleAngularDistribution(psi, Z, interpEnergies, ionEnergies, ionStates, assocLegendre, theta):
	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	repr = psi.GetRepresentation()
	repr.MultiplyIntegrationWeights(tempPsi)
	distr = psi.GetRepresentation().GetDistributedModel()

	data = tempPsi.GetData()
	population = []

	angularRank = 0
	angRepr = repr.GetRepresentation(angularRank)

	dE = diff(interpEnergies)[0]**2
	dTh = diff(theta)[0]**2 * outer(sin(theta), sin(theta)) * (2*pi)**2

	cg = pyprop.core.ClebschGordan()

	interpK = sqrt(2 * interpEnergies)
	interpCount = len(interpEnergies)

	thetaCount = assocLegendre.shape[1]
	stateCount = ionStates[0].shape[1]
	angularDistr = zeros((thetaCount, thetaCount, interpCount, interpCount), dtype=double)

	assocLegendre = array(assocLegendre, dtype=complex)

	pop = 0
	M = 0
	#for some reason there is no population in m shells != 0. I don't know why, but it saves a lot of comp. work
	#for m in range(-5,5+1): 
	for m in [0]:
		angularDistrProj = zeros(angularDistr.shape, dtype=complex)
		for l1, V1 in enumerate(ionStates):
			E1 = array(ionEnergies[l1])
			print "%i/%i" % (l1, len(ionStates))
		
			for l2, V2 in enumerate(ionStates):
				E2 = array(ionEnergies[l2])
		
				#filter out coupled spherical harmonic indices. this gives us a set of L's for the given l1, l2, M
				lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2  \
							and l2>=abs(m) and l1>=abs(m) and coupledIndex.M == M
				angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(psi, lfilter), dtype=int)
				coupledIndices = map(angRepr.Range.GetCoupledIndex, angularIndices)

				if len(angularIndices) == 0:
					continue
			
				#calculate projection on radial states
				def doProj():
					radialProj = CalculateProjectionRadialProductStates(l1, V1, l2, V2, data, angularIndices)
					return radialProj
				radialProj = doProj()

				#scale states with 1/dE_1 dE_2
				def GetDensity(curE):
					interiorSpacing = list(diff(curE)[1:])
					leftSpacing = (curE[1] - curE[0])
					rightSpacing = (curE[-1] - curE[-2])
					spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
					return 1.0 / sqrt(spacing)
				stateDensity = outer(GetDensity(E1), GetDensity(E2))

				#coulomb phases (-i)**(l1 + l2) * exp( sigma_l1 * sigma_l2 )
				phase1 = exp(-1.0j * array([GetCoulombPhase(l1, Z/curK) for curK in sqrt(2*E1)]))
				phase2 = exp(-1.0j * array([GetCoulombPhase(l2, Z/curK) for curK in sqrt(2*E2)]))
				phase = (-1.j)**(l1 + l2) * outer(phase1, phase2)

				#interpolate projection on equidistant energies and sum over L
				interpProj = zeros((interpCount, interpCount), dtype=complex)
				proj = zeros((len(E1), len(E2)), dtype=complex)
				for j in range(radialProj.shape[0]):
					curRadialProj = phase * stateDensity * radialProj[j,:,:]

					#interpolate in polar complex coordinates
					def dointerp():
						r = abs(curRadialProj)**2
						i = arctan2(imag(curRadialProj), real(curRadialProj))
						argr = cos(i)
						argi = sin(i)
						interpr = scipy.interpolate.RectBivariateSpline(E1, E2, r, kx=1, ky=1)(interpEnergies, interpEnergies)
						interpArgR = scipy.interpolate.RectBivariateSpline(E1, E2, argr, kx=1, ky=1)(interpEnergies, interpEnergies)
						interpArgI = scipy.interpolate.RectBivariateSpline(E1, E2, argi, kx=1, ky=1)(interpEnergies, interpEnergies)
						interpPhase = (interpArgR + 1.j*interpArgI) / sqrt(interpArgR**2 + interpArgI**2)
						curInterpProj = sqrt(maximum(interpr, 0)) * interpPhase
						return curInterpProj
					curInterpProj = dointerp()

					#Sum over L-shells
					idx = coupledIndices[j]
					interpProj += curInterpProj * cg(idx.l1, idx.l2, m, 0, idx.L, M)
					proj += curRadialProj * cg(idx.l1, idx.l2, m, 0, idx.L, M)
		
				#sum up over l1, l2, E1, E2 to get total double ion prob
				pop += sum(abs(proj / stateDensity)**2)
				#expand spherical harmonics and add to angular distr proj
				def doSum():
					AddAngularProjectionAvgPhi(angularDistrProj, assocLegendre, interpProj, l1, l2, m, M)
				doSum()

		#calculate projection for this m-shell
		curAngularDistr = real(angularDistrProj * conj(angularDistrProj))
		curPop = sum(sum(sum(curAngularDistr ,axis=3),axis=2) * dTh ) * dE
		print "for m = %i, curAngularDistr = %f" % (m, curPop)
		angularDistr += curAngularDistr

	#estimate total ion prop by integrating over theta1, theta2, E1, E2
	totalPop = sum(sum(sum(angularDistr ,axis=3),axis=2) * dTh ) * dE
	print "Double Ionization Probability (integrated) = %s" % (totalPop)
	print "Double Ionization Probability (summed)     = %s" % (pop)


	return angularDistr

def GetParallelMomentumDistribution(energy, th, dp):
	nE = len(energy)
	dp2 = zeros((nE*2, nE*2), dtype=double)
	dp2[nE:, nE:] = dp[0,0,:,:]
	dp2[:nE, :nE] = fliplr(flipud(dp[-1,-1,:,:]))
	dp2[:nE, nE:] = fliplr(dp[0,-1,:,:])
	dp2[nE:, :nE] = flipud(dp[-1,0,:,:])

	e2 = array([-x for x in reversed(energy)] + [x for x in energy]) 

	return e2, dp2
