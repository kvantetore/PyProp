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
	filenameTemplate = "raymond/stabilization_freq_5_scan_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0/stabilization_I_%i_kb20_dt_1e-02.h5"
	outputPrefix = "stabilization_scan2"
	intensity = r_[25:36]

	filenames = [filenameTemplate % i for i in intensity]
	absorb, totalIon, singleIon, doubleIon = RunSingleIonizationScan(filenames, outputPrefix)

	pylab.plot(intensity, totalIon, "-", label="Total Ion.")
	pylab.plot(intensity, singleIon, "--", label="Single Ion.")
	pylab.plot(intensity, doubleIon, ":", label="Double Ion.")
	xlabel("Intensity")
	title("Ionization Probability for w=5")
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

def RunSingleIonizationScan(filenames, outputPrefix):
	absorb, ion, singleIon = RunGetSingleIonizationProbability(filenames)

	absorb = array(absorb)
	totalIon = array(ion)
	singleIon = array(singleIon)
	doubleIon = totalIon - singleIon

	f = tables.openFile("%s.h5" % outputPrefix, "w")
	try:
		f.createArray(f.root, "absorb", absorb)
		f.createArray(f.root, "totalIon", totalIon)
		f.createArray(f.root, "singleIon", singleIon)
		f.createArray(f.root, "doubleIon", doubleIon)
	finally:
		f.close()

	return absorb, totalIon, singleIon, doubleIon


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


def RunGetSingleIonizationProbability(fileList):
	"""
	Calculates total and single ionization probability
	for a list of wavefunction files
	"""

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	boundEnergies, boundStates = GetBoundStates(config=conf)

	#Get single particle states
	isIonized = lambda E: E > 0.0
	isBound = lambda E: not isIonized(E)
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)

	def getIonProb(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		SymmetrizeWavefunction(psi, True)
		return GetSingleIonizationProbability(psi, boundStates, singleBoundStates, singleIonStates)

	return zip(*map(getIonProb, fileList))


def RunGetSingleIonizationEnergyDistribution(fileList):
	"""
	Calculates the energy distribution (dP/dE) of the 
	single ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	maxE = 15.
	dE = 0.1

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	boundEnergies, boundStates = GetBoundStates(config=conf)

	#Get single particle states
	isIonized = lambda E: E > 0.0
	isBound = lambda E: not isIonized(E)
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)

	#Calculate Energy Distribution (dP/dE)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStates, singleIonStates, singleIonEnergies, dE, maxE)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


def RunGetDoubleIonizationEnergyDistribution(fileList):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	maxE = 15.
	dE = 0.2

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	boundEnergies, boundStates = GetBoundStates(config=conf)

	#Get single particle states
	isIonized = lambda E: 0.0 < E <= maxE
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he+", isIonized, config=conf)

	#Calculate Energy Distribution (dP/dE1 dE2)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetDoubleIonizationEnergyDistribution(psi, boundStates, singleIonStates, singleIonEnergies, dE, maxE)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde



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
	RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	#calculate populations in product states containing bound he+ states
	#populations = GetPopulationSingleParticleStates(psi, singleBoundStates)
	populationsOld = GetPopulationProductStatesOld(psi, singleBoundStates, singleIonStates)
	populations = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)

	#Calculate single ionization probability
	lpop = [sum([p for i1, i2, p in pop]) for l1, l2, pop in populations]
	singleIonizationProbability = sum(lpop)

	lpopOld = [sum([p for i1, i2, p in pop]) for l1, l2, pop in populationsOld]
	singleIonizationProbabilityOld = sum(lpopOld)


	print "Absorbed Probability     = %s" % (absorbedProbability)
	print "Ioniziation Probability  = %s" % (ionizationProbability)
	print "Single Ionization Prob.  = %s, %s" % (singleIonizationProbability, singleIonizationProbabilityOld)
	print "Double Ionization Prob.  = %s" % (ionizationProbability - singleIonizationProbability)
	print "Single Ionization ratio  = %s" % (singleIonizationProbability/ionizationProbability)

	return absorbedProbability, ionizationProbability, singleIonizationProbability


def GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStates, singleIonStates, singleIonEnergies, dE, maxE):
	"""
	Calculates dP/dE for the single ionized part of the wavefunction, by first
	projecting away doubly bound states, and then projecting on single particle
	product states, and binning the probability by energy.
	"""

	#get absorbed prob
	absorbedProbability = 1.0 - real(psi.InnerProduct(psi))

	#remove boundstate projection
	RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	#calculate populations in product states containing bound he+ states
	populations = GetPopulationProductStatesOld(psi, singleBoundStates, singleIonStates)

	def getProbabilityL(startE, stopE, lPop, lEnergy):
		return sum([rPop for boundIndex, ionIndex, rPop in lPop if startE <= lEnergy[ionIndex] < stopE])

	def getEnergyGapProbability(startE, stopE):
		return sum([getProbabilityL(startE, stopE, lPop, singleIonEnergies[lIon]) for lBound, lIon, lPop in populations])

	E = r_[0:maxE:dE]
	dpde = map(getEnergyGapProbability, E, E+dE)

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
	RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	populations = GetPopulationProductStatesOld(psi, singleIonStates, singleIonStates)

	def getProbabilityL(startE1, stopE1, startE2, stopE2, lPop, lEnergy1, lEnergy2):
		return sum([rPop for i1, i2, rPop in lPop if (startE1 <= lEnergy1[i1] < stopE1) and (startE2 <= lEnergy2[i2] < stopE2)])

	def getEnergyGapProbability(startE1, stopE1, startE2, stopE2):
		return sum([getProbabilityL(startE1, stopE1, startE2, stopE2, lPop, singleEnergies[l1], singleEnergies[l2]) for l1, l2, lPop in populations])

	E = r_[0:maxE:dE]
	E1, E2 = meshgrid(E, E)
	startE1 = E1.flatten()
	stopE1 = startE1 + dE
	startE2 = E2.flatten()
	stopE2 = startE2 + dE

	dpde = array(map(getEnergyGapProbability, startE1, stopE1, startE2, stopE2)).reshape(len(E), len(E))

	return E, dpde



def GetSingleStatesFile(**args):
	radialGridPostfix = "_".join(GetRadialGridPostfix(**args))
	model = args.get("model", "he+")
	singleStatesFile = "output/singleelectron/eigenstates_sae_model_%s_%s.h5" % (model, radialGridPostfix)
	return singleStatesFile


def LoadSingleParticleStates(singleStatesFile):
	f = tables.openFile(singleStatesFile, "r")
	eigenvalues = []
	eigenvectors = []
	try:
		lList = f.root.RadialEig.l[:]
		for l in lList:
			node = f.getNode("/RadialEig/L%03i" % l)
			E = node.eigenvalues[:]
			V = node.eigenvectors[:]

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

	angularRank = 0
	angRepr = repr.GetRepresentation(0)
	angIndexGrid = repr.GetLocalGrid(angularRank)

	data = tempPsi.GetData()
	population = []

	m = 0
	clebschGordan = pyprop.core.ClebschGordan()
	l1, l2, L, M = zip(*map(lambda idx: angRepr.Range.GetCoupledIndex(int(idx)), angIndexGrid))
	m1 = (m,)*len(l1)
	m2 = array(M)-m
	cgList = map(clebschGordan, l1, l2, m1, m2, L, M)

	for l1, V1 in enumerate(singleStates1):
		print "%i/%i" % (l1, len(singleStates1))
		if V1.size == 0:
			continue

		for l2, V2 in enumerate(singleStates2):
			if V2.size == 0:
				continue

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, lfilter)
			#filter away all angular indices with zero clebsch-gordan coeff
			angularIndices = array(filter(lambda idx: abs(cgList[idx])>0, angularIndices), dtype=int32)
			if len(angularIndices) == 0:
				continue
		
			#Get the population for every combination of v1 and v2
			projV = CalculatePopulationRadialProductStates(l1, V1, l2, V2, data, angularIndices)
			cursum = sum([p for i1, i2, p in projV])
			print l1, l2, len(projV), cursum
			population.append((l1, l2, projV))

	return population


def GetPopulationProductStatesOld(psi, singleStates1, singleStates2):
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

	angularRank = 0
	angRepr = repr.GetRepresentation(0)
	angIndexGrid = repr.GetLocalGrid(angularRank)

	data = tempPsi.GetData()
	population = []

	m = 0
	clebschGordan = pyprop.core.ClebschGordan()
	l1, l2, L, M = zip(*map(lambda idx: angRepr.Range.GetCoupledIndex(int(idx)), angIndexGrid))
	m1 = (m,)*len(l1)
	m2 = array(M)-m
	cgList = map(clebschGordan, l1, l2, m1, m2, L, M)

	for l1, V1 in enumerate(singleStates1):
		print "%i/%i" % (l1, len(singleStates1))
		if V1.size == 0:
			continue

		for l2, V2 in enumerate(singleStates2):
			if V2.size == 0:
				continue

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, lfilter)
			#filter away all angular indices with zero clebsch-gordan coeff
			angularIndices = filter(lambda idx: abs(cgList[idx])>0, angularIndices)
			if len(angularIndices) == 0:
				continue
			
			def getPopulation(i1, i2):
				"""
				Calculate <v1(1), v2(2) | psi(1,2)> 
				"""
				#Sum over all local indices

				#sum coherently over L (product states in pure sph-harm)
				#getRadialProjection = lambda angIdx: dot(conj(V1[:,i1]), dot(conj(V2[:,i2]), data[angIdx, :, :]) ) * cgList[angIdx] 
				#popList = map(getRadialProjection, angularIndices)
				#pop = abs(sum(popList))**2

				##sum incoherently over L (product states in pure sph-harm) (supposedly wrong, but gives convincing results)
				#getRadialProjection = lambda angIdx: dot(conj(V1[:,i1]), dot(conj(V2[:,i2]), data[angIdx, :, :]) ) * cgList[angIdx] 
				#popList = map(getRadialProjection, angularIndices)
				#pop = sum(map(lambda x: abs(x)**2, popList))

				#sum incoherently over L (product states in coupled sph-harm)
				getRadialProjection = lambda angIdx: dot(conj(V2[:,i2]), dot(conj(V1[:,i1]), data[angIdx, :, :]) )  
				popList = map(getRadialProjection, angularIndices)
				pop = sum(map(lambda x: abs(x)**2, popList))

				#Sum over all processors
				pop = 2 * real(distr.GetGlobalSum(pop))
				return i1, i2, pop
		
			#Get the population for every combination of v1 and v2
			projV = map(getPopulation, *zip(*[(i1, i2) for i1 in range(V1.shape[1]) for i2 in range(V2.shape[1])]))
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

	angularRank = 0
	angRepr = repr.GetRepresentation(0)
	angIndexGrid = repr.GetLocalGrid(angularRank)

	data = tempPsi.GetData()
	population = []

	m = 0
	clebschGordan = pyprop.core.ClebschGordan()
	l1, l2, L, M = zip(*map(lambda idx: angRepr.Range.GetCoupledIndex(int(idx)), angIndexGrid))
	m1 = (m,)*len(l1)
	m2 = array(M)-m
	cgList = map(clebschGordan, l1, l2, m1, m2, L, M)

	for l1, V1 in enumerate(singleStates1):
		print "%i/%i" % (l1, len(singleStates1))

		for l2, V2 in enumerate(singleStates2):

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l2 == l1 and coupledIndex.l2 == l2
			angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, lfilter)
			
			def removePopulation(i1, i2):
				"""
				Calculate <v1(1), v2(2) | psi(1,2)> 
				"""
				#Sum over all local indices
				pop = 0
				for angIdx in angularIndices:
					cg = cgList[angIdx]
					if abs(cg) > 0:
						pop += real(abs( dot(conj(V1[:,i1]), dot(conj(V2[:,i2]), data[angIdx, :, :]) ) * cg )**2)
				#Sum over all processors
				pop = 2 * real(distr.GetGlobalSum(pop))
				return i1, i2, pop
		
			#Get the population for every combination of v1 and v2
			projV = map(getPopulation, *zip(*[(i1, i2) for i1 in range(V1.shape[1]) for i2 in range(V2.shape[1])]))
			population.append(projV)

	return population
