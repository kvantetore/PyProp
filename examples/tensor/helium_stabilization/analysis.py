
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
	Returns the local indices
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


#------------------------------------------------------------------------
#                        Product State Analysis
#------------------------------------------------------------------------

def RunSingleIonizationNikolopoulos():
	filenames = ["raymond/example_nikolopoulos/nikolopoulos_ionization_omega_165.h5_%04i.h5" % i for i in range(300)]
	absorb, ion, singleIon = RunGetSingleIonizationProbability()

def RunGetSingleIonizationProbability(wavefunctionFile):
	if isinstance(wavefunctionFile, list):
		fileList = wavefunctionFile
	elif isinstance(wavefunctionFile, str):
		fileList = [wavefunctionFile]
	else:
		raise Exception("wavefunctionFile should be list or str")

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	boundEnergies, boundStates = GetBoundStates(config=conf)

	#load single particle states
	singleStatesFile = GetSingleStatesFile(config=conf)
	lList, singleEnergies, singleStates = LoadSingleParticleStates(singleStatesFile)
	
	#filter bound he+ states
	doubleIonThreshold = 0.0
	getBoundStates = lambda E, V: transpose(array([V[:, i] for i in range(V.shape[1]) if E[i]<doubleIonThreshold]))
	getIonStates = lambda E, V: transpose(array([V[:, i] for i in range(V.shape[1]) if E[i]>doubleIonThreshold]))
	singleBoundStates = map(getBoundStates, singleEnergies, singleStates)
	singleIonStates = map(getIonStates, singleEnergies, singleStates)

	def getIonProb(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetSingleIonizationProbability(psi, boundStates, singleBoundStates, singleIonStates)

	if len(fileList[0]) == 1:
		return getIonProb(fileList[0])

	else:
		return zip(*map(getIonProb, fileList))


def	GetSingleIonizationProbability(psi, boundStates, singleBoundStates, singleIonStates):
	#get absorbed prob
	absorbedProbability = 1.0 - real(psi.InnerProduct(psi))

	#remove boundstate projection
	RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	#calculate populations in product states containing bound he+ states
	#populations = GetPopulationSingleParticleStates(psi, singleBoundStates)
	populations = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)

	#Calculate single ionization probability
	lpop = [sum([p[-1] for p in pop]) for pop in populations]
	singleIonizationProbability = sum(lpop)

	print "Absorbed Probability     = %s" % (absorbedProbability)
	print "Ioniziation Probability  = %s" % (ionizationProbability)
	print "Single Ionization Prob.  = %s" % (singleIonizationProbability)
	print "Double Ionization Prob.  = %s" % (ionizationProbability - singleIonizationProbability)
	print "Single Ionization ratio  = %s" % (singleIonizationProbability/ionizationProbability)

	return absorbedProbability, ionizationProbability, singleIonizationProbability


def GetSingleStatesFile(**args):
	radialGridPostfix = "_".join(GetRadialGridPostfix(**args))
	singleStatesFile = "output/singleelectron/eigenstates_sae_model_he+_%s.h5" % radialGridPostfix
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

def GetPopulationSingleParticleStates(psi, singleStates):
	"""
	Calculates projection of psi on a set of single electron states,
	and sums over all possible single particle states for the second electron

	P_i = sum_{j} < SingleState_i(2), j(1) | psi(1,2) >

	singleStates is a list of angular momentum states, containing an array 
	of radial states for the given angular momentum number such as generated
	by SetupRadialEigenstates in the Helium SAE example
	
	the projection is carried out for every such state, and the result 
	is returned in a similar structure
	"""

	raise Exception("This is not correct, we must first multiply S(r2) integrate r2, and then integrate ||. S(r1) .||^2")
	

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

	clebschGordan = pyprop.core.ClebschGordan()

	m = 0
	for l, V in enumerate(singleStates):
		#filter out coupled spherical harmonic indices corresponding to this l
		l2filter = lambda coupledIndex: coupledIndex.l2 == l
		angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, l2filter)
		
		def getPopulation(v0):
			"""
			gets the population of psi on v0 summed over particle 1
			"""
			#Sum over all local indices
			pop = 0
			for angIdx in angularIndices:
				globalAngIdx = int(angIndexGrid[angIdx])
				l1, l2, L, M = angRepr.Range.GetCoupledIndex(globalAngIdx)
				cg = clebschGordan(l1, l2, m, M-m, L, M)
				if abs(cg) > 0:
					for r1Idx in range(data.shape[1]):
						pop += real(abs( dot(conj(v0), data[angIdx, r1Idx, :]) * cg )**2)
			#Sum over all processors
			pop = real(distr.GetGlobalSum(pop))
			return pop
	
		#Get the population for every state in this l-shell
		if V.size > 0:
			projV = map(getPopulation, [V[:, i] for i in range(V.shape[1])])
		else:	
			projV = []
		population.append(projV)

	return population


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
			lfilter = lambda coupledIndex: coupledIndex.l2 == l1 and coupledIndex.l2 == l2 
			angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, lfilter)
			#filter away all angular indices with zero clebsch-gordan coeff
			angularIndices = filter(lambda idx: abs(cgList[idx])>0, angularIndices)
			
			def getPopulation(i1, i2):
				"""
				Calculate <v1(1), v2(2) | psi(1,2)> 
				"""
				#Sum over all local indices
				getRadialProjection = lambda angIdx: dot(conj(V1[:,i1]), dot(conj(V2[:,i2]), data[angIdx, :, :]) ) * cgList[angIdx] 
				popList = map(getRadialProjection, angularIndices)

				pop = abs(sum(popList))**2
				#pop = sum(map(lambda x: abs(x)**2, popList))
				
				#Sum over all processors
				pop = 2 * real(distr.GetGlobalSum(pop))
				return i1, i2, pop
		
			#Get the population for every combination of v1 and v2
			projV = map(getPopulation, *zip(*[(i1, i2) for i1 in range(V1.shape[1]) for i2 in range(V2.shape[1])]))
			population.append(projV)

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
