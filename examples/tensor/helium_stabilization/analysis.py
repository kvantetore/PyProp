
def RunSphericalHarmonicDistribution(wavefunctionFile):
	rightPsi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
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


#------------------------------------------------------------------------
#                        Product State Analysis
#------------------------------------------------------------------------

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
	localStart = distr.GetLocalStartIndex(nL, angularRank)
	localEnd = localStart + nL
	isLocal = lambda i: localStart <= i < localEnd
	globalIndices = filter(isLocal, r_[:nL])

	#Filter indices where coupled index is corresponding to supplied filter
	curFilter = lambda i: coupledIndexFilter(repr.Range.GetCoupledIndex(i))
	globalFilteredIndices = filter(l1Filter, globalIndices)

	#map global to local indices
	globalToLocal = lambda i: i - localStart
	localFilteredIndices = map(globalToLocal, globalFilteredIndices)

	return localFilteredIndices

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

	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	tempPsi.MultiplyIntegrationWeights()

	data = tempPsi.GetData()
	distr = psi.GetRepresentation().GetDistributedModel()
	population = []

	for l, V in range(len(singleStates)):
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
				for r1Idx in range(data.shape[1]):
					pop += abs( dot(conj(v0), data[angIdx, r1Idx, :]) )**2
			#Sum over all processors
			pop = distr.GlobalSum(pop)
			return pop
		
		#Get the population for every state in this l-shell
		projV = map(getPopulation, [V[:, i] for i in range(V.shape[1])])
		population.append(projV)

	return population
