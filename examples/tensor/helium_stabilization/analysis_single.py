
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
		antiProb = GetSingleIonizationProbability(anti, boundStates, singleBoundStates, singleIonStates)
		return (symProb, antiProb)

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


def RunGetSingleIonizationBoundEnergyDistribution(fileList, removeBoundStates=True):
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
		return GetSingleIonizationBoundEnergyDistribution(psi, boundStates, singleBoundStates, singleBoundEnergies, singleIonStates, singleIonEnergies, dE, maxE)

	ionE, boundE, dpde = zip(*map(getdPdE, fileList))
	return ionE[0], boundE[0], dpde





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
	return absorbedProbability, ionizationProbability, singleIonizationProbability, ionizationProbability-singleIonizationProbability


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


def GetSingleIonizationBoundEnergyDistribution(psi, boundStates, singleBoundStates, singleBoundEnergies, singleIonStates, singleIonEnergies, dE, maxE):
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
	boundCount = sum([lBound.shape[1] for lBound in singleBoundStates])
	E = r_[0:maxE:dE]
	dpde = []
	boundE = []

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
			dpde += [interp(E, singleIonEnergies[lIon][:-1], curPop)]
			boundE += [singleBoundEnergies[lBound][iBound]]

	idx = argsort(boundE)
	dpde = array(dpde)[idx]
	boundE = array(boundE)[idx]

	outDp = [dpde[0,:]]
	outBoundE= [boundE[0]]
	for i in range(1,boundCount):
		if abs(boundE[i]-boundE[i-1]) > 1e-8:
			outDp += [dpde[i,:]]
			outBoundE += [boundE[i]]
		else:
			outDp[-1] += dpde[i,:]
		
	return E, array(outBoundE), array(outDp)




