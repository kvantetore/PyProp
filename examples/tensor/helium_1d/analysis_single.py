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
		#boundStates = [GetGroundstate(config=conf)]
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: E > 0.0
	isBound = lambda E: not isIonized(E)
	singleIonEnergies, singleIonStates = GetEnergyFilteredSingleParticleStates(isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetEnergyFilteredSingleParticleStates(isBound, config=conf)

	def getIonProb(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
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
	
	maxE = 5.
	dE = 0.00025

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
		#boundStates = [GetGroundstate(config=conf)]

	else:
		boundEnergies = boundStates = None


	#Get single particle states
	isIonized = lambda E: 0.0 <= E <= maxE
	isBound = lambda E: E < 0.0
	singleIonEnergies, singleIonStates = GetEnergyFilteredSingleParticleStates(isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetEnergyFilteredSingleParticleStates(isBound, config=conf)

	#Calculate Energy Distribution (dP/dE)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStates, singleIonStates, singleIonEnergies, dE, maxE)

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
	if boundStates != None:
		RemoveBoundStateProjection(psi, boundStates)
	ionizationProbability = real(psi.InnerProduct(psi))

	##calculate sum of populations in singly ionized product states
	popList = GetPopulationProductStates(psi, singleBoundStates, singleIonStates)
	
	populations = [p for (i1,i2,p) in popList]
	singleIonizationProbability = sum(populations)
	
	print "Absorbed Probability     = %s" % (absorbedProbability)
	print "Ioniziation Probability  = %s" % (ionizationProbability)
	print "Single Ionization Prob.  = %s" % (singleIonizationProbability)
	print "Double Ionization Prob.  = %s" % (ionizationProbability - singleIonizationProbability)
	print "Single Ionization ratio  = %s" % (singleIonizationProbability/ionizationProbability)

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
	#dpde = zeros(len(E), dtype=double)
	dpde = zeros(len(singleIonEnergies)-1, dtype=double)

	#In order to create an approx to dp/de_ion, we interpolate for each iBound,
	#and add coherently to the dpde array

	#number of states
	nBound = shape(singleBoundStates)[0]
	nIon = shape(singleIonStates)[0]

	#iterate over all bound states
	pop = array([d[2] for d in populations]).reshape(nBound, nIon)
	for iBound in range(nBound):
		#interpolate over ionized populations
		curPop = pop[iBound, :-1] / diff(singleIonEnergies)
		dpde += curPop
		#dpde += interp(E, singleIonEnergies[:-1], curPop)
		
	#return E, interp(E, singleIonEnergies[:-1], dpde)
	return singleIonEnergies[:-1], dpde

