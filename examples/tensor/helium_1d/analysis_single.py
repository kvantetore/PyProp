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
		boundEnergies = [0]
		gsPsi = pyprop.CreateWavefunctionFromFile(GroundstateFilename(config=conf), "/wavefunction")
		boundStates = [gsPsi]

	#Get single particle states for odd/even symmetry
	isEven = lambda v: not (abs(sum([a+b for a,b in zip(v, reversed(v))])) < 1e-12)
	isIonizedEven = lambda en, v: (0.0 <= en) and isEven(v)
	isIonizedOdd = lambda en, v: (0.0 <= en) and not isEven(v)
	isBoundEven = lambda en, v: (en < 0.0) and isEven(v)
	isBoundOdd = lambda en, v: (en < 0.0) and not isEven(v)

	singleIonEnergiesEven, singleIonStatesEven = GetFilteredSingleParticleStates(isIonizedEven, config=conf)
	singleIonEnergiesOdd, singleIonStatesOdd = GetFilteredSingleParticleStates(isIonizedOdd, config=conf)
	singleBoundEnergiesEven, singleBoundStatesEven = GetFilteredSingleParticleStates(isBoundEven, config=conf)
	singleBoundEnergiesOdd, singleBoundStatesOdd = GetFilteredSingleParticleStates(isBoundOdd, config=conf)

	def getIonProb(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		sym, anti = GetSymmetrizedWavefunction(psi)
		symProb = GetSingleIonizationProbability(sym, boundStates, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd)
		antiProb = GetSingleIonizationProbability(anti, boundStates, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd)
		return (symProb, antiProb)

	print
	print len(fileList)
	print
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
	else:
		boundEnergies = boundStates = None

	#Get single particle states for odd/even symmetry
	isEven = lambda v: not (abs(sum([a+b for a,b in zip(v, reversed(v))])) < 1e-12)
	isIonizedEven = lambda en, v: (0.0 <= en <= maxE) and isEven(v)
	isIonizedOdd = lambda en, v: (0.0 <= en <= maxE) and not isEven(v)
	isBoundEven = lambda en, v: (en < 0.0) and isEven(v)
	isBoundOdd = lambda en, v: (en < 0.0) and not isEven(v)

	singleIonEnergiesEven, singleIonStatesEven = GetFilteredSingleParticleStates(isIonizedEven, config=conf)
	singleIonEnergiesOdd, singleIonStatesOdd = GetFilteredSingleParticleStates(isIonizedOdd, config=conf)
	singleBoundEnergiesEven, singleBoundStatesEven = GetFilteredSingleParticleStates(isBoundEven, config=conf)
	singleBoundEnergiesOdd, singleBoundStatesOdd = GetFilteredSingleParticleStates(isBoundOdd, config=conf)

	#Calculate Energy Distribution (dP/dE)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		return GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd, dE, maxE)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


def GetSingleIonizationProbability(psi, boundStates, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd):
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

	#calculate populations in product states containing bound 1D states
	popListEvenOdd, popListOddEven = GetPopulationProductStates(psi, singleIonStatesEven, singleIonStatesOdd, singleIonStatesOdd, singleIonStatesEven)
	
	populations = [p for (i1,i2,p) in popListEvenOdd] + [p for (i1,i2,p) in popListOddEven]
	#singleIonizationProbability = sum(populations)
	doubleIonizationProbability = sum(populations)
	singleIonizationProbability = ionizationProbability - doubleIonizationProbability
	
	print "Absorbed Probability     = %s" % (absorbedProbability)
	print "Ioniziation Probability  = %s" % (ionizationProbability)
	print "Single Ionization Prob.  = %s" % (singleIonizationProbability)
	print "Double Ionization Prob.  = %s" % (ionizationProbability - singleIonizationProbability)
	print "Single Ionization ratio  = %s" % (singleIonizationProbability/ionizationProbability)

	return absorbedProbability, ionizationProbability, singleIonizationProbability, ionizationProbability-singleIonizationProbability


def GetSingleIonizationEnergyDistribution(psi, boundStates, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd, dE, maxE):
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
	populationsEvenOdd, populationsOddEven = GetPopulationProductStates(psi, singleBoundStatesEven, singleBoundStatesOdd, singleIonStatesEven, singleIonStatesOdd)

	#Create an array for dpde
	E = r_[0:maxE:dE]
	dpde = zeros(len(E), dtype=double)
	#dpde = zeros(len(singleIonEnergiesEven)-1, dtype=double)

	#In order to create an approx to dp/de_ion, we interpolate for each iBound,
	#and add coherently to the dpde array

	#number of states
	nBoundEven = shape(singleBoundStatesEven)[0]
	nBoundOdd = shape(singleBoundStatesOdd)[0]
	nIonEven = shape(singleIonStatesEven)[0]
	nIonOdd = shape(singleIonStatesOdd)[0]

	#iterate over all bound states (even-odd)
	pop = array([d[2] for d in populationsEvenOdd]).reshape(nBoundEven, nIonOdd)
	for iBound in range(nBoundEven):
		#interpolate over ionized populations
		curPop = pop[iBound, :-1] / diff(singleIonEnergiesOdd)
		dpde += interp(E, singleIonEnergiesOdd[:-1], curPop)

	#iterate over all bound states (odd-even)
	pop = array([d[2] for d in populationsOddEven]).reshape(nBoundOdd, nIonEven)
	for iBound in range(nBoundOdd):
		#interpolate over ionized populations
		curPop = pop[iBound, :-1] / diff(singleIonEnergiesEven)
		dpde += interp(E, singleIonEnergiesEven[:-1], curPop)		
	
	return E, dpde

