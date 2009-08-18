maxEnergy = 1.0
energyRes = 0.001

def RunGetDoubleIonizationEnergyDistribution(fileList, removeBoundStates=True, removeSingleIonStates=False):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states for odd/even symmetry
	isEven = lambda v: not (abs(sum([a+b for a,b in zip(v, reversed(v))])) < 1e-12)
	isIonizedEven = lambda en, v: (0.0 <= en <= maxEnergy) and isEven(v)
	isIonizedOdd = lambda en, v: (0.0 <= en <= maxEnergy) and not isEven(v)
	isBoundEven = lambda en, v: (en < 0.0) and isEven(v)
	isBoundOdd = lambda en, v: (en < 0.0) and not isEven(v)

	singleIonEnergiesEven, singleIonStatesEven = GetFilteredSingleParticleStates(isIonizedEven, config=conf)
	singleIonEnergiesOdd, singleIonStatesOdd = GetFilteredSingleParticleStates(isIonizedOdd, config=conf)
	singleBoundEnergiesEven, singleBoundStatesEven = GetFilteredSingleParticleStates(isBoundEven, config=conf)
	singleBoundEnergiesOdd, singleBoundStatesOdd = GetFilteredSingleParticleStates(isBoundOdd, config=conf)

	#Calculate Energy Distribution (dP/dE1 dE2)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)

		if removeSingleIonStates:
			RemoveProductStatesProjection(psi, singleBoundStatesEven, singleIonStatesOdd)
			RemoveProductStatesProjection(psi, singleBoundStatesOdd, singleIonStatesEven)
			RemoveProductStatesProjection(psi, singleIonStatesEven, singleBoundStatesOdd)
			RemoveProductStatesProjection(psi, singleIonStatesOdd, singleBoundStatesEven)

		return GetDoubleIonizationProbability(psi, boundStates, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd, energyRes, maxEnergy)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


def GetDoubleIonizationProbability(psi, boundStates, singleIonStatesEven, singleIonEnergiesEven, singleIonStatesOdd, singleIonEnergiesOdd, dE, maxE):
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

	populationsEvenEven, populationsOddOdd = GetPopulationProductStates(psi, singleIonStatesEven, singleIonStatesOdd, singleIonStatesOdd, singleIonStatesEven)

	E = r_[0:maxE:dE]
	dpde = zeros((len(E), len(E)), dtype=double)

	#number of states
	nIonEven = shape(singleIonStatesEven)[0]
	nIonOdd = shape(singleIonStatesOdd)[0]

	def calculateDpDe(popList, E, nStates):
		pop = array([d[2] for d in popList]).reshape(nStates, nStates)
		meshE1, meshE2 = meshgrid(E[:-1], E[:-1])
		pop[:-1,:-1] /= outer(diff(E), diff(E))
	
		#2d interpolation over all states with this symmetry
		interpolator = scipy.interpolate.RectBivariateSpline(E[:-1], E[:-1], pop[:-1, :-1], kx=1, ky=1)
		return interpolator

	#Calculate dpde for both symmetry combinations
	dpde += calculateDpDe(populationsEvenEven, singleIonEnergiesEven, nIonEven)(E, E)
	dpde += calculateDpDe(populationsOddOdd, singleIonEnergiesOdd, nIonOdd)(E, E)
	
	return E, dpde


def GetDoubleParallelMomentumDistribution(energy, th, dp):
	nE = len(energy)
	dp2 = zeros((nE*2, nE*2), dtype=double)
	dp2[nE:, nE:] = dp[0,0,:,:]
	dp2[:nE, :nE] = fliplr(flipud(dp[-1,-1,:,:]))
	dp2[:nE, nE:] = fliplr(dp[0,-1,:,:])
	dp2[nE:, :nE] = flipud(dp[-1,0,:,:])

	e2 = array([-x for x in reversed(energy)] + [x for x in energy]) 

	return e2, dp2
