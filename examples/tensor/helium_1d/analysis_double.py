maxEnergy = 1.0
energyRes = 0.001

def RunGetDoubleIonizationEnergyDistribution(fileList, removeBoundStates=True, removeSingleIonStates=False):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	#maxE = 4.
	#dE = 0.005

	#load wavefunction
	conf = pyprop.LoadConfigFromFile(fileList[0])

	#load bound states
	if removeBoundStates:
		boundEnergies, boundStates = GetBoundStates(config=conf)
	else:
		boundEnergies = boundStates = None

	#Get single particle states
	isIonized = lambda E: 0.0 < E
	isFilteredIonized = lambda E: 0.0 < E < maxEnergy
	isBound = lambda E: E <= 0.
	singleIonEnergies, singleIonStates = GetEnergyFilteredSingleParticleStates(isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetEnergyFilteredSingleParticleStates(isBound, config=conf)
	doubleIonEnergies, doubleIonStates = GetEnergyFilteredSingleParticleStates(isFilteredIonized, config=conf)

	#Calculate Energy Distribution (dP/dE1 dE2)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)

		if removeSingleIonStates:
			RemoveProductStatesProjection(psi, singleBoundStates, singleIonStates)
			RemoveProductStatesProjection(psi, singleIonStates, singleBoundStates)

		return GetDoubleIonizationEnergyDistribution(psi, boundStates, doubleIonStates, doubleIonEnergies, energyRes, maxEnergy)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


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

	#number of states 
	n1 = len(singleIonStates)
	n2 = len(singleIonStates)

	#scale states with 1/dE_1 dE_2
	sliceStart = 1
	sliceStep = 2
	pop = array([d[2] for d in populations]).reshape(n1, n2)[sliceStart::sliceStep, sliceStart::sliceStep]
	E1 = singleEnergies[sliceStart::sliceStep]
	E2 = singleEnergies[sliceStart::sliceStep]
	meshE1, meshE2 = meshgrid(E1[:-1], E2[:-1])
	pop[:-1,:-1] /= outer(diff(E1), diff(E2))
	
	#2d interpolation over all states in this shell
	interpolator = scipy.interpolate.RectBivariateSpline(E1[:-1], E2[:-1], pop[:-1, :-1], kx=1, ky=1)
	dpde += interpolator(E, E)
	
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
