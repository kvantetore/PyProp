maxEnergy = 1.0
energyRes = 0.01

def RunGetDoubleIonizationEnergyDistribution(fileList, removeBoundStates=True, removeSingleIonStates=False):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""
	
	#maxE = 4.
	#dE = 0.005
	lmax = 10

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
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)
	#doubleIonEnergies, doubleIonStates = GetFilteredSingleParticleStates("he+", isFilteredIonized, config=conf)

	#Get energy normalized Coulomb waves
	psi = pyprop.CreateWavefunctionFromFile(fileList[0])
	doubleIonEnergies, doubleIonStates = SetupRadialCoulombStatesEnergyNormalized(psi, -2, maxEnergy, energyRes)

	#Calculate Energy Distribution (dP/dE1 dE2)
	def getdPdE(filename):
		psi = pyprop.CreateWavefunctionFromFile(filename)

		if removeSingleIonStates:
			RemoveProductStatesProjection(psi, singleBoundStates, singleIonStates)
			RemoveProductStatesProjection(psi, singleIonStates, singleBoundStates)

		return GetDoubleIonizationEnergyDistribution(psi, boundStates, doubleIonStates, doubleIonEnergies, energyRes, maxEnergy)

	E, dpde = zip(*map(getdPdE, fileList))
	return E[0], dpde


def RunGetDoubleIonizationAngularDistribution(fileList, dE=None, dTheta=None, removeBoundStates=True, removeSingleIonStates=False):
	"""
	Calculates the double differential energy distribution (dP/dE1 dE2) of the 
	doubly ionized continuum for a list of wavefunction 
	files by projecting onto products of single particle states.
	"""

	Z = 2
	lmax = 10
	if dTheta == None:
		dTheta = pi/19
	thetaCount = pi/dTheta+1
	theta = linspace(0, pi, thetaCount)

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
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)
	doubleIonEnergies, doubleIonStates = GetFilteredSingleParticleStates("he+", isFilteredIonized, config=conf)

	if dE == None:
		#dE = maximum(min(diff(singleIonEnergies[0])), 0.005)
		dE = energyRes
	interpE = r_[dE:maxEnergy:dE]
	interpK = sqrt(interpE * 2)

	#Get spherical harmonics
	assocLegendre = GetAssociatedLegendrePoly(lmax, theta)

	angDistr = []
	for i, filename in enumerate(fileList):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		#sym, anti = GetSymmetrizedWavefunction(psi)
		#print "norm_orig = %s" % psi.GetNorm()**2
		#print "norm_sym  = %s" % sym.GetNorm()**2
		#psi = sym
	
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

		#calculate angular distr for double ionized psi
		angDistr.append(GetDoubleAngularDistribution(psi, Z, interpE, doubleIonEnergies, doubleIonStates, assocLegendre, theta))

	return interpE, theta, angDistr


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
	for m in range(-5,5+1): 
		#if m == 0: continue

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
				phase1 = exp(1.0j * array([GetCoulombPhase(l1, -Z/curK) for curK in sqrt(2*E1)]))
				phase2 = exp(1.0j * array([GetCoulombPhase(l2, -Z/curK) for curK in sqrt(2*E2)]))
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
					interpProj += curInterpProj * cg(idx.l1, idx.l2, m, M-m, idx.L, M)
					proj += curRadialProj * cg(idx.l1, idx.l2, m, M-m, idx.L, M)
		
				#sum up over l1, l2, E1, E2 to get total double ion prob
				pop += sum(abs(proj / stateDensity)**2)
				#expand spherical harmonics and add to angular distr proj
				def doSum():
					AddDoubleAngularProjectionAvgPhi(angularDistrProj, assocLegendre, interpProj, l1, l2, m, M)
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

def GetDoubleParallelMomentumDistribution(energy, th, dp):
	nE = len(energy)
	dp2 = zeros((nE*2, nE*2), dtype=double)
	dp2[nE:, nE:] = dp[0,0,:,:]
	dp2[:nE, :nE] = fliplr(flipud(dp[-1,-1,:,:]))
	dp2[:nE, nE:] = fliplr(dp[0,-1,:,:])
	dp2[nE:, :nE] = flipud(dp[-1,0,:,:])

	e2 = array([-x for x in reversed(energy)] + [x for x in energy]) 

	return e2, dp2
