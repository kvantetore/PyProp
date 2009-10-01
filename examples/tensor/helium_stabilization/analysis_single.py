
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
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("h", isIonized, config=conf)
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
	
	maxE = 14.
	dE = 0.005

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
	singleIonEnergies, singleIonStates = GetFilteredSingleParticleStates("he", isIonized, config=conf)
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


def RunGetSingleIonizationAngularDistribution(fileList, dE=None, dTheta=None, removeBoundStates=True):
	"""
	Calculates the double differential energy distribution (dP/dE dOmega) of the 
	singly ionized continuum for a list of wavefunctions 
	files by projecting onto products of single particle states.
	"""

	#single ionized continuum is ~ hydrogen continuum
	lmax=10
	Z = 2-1
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
	isFilteredIonized = lambda E: 0.0 < E < 15
	isBound = lambda E: E <= 0.
	ionEnergies, ionStates = GetFilteredSingleParticleStates("h", isIonized, config=conf)
	singleBoundEnergies, singleBoundStates = GetFilteredSingleParticleStates("he+", isBound, config=conf)

	if dE == None:
		dE = maximum(min(diff(ionEnergies[0])), 0.05)
	interpE = r_[dE:15:dE]
	interpK = sqrt(interpE * 2)

	#Get spherical harmonics
	assocLegendre = GetAssociatedLegendrePoly(lmax, theta)

	angDistr = []
	for i, filename in enumerate(fileList):
		psi = pyprop.CreateWavefunctionFromFile(filename)
		sym, anti = GetSymmetrizedWavefunction(psi)
		print "norm_orig = %s" % psi.GetNorm()**2
		print "norm_sym  = %s" % sym.GetNorm()**2
		psi = sym
	
		#get absorbed prob
		absorbedProbability = 1.0 - real(psi.InnerProduct(psi))
		
		#remove boundstate projection
		if boundStates != None:
			RemoveBoundStateProjection(psi, boundStates)
		ionizationProbability = real(psi.InnerProduct(psi))

		#calculate angular distr for double ionized psi
		angDistr.append(GetSingleAngularDistribution(psi, Z, interpE, singleBoundEnergies, singleBoundStates, ionEnergies, ionStates, assocLegendre, theta))

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

	#Calculate single ionization probability, factor 2
	#because of implicit symmetrization of prod. state
	lpop = [sum([p for i1, i2, p in pop]) for l1, l2, pop in populations]
	singleIonizationProbability = 2.0 * sum(lpop)
	
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

		#iterate over all bound states in this l-shell, factor 2
		#because of implicit symmetrization of prod. state
		pop = 2 * array([d[2] for d in lPop]).reshape(nBound, nIon)
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

		#iterate over all bound states in this l-shell, factor 2
		#because of implicit symmetrization of prod. state
		pop = 2 * array([d[2] for d in lPop]).reshape(nBound, nIon)
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


def GetSingleAngularDistribution(psi, Z, interpEnergies, boundEnergies, boundStates, ionEnergies, ionStates, assocLegendre, theta):
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

	dE = diff(interpEnergies)[0]
	dTh = diff(theta)[0] * sin(theta) * 2*pi

	cg = pyprop.core.ClebschGordan()

	interpK = sqrt(2 * interpEnergies)
	interpCount = len(interpEnergies)

	thetaCount = assocLegendre.shape[1]
	stateCount = ionStates[0].shape[1]
	angularDistr = zeros((thetaCount, interpCount), dtype=double)

	assocLegendre = array(assocLegendre, dtype=complex)

	pop = 0
	M = 0
	#for some reason there is no population in m shells != 0. I don't know why, but it saves a lot of comp. work
	for m in range(-5,5+1): 
	#for m in [0]:
		for l1, V1 in enumerate(boundStates):
			E1 = array(boundEnergies[l1])
			unitVector1 = ones(len(E1))
			print "%i/%i" % (l1, len(boundStates))
			
			angularDistrProj = zeros((len(E1),) + angularDistr.shape, dtype=complex)
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
				stateDensity = outer(unitVector1, GetDensity(E2))

				#coulomb phases (-i)**(l) * exp( sigma_l )
				phase2 = (-1.j)**(l2) * exp(1.0j * array([GetCoulombPhase(l2, -Z/curK) for curK in sqrt(2*E2)]))
				phase = outer(unitVector1, phase2)

				#interpolate projection on equidistant energies in E2 and sum over L
				interpProj = zeros((len(E1), interpCount), dtype=complex)
				proj = zeros((len(E1), len(E2)), dtype=complex)
				for j in range(radialProj.shape[0]):
					curRadialProj = phase * stateDensity * radialProj[j,:,:]

					#interpolate in polar complex coordinates
					def dointerp():
						r = abs(curRadialProj)**2
						i = arctan2(imag(curRadialProj), real(curRadialProj))
						argr = cos(i)
						argi = sin(i)
						interpr = zeros((curRadialProj.shape[0], len(interpEnergies)), dtype=double)
						interpArgR = zeros((curRadialProj.shape[0], len(interpEnergies)), dtype=double)
						interpArgI = zeros((curRadialProj.shape[0], len(interpEnergies)), dtype=double)
						for idx in xrange(curRadialProj.shape[0]):
							interpr[idx,:] = scipy.interpolate.interp1d(E2, r[idx,:])(interpEnergies)
							interpArgR[idx,:] = scipy.interpolate.interp1d(E2, argr[idx,:])(interpEnergies)
							interpArgI[idx,:] = scipy.interpolate.interp1d(E2, argi[idx,:])(interpEnergies)
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
					AddSingleAngularProjectionAvgPhi(angularDistrProj, assocLegendre, interpProj, l2, M-m)
				doSum()

			#calculate projection for this l(bound), m-shell
			curAngularDistr = sum(real(angularDistrProj * conj(angularDistrProj)), axis=0)
			curPop = sum(sum(curAngularDistr, axis=1) * dTh ) * dE
			print "for m = %i, l1 = %i, curAngularDistr = %.10f" % (m, l1, curPop)
			angularDistr += curAngularDistr

	#estimate total ion prop by integrating over theta1, theta2, E1, E2
	totalPop = sum(sum(angularDistr, axis=1) * dTh ) * dE
	print "Double Ionization Probability (integrated) = %s" % (totalPop)
	print "Double Ionization Probability (summed)     = %s" % (pop)


	return angularDistr


