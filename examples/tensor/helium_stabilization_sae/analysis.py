#---------------------------------------------------------------------------------------
#            Eigenstate Analysis
#---------------------------------------------------------------------------------------

def SaveRadialEigenstates(**args):
	postfix = "_".join(GetModelPostfix(**args) + GetRadialGridPostfix(**args))
	outputFile = "eigenstates/eigenstates_sae_%s.h5" % postfix

	#Setup problem
	prop = SetupProblem(**args)

	#Setp eigenvalues and eigenvectors
	eigenValues, eigenVectors = SetupRadialEigenstates(prop)
	lCount = len(eigenValues)

	#Save eigenvalues and eigenvectors to file
	if outputFile != None:
		f = tables.openFile(outputFile, "w")
		try:
			#Create main group
			eigGroup = f.createGroup(f.root, "RadialEig")

			#save config object
			eigGroup._v_attrs.configObject = prop.Config.cfgObj
			f.createArray(eigGroup, "l", r_[:lCount])

			#Save each l-shell in a separate group
			for l, (E, V) in enumerate(zip(eigenValues, eigenVectors)):
				lGroup = f.createGroup(eigGroup, "L%03i" % l)
				lGroup._v_attrs.l = l
				f.createArray(lGroup, "eigenvalues", E)
				f.createArray(lGroup, "eigenvectors", V)

		finally:
			f.close()


def SetupRadialBesselStatesEnergyNormalized(prop, Z, Emax, dE):
	E = r_[dE:Emax:dE]
	k = sqrt(E*2)

	psi = prop.psi
	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)

	
	#Setup Radial Waves
	S = SetupOverlapMatrix(prop)
	coulWaves = map(lambda curL: array(map(lambda curK: sqrt(2*dE/pi/curK)*GetRadialBesselWaveBSplines(Z, int(curL), curK, bspline, S), k)).transpose(), l)
	energies = [E]*len(l)

	return energies, coulWaves


def SetupRadialCoulombStatesMomentumNormalized(prop, Z, kmax, dk):
	k = r_[dk:kmax:dk]
	E = k**2 / 2

	psi = prop.psi
	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)
	
	#Setup Radial Waves
	S = SetupOverlapMatrix(prop)
	coulWaves = map(lambda curL: array(map(lambda curK: sqrt(2*dk/pi)*GetRadialCoulombWaveBSplines(Z, int(curL), curK, bspline, S), k)).transpose(), l)
	energies = [E]*len(l)

	return energies, coulWaves


def SetupRadialCoulombStatesEnergyNormalized(prop, Z, Emax, dE):
	E = r_[dE:Emax:dE]
	k = sqrt(E*2)

	psi = prop.psi
	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)
	
	#Setup Radial Waves
	S = SetupOverlapMatrix(prop)
	coulWaves = map(lambda curL: array(map(lambda curK: sqrt(2*dE/pi/curK)*GetRadialCoulombWaveBSplines(Z, int(curL), curK, bspline, S), k)).transpose(), l)
	energies = [E]*len(l)

	return energies, coulWaves


def SetupRadialEigenstates(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	SAE Helium system.

	The eigenvalues are found by setting up a radial matrix for each l-value
	and using the generalized eigenvalue solver in scipy to find all
	eigenvalues and vectors. 
	
	eigenvalues is a list of 1-d eigenvalue arrays. Each array corresponding
	to a
	"""
	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)
	bspl = prop.psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	eigenValues = []
	eigenVectors = []

	lCount = prop.psi.GetData().shape[0]
	eigenvectorScaling = 1

	for l in range(lCount):
		l = int(l)
		M = SetupRadialMatrix(prop, potentialIndices, l)

		E, V = scipy.linalg.eig(a=M, b=S)

		idx = argsort(real(E))
		E = real(E[idx])
		eigenValues.append(E)

		#Sort and normalize eigenvectors
		sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
		V = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()
		eigenVectors.append(V)

		#assure correct phase convention (first oscillation should start out real positive)
		for i, curE in enumerate(E):
			bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
			phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
			V[:,i] *= exp(-1.0j * phase)

	return eigenValues, eigenVectors

def SetRadialEigenstate(psi, eigenVectors, n, l, sourceScaling=0., destScaling=1.0):
	"""
	Sets psi to an eigenvector from a list of eigenvectors as calculated by
	SetupRadialEigenstates()

	n, l is the quantum number of the state to set in psi
	"""
	radialIndex = n - l - 1
	vec = eigenVectors[l][:, radialIndex]
	psi.GetData()[:] *= sourceScaling
	psi.GetData()[l, :] += destScaling * vec


def CalculateRadialCorrelation(psi, eigenVectors, n, l, overlap):
	radialIndex = n - l - 1
	psi.Clear()
	vec = eigenVectors[l][:, radialIndex]
	return sum( conj(vec) * dot(overlap, psi.GetData()[l, :]) )
	

def CalculateEnergyDistribution(psi, eigenValues, eigenVectors, overlap, dE=None):
	if not dE:
		dE = min(diff(sorted(eigenValues[0])))
	minE = 0
	maxE = eigenValues[0][3*len(eigenValues[0])/4]
	E = r_[minE:maxE:dE]
	energyDistr = zeros(len(E), dtype=double)
	#energyDistr = []

	totalIon = 0 
	for l, (curE, curV) in enumerate(zip(eigenValues, eigenVectors)):
		idx = where(curE>0)[0]
		curE = curE[idx]
		curV = curV[:,idx]
		idx = argsort(curE)
		curE = curE[idx]
		curV = curV[:,idx]

		#Get projection on eigenstates
		psiSlice = psi.GetData()[l, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = abs(dot(conj(curV[:,idx].transpose()), overlapPsi))**2
		totalIon += sum(proj)

		#Calculate density of states
		interiorSpacing = list((diff(curE[:-1]) + diff(curE[1:])) / 2.)
		leftSpacing = (curE[1] - curE[0]) / 2.
		rightSpacing = (curE[-1] - curE[-2]) / 2.
		spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
		density = 1.0 / spacing

		#Interpolate to get equispaced dP/dE
		#energyDistr.append( scipy.interp(E, curE, proj * density, left=0, right=0) )
		energyDistr +=  scipy.interp(E, curE, proj * density, left=0, right=0) 

	totalIon2 = sum(sum(array(energyDistr), axis=0)) * dE
	print totalIon
	print totalIon2
	interpolateError = totalIon - totalIon2 

	return E, energyDistr

def CalculateAngularDistribution(psi, eigenValues, eigenVectors, overlap, interpMethod="polar-square"):
	dE = maximum(min(diff(eigenValues[0])), 0.1)
	minE = dE
	maxE = 16 #eigenValues[0][-1] #eigenValues[0][3*len(eigenValues[0])/4]
	E = r_[minE:maxE:dE]
	energyDistr = []

	Z = 1
	thetaCount = 10
	lmax = len(eigenValues) -1
	theta = linspace(0, pi, thetaCount)
	leg = GetLegendrePoly(lmax, theta)

	angularDistrProj = zeros((thetaCount, len(E)), dtype=complex)

	for l, (curE, curV) in enumerate(zip(eigenValues, eigenVectors)):
		idx = where(curE>=0)[0]
		idx = idx[where(curE[idx]<=maxE)[0]]
		
		curE = curE[idx]
		curk = sqrt(curE/2)

		#Phase for outgoing waves
		sigma = array([GetCoulombPhase(l, Z / k) for k in curk])
		phase = (-1.j)**l * exp(1.j * sigma)

		#Get projection on eigenstates
		psiSlice = psi.GetData()[l, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = dot(conj(curV[:,idx].transpose()), overlapPsi)

		#Calculate density of states
		#interiorSpacing = list((diff(curE[:-1]) + diff(curE[1:])) / 2.)
		interiorSpacing = list(diff(curE)[1:])
		leftSpacing = (curE[1] - curE[0])
		rightSpacing = (curE[-1] - curE[-2])
		spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
		density = 1.0 / sqrt(spacing)

		#Interpolate to get equispaced dP/dE
		partialProj = phase * proj * density 
		if interpMethod == "polar-square":
			#interpolate in complex polar coordinates
			r = abs(partialProj)**2
			i = arctan2(imag(partialProj), real(partialProj))
			argr = cos(i)
			argi = sin(i)
			interpr = scipy.interpolate.UnivariateSpline(curE, r, s=0)(E)
			interpArgR = scipy.interpolate.UnivariateSpline(curE, argr, s=0)(E)
			interpArgI = scipy.interpolate.UnivariateSpline(curE, argi, s=0)(E)
			#plot(interpArgR)
			#plot(interpArgI)
			interpPhase = (interpArgR + 1.j*interpArgI) / sqrt(interpArgR**2 + interpArgI**2)
			interpProj = sqrt(maximum(interpr, 0)) * interpPhase
		elif interpMethod == "polar-square-old":
			#interpolate in complex polar coordinates
			r = abs(partialProj)**2
			i = arctan2(imag(partialProj), real(partialProj))
			interpr = scipy.interpolate.UnivariateSpline(curE, r, s=0)(E)
			interpi = scipy.interpolate.UnivariateSpline(curE, i, s=100)(E)
			plot(curE, i)
			plot(E, interpi)
			#interpr = scipy.interp(E, curE, r, right=0, left=0)
			#interpi = scipy.interp(E, curE, i, right=0, left=0)
			interpProj = sqrt(maximum(interpr,0)) * exp(1.0j * interpi)
		elif interpMethod == "polar":
			r = abs(partialProj)
			i = arctan2(imag(partialProj), real(partialProj))
			interpr = scipy.interp(E, curE, r, right=0, left=0)
			interpi = scipy.interp(E, curE, i, right=0, left=0)
			interpProj = interpr * exp(1.0j * interpi)
		elif interpMethod == "cartesian":
			#interpolate in complex cartesian coordinates
			r = real(partialProj)
			i = imag(partialProj)
			interpr = scipy.interp(E, curE, r, left=0, right=0) 
			interpi = scipy.interp(E, curE, i, left=0, right=0)
			interpProj= interpr + 1.j*interpi
		else:
			raise Exception("unknown interpMethod %s" % interpMethod)

		print sum(abs(interpProj)**2) * dE, sum(abs(partialProj/density)**2)
		angularDistrProj += outer(leg[l,:], interpProj)

	return E, theta, abs(angularDistrProj)**2


def GetLegendrePoly(lmax, theta):
	leg = []
	m = 0
	for l in range(lmax+1):
		leg.append(scipy.special.sph_harm(m, l, 0, theta))
	return array(leg)


def CalculateBoundDistribution(psi, eigenValues, eigenVectors, overlap, boundThreshold=0):
	boundDistr = []
	boundE = []
	boundV = []
	boundTotal = 0

	for l, (curE, curV) in enumerate(zip(eigenValues, eigenVectors)):
		boundIdx = where(curE < boundThreshold)[0]

		#Get projection on eigenstates
		psiSlice = psi.GetData()[l, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = dot(conj(curV[:,boundIdx].transpose()), overlapPsi)

		#Interpolate to get equispaced dP/dE
		boundDistr.append( proj )
		boundE.append( curE[boundIdx] )
		boundV.append(curV[:, boundIdx])
		boundTotal += sum(abs(proj)**2)

	return boundE, boundV, boundDistr, boundTotal
	

def RemoveBoundDistribution(psi, E, V, S):
	boundE, boundV, boundDistr, boundTotal = CalculateBoundDistribution(psi, E, V, S)
	
	for l, (curE, curV, curDist) in enumerate(zip(boundE, boundV, boundDistr)):
		psi.GetData()[l, :] -= dot(curV, curDist)

	return boundE, boundV, boundDistr, boundTotal


def SetupRadialMatrix(prop, whichPotentials, angularIndex):
	matrixSize = prop.psi.GetData().shape[1]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		angularBasisPairs = potential.BasisPairs[0]
		idx = [idx for idx, (i,j) in enumerate(zip(angularBasisPairs[:,0], angularBasisPairs[:,1])) if i==j==angularIndex]
		if len(idx) != 1:
			raise "Invalid angular indices %s" % idx
		idx = idx[0]

		basisPairs = potential.BasisPairs[1]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[idx, i]

	return matrix


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupRadialMatrix(prop, [overlap], 0)
	return matrix
	

#---------------------------------------------------------------------------------------
#            Coulomb Wave Analysis
#---------------------------------------------------------------------------------------

def GetRadialCoulombWave(Z, l, k, rmax=100):
	#Get the Coulomb function in grid space
	dr = 0.01
	r = r_[dr:rmax:dr]
	wav = zeros(len(r), dtype=double)
	SetRadialCoulombWave(Z, l, k, r, wav)

	return r, wav	


def GetRadialCoulombWaveBSplines(Z, l, k, bsplineObj, S):
	#Get the Coulomb function in grid space
	r = bsplineObj.GetQuadratureGridGlobal()
	wav = zeros(len(r), dtype=double)
	SetRadialCoulombWave(Z, l, k, r, wav)
	cplxWav = array(wav, dtype=complex)

	#get bspline coeffs
	coeff = zeros(bsplineObj.NumberOfBSplines, dtype=complex)
	bsplineObj.ExpandFunctionInBSplines(cplxWav, coeff)

	return coeff

def GetRadialBesselWaveBSplines(Z, l, k, bsplineObj, S):
	besselFunc = lambda k, r: sin(k*r + Z/k * log(k*r) - l*pi/2 + GetCoulombPhase(l, Z / k))

	#Get the Coulomb function in grid space
	r = bsplineObj.GetQuadratureGridGlobal()
	wav = array(besselFunc(k, r), dtype=complex)

	#get bspline coeffs
	coeff = zeros(bsplineObj.NumberOfBSplines, dtype=complex)
	bsplineObj.ExpandFunctionInBSplines(wav, coeff)

	return coeff


def CalculateDpDk(Z, k, psi, S):
	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)
	
	#Setup Radial Waves
	coulWaves = map(lambda curL: map(lambda curK: GetRadialCoulombWaveBSplines(Z, int(curL), curK, bspline, S), k), l)

	#calculate dpdk
	dpdk = zeros(len(k), dtype=double)
	for lIdx, (curL, wavList) in enumerate(zip(l, coulWaves)):
		for kIdx, (curK, wav) in enumerate(zip(k, wavList)):
			dpdk[kIdx] += abs(dot(conj(wav), dot(S, psi.GetData()[lIdx, :])))**2

	return dpdk


def PlotEnergyDistributionMorten():
	E0, dpde = numpy.loadtxt("dpde_morten.txt").transpose()
	plot(E0, dpde, label="morten")
	return E0, dpde
