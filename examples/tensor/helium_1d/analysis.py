#------------------------------------------------------------------------
#                Symmetrization related functions
#------------------------------------------------------------------------
def GetSymmetrizedWavefunction(psi):
	"""
	Symmetrizes and anti-symmetrizes the wavefunction with respect to
	particle exchange.

	Returns a tuple of the symmetrized and anti-symmetrized wavefunction 
	(symPsi, antiSymPsi)
	"""
	exchgPsi = psi.Copy()
	exchgPsi.GetData()[:] = exchgPsi.GetData().transpose()

	#create symmetrized wavefunction
	symPsi = psi.Copy()
	symPsi.GetData()[:] += exchgPsi.GetData()
	symPsi.GetData()[:] *= 0.5
	
	antiSymPsi = exchgPsi
	antiSymPsi.GetData()[:] -= psi.GetData()
	antiSymPsi.GetData()[:] *= 0.5

	return symPsi, antiSymPsi


def SymmetrizeWavefunction(psi, symmetrize):
	if symmetrize:
		symfactor = 1
	else:
		symfactor = -1

	exchgPsi = psi.GetData().transpose()

	#create symmetrized wavefunction
	exchgPsi.GetData()[:] *= symfactor
	psi.GetData()[:] += exchgPsi.GetData()
	psi.GetData()[:] *= 0.5


#------------------------------------------------------------------------
#                        Setting up 1D eigenstates
#------------------------------------------------------------------------
def SetupEigenstates1D(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	1D Hydrogen system.

	eigenvalues is a list of 1-d eigenvalue arrays. Each array corresponding
	to a
	"""
	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)
	bspl = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	eigenValues = []
	eigenVectors = []

	eigenvectorScaling = 1

	M = SetupPotentialMatrix(prop, potentialIndices)

	E, V = scipy.linalg.eig(a=M, b=S)

	idx = argsort(real(E))
	E = real(E[idx])
	eigenValues = E

	#Sort and normalize eigenvectors
	sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
	eigenVectors = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()

	#assure correct phase convention (first oscillation should start out real positive)
	for i, curE in enumerate(E):
		bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
		phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
		eigenVectors[:,i] *= exp(-1.0j * phase)

	return eigenValues, eigenVectors


def SaveEigenstates1D(**args):
	postfix = GetGridPostfix(**args)
	outputFile = "eigenstates/eigenstates_1D_%s.h5" % ("_".join(postfix))

	#Setup problem
	prop = SetupProblem(**args)

	#Setp eigenvalues and eigenvectors
	eigenValues, eigenVectors = SetupEigenstates1D(prop)

	#Save eigenvalues and eigenvectors to file
	if outputFile != None:
		f = tables.openFile(outputFile, "w")
		try:
			#Create main group
			eigGroup = f.createGroup(f.root, "Eig")

			#save config object
			eigGroup._v_attrs.configObject = prop.Config.cfgObj

			#save eigenvalues and eigenstates
			f.createArray(eigGroup, "eigenvalues", eigenValues)
			f.createArray(eigGroup, "eigenvectors", eigenVectors)

		finally:
			f.close()


def SetupPotentialMatrix(prop, whichPotentials):
	matrixSize = prop.psi.GetData().shape[0]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		basisPairs = potential.BasisPairs[0]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[i]

	return matrix


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	matrix = SetupPotentialMatrix(prop, [overlap])
	return matrix


#------------------------------------------------------------------------
#        Load Single Particle States
#------------------------------------------------------------------------


def GetSingleStatesFile(**args):
	locations = ["eigenstates"]
	gridPostfix = "_".join(GetGridPostfix(**args))
	singleStatesFile = filter(os.path.exists, ["%s/eigenstates_1D_%s.h5" % (loc, gridPostfix) for loc in locations])
	if len(singleStatesFile) == 0:
		raise Exception("Could not find single-particle data files!")
	return singleStatesFile[0]


def LoadSingleParticleStates(singleStatesFile):
	f = tables.openFile(singleStatesFile, "r")

	#Create BSpline object
	cfgSection = pyprop.Section("BsplineRepresentation", f.root.Eig._v_attrs.configObject)
	bspl = pyprop.BSPLINE()
	bspl.ApplyConfigSection(cfgSection)

	#setup grid to check for correct phase convention
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	try:
		node = f.getNode("/Eig")
		E = node.eigenvalues[:]
		V = node.eigenvectors[:]

		#assure correct phase convention (first oscillation should start out real positive)
		#for i, curE in enumerate(E):
		#	bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
		#	phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
		#	V[:,i] *= exp(-1.0j * phase)

	finally:
		f.close()

	return E, V


def GetEnergyFilteredSingleParticleStates(stateFilter, **args):
	"""
	Returns the single particle states and energies corresponding to the 1D model and 
	an energy filter specified by stateFilter
	"""

	#load single particle states
	singleStatesFile = GetSingleStatesFile(**args)
	singleEnergies, singleStates = LoadSingleParticleStates(singleStatesFile)
	
	#filter states 
	filteredStates = [singleStates[:,i] for (i,E) in enumerate(singleEnergies) if stateFilter(E)]
	filteredEnergies = filter(stateFilter, singleEnergies)

	return filteredEnergies, filteredStates	


def SetupRadialCoulombStatesEnergyNormalized(psi, Z, Emax, dE, lmax):
	E = r_[dE:Emax:dE]
	k = sqrt(E*2)

	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)
	rcount = psi.GetRepresentation().GetRepresentation(1).GetFullShape()[0]
	
	#Setup Radial Waves
	states = []
	for l in range(lmax+1):
		V = zeros((rcount, len(k)), dtype=complex)
		for i, curk in enumerate(k):
			coeff = GetRadialCoulombWaveBSplines(Z, l, curk, bspline)
			V[:,i] = sqrt(2*dE/pi/curk) * coeff
		states.append(V)
	
	energies = [E]*(lmax+1)

	return energies, states


def GetSingleParticleCoulombStates(Z, dk, mink, maxk, lmax, radialRepr):
	"""
	Gets coulomb wave functions for every k between mink and maxk (in dk steps), 
	for every l up to (and including) lmax evaluated in bsplines.

	The structure returned is similar to that of LoadSingleParticleStates
	"""
	bspl = radialRepr.GetBSplineObject()
	k = r_[mink:maxk:dk]
	rcount = radialRepr.GetFullShape()[0]
	
	states = []
	for l in range(lmax+1):
		V = zeros((rcount, len(k)), dtype=complex)
		for i, curk in enumerate(k):
			coeff = GetRadialCoulombWaveBSplines(Z, l, curk, bspl)
			V[:,i] = coeff
		states.append(V)

	return [k]*(lmax+1), states


def GetRadialCoulombWaveBSplines(Z, l, k, bsplineObj):
	#Get the Coulomb function in grid space
	r = bsplineObj.GetQuadratureGridGlobal()
	wav = zeros(len(r), dtype=double)
	SetRadialCoulombWave(Z, l, k, r, wav)
	cplxWav = array(wav, dtype=complex)

	#get bspline coeffs
	coeff = zeros(bsplineObj.NumberOfBSplines, dtype=complex)
	bsplineObj.ExpandFunctionInBSplines(cplxWav, coeff)

	return coeff


def SetSymmetrizedProductState(psi, stateIdx1, stateIdx2, **args):
	"""
	Create a symmetrized 2D product state from two 1D eigenstates
	"""
	
	E, V = GetEnergyFilteredSingleParticleStates(lambda x: True, **args)

	singleState1 = V[stateIdx1]
	singleState2 = V[stateIdx2]

	psi.Clear()
	psi.GetData()[:] = outer(singleState1, singleState2)
	psi.GetData()[:] += outer(singleState2, singleState1)
	psi.Normalize()


#------------------------------------------------------------------------
#        Calculations on general product state combinations
#------------------------------------------------------------------------
def CalculatePopulationProductStates(V1, V2, psiData):
	numBsplines = psiData.shape[0]
	numStates1 = shape(V1)[0]
	numStates2 = shape(V2)[0]
	#print "numStates1 = %s, numStates2 = %s" % (numStates1, numStates2)
	tempData = zeros((numBsplines, numStates2), dtype=complex)
	populations = zeros((numStates1, numStates2), dtype=complex)

	#Project on V2 states
	#for i, v2 in enumerate(V2):
	tempData[:] = transpose(dot(conj(V2), psiData[:]))

	#Project on V1 states
	#for j, v1 in enumerate(V1):
	populations[:] = dot(conj(V1), tempData[:])

	#Get absolute square
	populations *= conj(populations)

	popList = []
	for i1 in range(numStates1):
		for i2 in range(numStates2):
			popList += [[i1, i2, 2 * real(populations[i1,i2])]]

	return popList
		

def CalculatePopulationProductStates2(V1, V2, psi):
	"""
	Calculate projection onto all possible symmetric product
	state combinations of V1 and V2 by explicitly forming the
	product states
	"""
	tmpPsi = psi.Copy()
	popList = []

	for i, v1 in enumerate(V1):
		for j, v2 in enumerate(V2):
			#print i, j
			#sys.stdout.flush()
			tmpPsi.GetData()[:] = outer(v1, v2)
			tmpPsi.GetData()[:] += outer(v2, v1)
			tmpPsi.Normalize()
			
			p = abs(psi.InnerProduct(tmpPsi))**2

			popList += [[i, j, p]]	

	return popList


def GetPopulationProductStates(psi, singleStates1, singleStates2):
	"""
	Calculates the population of psi in a set of single electron product states

	P_i = 2|< SingleState1_i(1), SingleState2_j(2) | psi(1,2) >|^2

	singleStates 1 and 2 are lists of 1D eigenstates created by SetupEigenstates
	
	the projection is carried out for every combination of singlestate1 and singlestate2
	"""

	#Make a copy of the wavefunction and multiply 
	#integration weights and overlap matrix
	tempPsi = psi.Copy()
	repr = psi.GetRepresentation()
	repr.MultiplyIntegrationWeights(tempPsi)
	data = tempPsi.GetData()

	#Get the population for every combination of v1 and v2
	popList = CalculatePopulationProductStates(singleStates1, singleStates2, data)

	return popList


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

	data = tempPsi.GetData()
	population = []

	for l1, V1 in enumerate(singleStates1):
		print "%i/%i" % (l1, len(singleStates1))

		for l2, V2 in enumerate(singleStates2):

			#filter out coupled spherical harmonic indices corresponding to this l
			lfilter = lambda coupledIndex: coupledIndex.l1 == l1 and coupledIndex.l2 == l2 
			angularIndices = array(GetLocalCoupledSphericalHarmonicIndices(psi, lfilter), dtype=int32)
			if len(angularIndices) == 0:
				continue
		
			#Remove projection for every combination of v1 and v2
			projV = RemoveProjectionRadialProductStates(l1, V1, l2, V2, data, angularIndices, psi.GetData())


	return population

