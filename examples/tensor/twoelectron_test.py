import scipy.linalg

def SetupBigMatrix2D(prop, whichPotentials):
	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	
	#Allocate the hamilton matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" % (matrixSize, matrixSize, matrixSize**2 * 16 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype="complex")

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		
		#GetCoupledIndex = prop.psi.GetRepresentation().GetRepresentation(2).Range.GetCoupledIndex
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]

		for i, (x0,x0p) in enumerate(basisPairs0):
			for j, (x1,x1p) in enumerate(basisPairs1):
				indexLeft = (x1 * Count0) + x0
				indexRight = (x1p * Count0) + x0p
				BigMatrix[indexLeft, indexRight] += potential.PotentialData[i, j]

	return BigMatrix


def SetupBigMatrix(prop, whichPotentials):
	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	
	#Allocate the hamilton matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" % (matrixSize, matrixSize, matrixSize**2 * 16 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype="complex")

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs1 = potential.BasisPairs[0]
		basisPairs2 = potential.BasisPairs[1]
		basisPairs3 = potential.BasisPairs[2]
		
		#GetCoupledIndex = prop.psi.GetRepresentation().GetRepresentation(2).Range.GetCoupledIndex
		r1Count = prop.psi.GetData().shape[0]
		r2Count = prop.psi.GetData().shape[1]
		lCount = prop.psi.GetData().shape[2]

		for i, (A, Ap) in enumerate(basisPairs3):
			for j, (r1,r1p) in enumerate(basisPairs1):
				for k, (r2,r2p) in enumerate(basisPairs2):
					#indexLeft = (A*r1Count*r2Count) + (r1 * r2Count) + r2
					#indexRight = (Ap*r1Count*r2Count) + (r1p * r2Count) + r2p
					indexLeft = (r1 * r2Count*lCount) + (r2*lCount) + A
					indexRight = (r1p * r2Count*lCount) + (r2p*lCount) + Ap
					BigMatrix[indexLeft, indexRight] += potential.PotentialData[j, k, i]

	return BigMatrix


def CheckSymmetryFieldMatrix(whichPotential):
	prop = SetupProblem(configFile="config_helium.ini")
	BigMatrix = SetupBigMatrix(prop, whichPotential)
	
	#Check if matrix is hermittian
	hermDeviation = numpy.max(numpy.abs(BigMatrix - conj(BigMatrix.transpose())))
	print "Deviation from hermiticity (max norm) = %f" % hermDeviation

	return prop, BigMatrix



def GetTwoElectronEnergies(L=0, lmax=3):
	"""
	Get energies and eigenstates by direct diagonalization of L-subspace matrix
	"""
	#Coupled spherical index iterator based on given lmax and L
	index_iterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
	
	#Set up problem
	prop = SetupProblem(configFile="config_helium_local.ini", index_iterator=index_iterator)

	#Set up hamilton and overlap matrices
	HamiltonMatrix = SetupBigMatrix(prop, [0,1])
	OverlapMatrix = SetupBigMatrix(prop, [2])

	#Calculate generalized eigenvalues and eigenvectors
	print "Calculating generalized eigenvalues and eigenvectors..."
	sys.stdout.flush()
	E, V = scipy.linalg.eig(HamiltonMatrix, b=OverlapMatrix)

	return prop, HamiltonMatrix, OverlapMatrix, E, V
