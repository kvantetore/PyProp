try:
	import scipy.linalg
except:
	pass
import tables

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

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		basisPairs2 = potential.BasisPairs[2]
		
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]
		Count2 = prop.psi.GetData().shape[2]

		for i, (x0,x0p) in enumerate(basisPairs0):
			for j, (x1,x1p) in enumerate(basisPairs1):
				for k, (x2,x2p) in enumerate(basisPairs2):
					indexLeft = x2 + (x1 * Count2) + (x0 * Count1 * Count2) 
					indexRight = x2p + (x1p * Count2) + (x0p * Count1 * Count2) 
					BigMatrix[indexLeft, indexRight] += potential.PotentialData[i, j, k]

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



def TestFieldCoupling(config="config_2e_fieldtest.ini", outputs=100):
	"""
	Testing the two-electron electric field velocity coupling. The electron-electron
	interaction is neglected, thus we may compare with known one-electron results.
	"""

	solver = FindEigenvalues(config=config)

	prop = SetupProblem(silent=False, imtime=False, additionalPotentials = \
		["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", \
		"LaserPotentialVelocity"], config=config)
	#prop = SetupProblem(silent=False, imtime=False, configFile=config)
	
	prop.psi.Clear()
	solver.SetEigenvector(prop.psi, 0)
	prop.psi.Normalize()
	initPsi = prop.psi.Copy()

	timeList = []
	corrList = []
	normList = []
	for t in prop.Advance(outputs):
		n = prop.psi.GetNorm()
		if n != n:
			return prop
		c = abs(prop.psi.InnerProduct(initPsi))**2
		timeList.append(t)
		normList.append(n)
		corrList.append(c)
		if pyprop.ProcId == 0:
			print "t = %.14f, N(t) = %.14f, P(t) = %.14f" % (t, n, c)
	
	prop.Propagator.PampWrapper.PrintStatistics()
		
	c = abs(prop.psi.InnerProduct(initPsi))**2
	print "Final Correlation = %f" % c

	outFile = "twoelectron_field_test.h5"
	prop.SaveWavefunctionHDF(outFile, "/wavefunction")
	if pyprop.ProcId == 0:
		h5file = tables.openFile(outFile, "r+")
		try:
			h5file.createArray("/", "SampleTimes", timeList)
			h5file.createArray("/", "Norm", normList)
			h5file.createArray("/", "InitialCorrelation", corrList)
		finally:
			h5file.close()

