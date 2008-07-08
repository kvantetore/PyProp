
class ArpackSolver:
	"""
	Pyprop wrapper for ARPACK (http://www.caam.rice.edu/software/ARPACK/)
	ARPACK is a package for finding a subset of the eigenvalues of a large
	(preferably sparse) matrix. 

	It works by using the matrix-vector product functionality of the Propagator

	"""

	def __init__(self, prop):
		self.BaseProblem = prop
		self.Rank = prop.psi.GetRank()

		configSection = prop.Config.Arpack
		silent = False
		if hasattr(prop.Config.Propagation, "silent"):
			silent = prop.Config.Propagation.silent 
		configSection.krylov_statistics = not silent

		#Create a copy of the wavefunction to calculate H|psi>
		if silent:
			Redirect.Enable(True)

		self.TempPsi = prop.psi.CopyDeep()

		if silent:
			Redirect.Disable()

		#Set up Arpack Solver
		self.Solver = CreateInstanceRank("core.krylov_ArpackPropagator", self.Rank)
		if not hasattr(configSection, "krylov_statistics"):
			silent = False
			if hasattr(prop.Config.Propagation, "silent"):
				silent = prop.Config.Propagation.silent 
			configSection.krylov_statistics = not silent
				
		prop.Config.Arpack.Apply(self.Solver)

		if not silent:
			memoryEstimate = self.Solver.GetMemoryEstimate(prop.psi)

		self.Solver.Setup(prop.psi)

	def Solve(self):
		psi = self.BaseProblem.psi;
		tempPsi = self.TempPsi

		#Run the Arnoldi iterations
		self.Solver.Solve(self.__MatVecCallback, psi, tempPsi)

	def __MatVecCallback(self, psi, tempPsi):
		tempPsi.GetData()[:] = 0
		self.BaseProblem.Propagator.MultiplyHamiltonian(tempPsi, 0, 0)

	def GetEigenvalues(self):
		"""
		Returns the real part of all the converged eigenvalues from arpack
	    """
		return self.Solver.GetEigenvalues().real

	def GetEigenvectors(self):
		"""
		Returns all eigenvectors as a N by M matrix, where 
		N is the number of converged eigenvalues, and M is the size of the
		wavefunction.
		The eigenvectors is normalized in the vector 2-norm, and can therefore not be
		expected to be normalized in the grid norm. Assign it to a wavefunction and call 
		psi.Normalize() to get a normalized eigenstate:

		eigenvectorIndex = 0
		eigenvectors = solver.GetEigenvectors()
		shape = psi.GetData().shape
		psi.GetData()[:] = numpy.reshape(eigenvectors[eigenvectorIndex, :])
		psi.Normalize()
		"""
		return self.Solver.GetEigenvectors()

	def SetEigenvector(self, psi, eigenvectorIndex, normalize=True):
		"""
		Sets psi to the eigenvector specified by eigenvetorIndex
		if normalize == True, psi will be normalized
		"""
		eigenvectors = self.GetEigenvectors()
		shape = psi.GetData().shape
		psi.GetData()[:] = numpy.reshape(eigenvectors[eigenvectorIndex, :], shape)
		if normalize:
			psi.Normalize()


