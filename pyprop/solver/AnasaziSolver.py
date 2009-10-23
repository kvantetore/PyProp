import sys

from pyprop import CreateInstanceRank, ProcId
from pyprop.core import AnasaziSolver_1, AnasaziSolver_2, AnasaziSolver_3


class AnasaziSolver:
	"""
	Pyprop wrapper for Anasazi, a Trilinos eigensolver package.
	
	It works by using the matrix-vector product functionality of the Problem-object
	supplied in the constructor
	"""

	def __init__(self, prop, preconditioner = None):
		self.BaseProblem = prop
		self.Rank = prop.psi.GetRank()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = prop.psi.CopyDeep()

		#Set up Anasazi Solver
		self.Solver = CreateInstanceRank("AnasaziSolver", self.Rank, globals())
		configSection = prop.Config.Anasazi
		configSection.Apply(self.Solver)
		self.Solver.Setup(prop.psi)

		#Check for inverse iterations
		useInverseIterations = False
		if hasattr(configSection, "inverse_iterations"):
			useInverseIterations = configSection.inverse_iterations

		#Do we have a generalized eigenvalue problem?
		generalizedEigenvalueProblem = configSection.generalized_eigenvalue_problem
	
		self.Debug = False
		if hasattr(configSection, "krylov_debug"):
			if configSection.krylov_debug == True:	
				self.Debug = True

		self.CounterOn = False
		if hasattr(configSection, "counter_on"):
			if configSection.counter_on == True:
				self.CounterOn = True

		if preconditioner:
			self.Preconditioner = preconditioner
		elif hasattr(configSection, "preconditioner"):
			config = configSection.Config
			preconditionerName = configSection.preconditioner
			if preconditionerName:
				preconditionerSection = config.GetSection(preconditionerName)
				preconditioner = preconditionerSection.type(self.TempPsi)
				preconditionerSection.Apply(preconditioner)
				self.Preconditioner = preconditioner
				self.Preconditioner.SetHamiltonianScaling(1.0)
				self.Preconditioner.SetOverlapScaling(0.0)
				self.Preconditioner.Setup(self.BaseProblem.Propagator)

		eigenvalueCount = configSection.krylov_eigenvalue_count
		maxIterations = configSection.krylov_max_iteration_count
		self.TotalMaxIterations = 1
				
		#Check for user-defined matrix-vector product and generalized eigenvalue problem
		self.ApplyMatrix = self.BaseProblem.Propagator.MultiplyHamiltonian
		if generalizedEigenvalueProblem and not configSection.krylov_method == "KrylovSchur":
			self.ApplyMatrix = self.BaseProblem.Propagator.BasePropagator.MultiplyHamiltonianNoOverlap
		if hasattr(configSection, "matrix_vector_func"):
			if configSection.matrix_vector_func:
				self.ApplyMatrix = configSection.matrix_vector_func
		else:
			if useInverseIterations:
				raise Exception("Inverse Iterations must be used together with a userdefined matrix_vector_func")

	def Solve(self):
		psi = self.BaseProblem.psi;
		tempPsi = self.TempPsi

		self.Count = 0

		#Run the Arnoldi iterations
		if hasattr(self, "Preconditioner"):
			precCallback = self.__PreconditionCallback
		else:
			precCallback = None

		self.Solver.Solve(self.__MatVecCallback, precCallback, self.__OverlapCallback, psi, tempPsi)


	def __PreconditionCallback(self, srcPsi, dstPsi):
		dstPsi.GetData()[:] = srcPsi.GetData()
		self.Preconditioner.Solve(dstPsi)


	def __OverlapCallback(self, srcPsi, dstPsi):
		dstPsi.GetData()[:] = srcPsi.GetData()
		dstPsi.GetRepresentation().MultiplyOverlap(dstPsi)


	def __MatVecCallback(self, psi, tempPsi):
		self.ApplyMatrix(psi, tempPsi, 0, 0)
		self.Count += 1


	def GetEigenvalues(self):
		"""
		Returns the real part of all the converged eigenvalues from Anasazi
	    """
		return self.Solver.GetEigenvalues().real.copy()


	def GetEigenvector(self, index):
		"""
		Returns an eigenvector as a 1d numpy array.
		The eigenvectors are normalized in the vector 2-norm, and can therefore not be
		expected to be normalized in the grid norm. Assign it to a wavefunction and call 
		psi.Normalize() to get a normalized eigenstate:

		eigenvectorIndex = 0
		shape = psi.GetData().shape
		psi.GetData()[:] = numpy.reshape(solver.GetEigenvector(eigenvectorIndex)
		psi.Normalize()
		"""
		return self.Solver.GetEigenvector(index)


	def SetEigenvector(self, psi, eigenvectorIndex, normalize=True):
		"""
		Sets psi to the eigenvector specified by eigenvetorIndex
		if normalize == True, psi will be normalized
		"""
		shape = psi.GetData().shape
		psi.GetData()[:] = numpy.reshape(self.GetEigenvector(eigenvectorIndex), shape)
		if normalize:
			psi.Normalize()

