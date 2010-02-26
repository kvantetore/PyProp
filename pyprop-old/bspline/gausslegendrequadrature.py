from numpy import *

class GaussQuadratureRule:
	"""
	Note: Only Gauss-Legendre quadrature available.
	"""

	def __init__(self, **args):

		self.RuleOrder = args['ruleOrder']
		self.Nodes = 0
		self.Weights = 0

		self.SetupJacobiMatrix()
		self.ComputeNodesAndWeights()


	def SetupJacobiMatrix(self):
		"""
		Set up (tridiagonal) Jacobi matrix for Legendre polynomials.
		In this case alpha = 0.
		"""

		beta = [1.0 / (2.0 * sqrt(1.0 - (2.0 * n)**-2)) for n in range(1, self.RuleOrder)] 
		self.JacobiMatrix = diag(array(beta), -1) + diag(array(beta), 1)


	def ComputeNodesAndWeights(self):
		"""
		Compute nodes and weights by finding eigenvalues and eigenvectors
		of the Jacobi matrix.
		"""
		
		# Compute nodes and weights from Jacobi matrix
		E, V = linalg.eig(self.JacobiMatrix)

		# Sort nodes in ascending order
		I = argsort(E)

		# Store nodes and weights
		self.Nodes = E[I]
		self.Weights = 2.0 * V[0,I]**2
