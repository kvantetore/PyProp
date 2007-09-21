

#----------------------------------------------------------------------------------------------------
# Exponential Finite Difference Propagator
#----------------------------------------------------------------------------------------------------
class ExponentialFiniteDifferencePropagator(PropagatorBase):
	"""
	Propagates the wavefunction by using the exponential finite difference scheme. 
	The idea behind the scheme is to approximate the 1d laplace operator by finite difference
	This gives a tri-diagonal hamiltonian. This hamiltonian is then splitted into two block
	diagonal parts (A and B), where each block is a 2x2 matrix.
	We now use the split step scheme to approximate exp(A + B) by exp(A)exp(B) (or some higher
	order scheme). exp(A) and exp(B) can be found by diagonalizing each block for itself.
	This gives an algorithm with complexity O(n) which can be VERY efficiently parallelized.
	...

	The finite difference evaluation is implemented in a special potential evaluator
	found in core/finitediff/exponentialfinitedifference.h.

	An alternative to the exponential finite difference propagator, is the crank-nicholson finte
	difference scheme. To use this propagator, use a different potential evaulator: 
	core/finitediff/cranknicholson.h
	"""
	
	def __init__(self, psi):
		PropagatorBase.__init__(self, psi)

	def ApplyConfig(self, config):
		print config.Propagation.__dict__
		PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		PropagatorBase.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		#Set up potentials
		self.SetupPotential(dt)
		
	def AdvanceStep(self, t, dt):
		#Apply potential
		self.ApplyPotential(t, dt)
		#self.AdvanceStepKrylov(t, dt)

	def MultiplyHamiltonian(self, destPsi, t, dt):
		self.MultiplyPotential(destPsi, t, dt)


