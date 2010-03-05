
#----------------------------------------------------------------------------------------------------
# Propagator for Vector Representation
#----------------------------------------------------------------------------------------------------
class VectorPropagator(PropagatorBase):
	def __init__(self, psi):
		self.Psi = psi
		PropagatorBase.__init__(self, psi)

	def ApplyConfig(self, config): 
		PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		PropagatorBase.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		PropagatorBase.SetupStep(self, dt)

	def AdvanceStep(self, t, dt):
		raise NotImplementedException("Only MultiplyHamiltonian-base propagators are implemented for VectorPropagator")

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		self.MultiplyPotential(srcPsi, destPsi, t, dt)


