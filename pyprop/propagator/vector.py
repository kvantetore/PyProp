import pyprop.propagator.base as base
from pyprop.exceptions import NotImplementedException

#----------------------------------------------------------------------------------------------------
# Propagator for Vector Representation
#----------------------------------------------------------------------------------------------------
class VectorPropagator(base.PropagatorBase):
	def __init__(self, psi):
		self.Psi = psi
		base.PropagatorBase.__init__(self, psi)

	def ApplyConfig(self, config): 
		base.PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		base.PropagatorBase.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		base.PropagatorBase.SetupStep(self, dt)

	def AdvanceStep(self, t, dt):
		raise NotImplementedException("Only MultiplyHamiltonian-base propagators are implemented for VectorPropagator")

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		self.MultiplyPotential(srcPsi, destPsi, t, dt)


