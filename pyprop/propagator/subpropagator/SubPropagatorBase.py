

class SubPropagatorBase:
	"""
	Base class for all subpropagators, that is all
	propagators that have responsibility for propagating
	a single rank
	"""

	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

	def ApplyConfigSection(self, configSection):
		pass

	def SetupStep(self, dt):
		raise Exception("SetupStep not implemented")

	def AdvanceStep(self, t, dt):
		raise Exception("AdvanceStep not implemented")

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		raise Exception("MultiplyHamiltonian not implemented")

	def SetupStepConjugate(self, dt):
		pass
	
	def AdvanceStepConjugate(self, t, dt):
		self.AdvanceStep(t, dt)
	
	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass


