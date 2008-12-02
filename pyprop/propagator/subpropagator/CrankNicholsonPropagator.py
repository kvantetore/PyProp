class CrankNicholsonPropagator(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		self.Propagator = CreateInstanceRank("core.CrankNicholsonPropagator", psi.GetRank())

	def ApplyConfigSection(self, configSection):
		configSection.transform_rank = self.TransformRank
		configSection.Apply(self.Propagator)
		self.ConfigSection = configSection

	def SetupStep(self, dt):
		self.Propagator.SetupStep(self.psi, dt)

	def AdvanceStep(self, t, dt):
		self.Propagator.AdvanceStep(self.psi, dt)

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		self.Propagator.MultiplyKineticEnergyOperator(self.psi, dstPsi)

	def SetupStepConjugate(self, dt):
		pass		
		
	def AdvanceStepConjugate(self, t, dt):
		self.AdvanceStep(t, dt)
	
	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass

	def SupportsParallelPropagation(self):
		return True
	
	def ForwardTransform(self, psi):
		pass

	def InverseTransform(self, psi):
		pass

