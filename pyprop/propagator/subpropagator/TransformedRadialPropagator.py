
class TransformedRadialPropagator(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		rank = psi.GetRank()	
		self.Propagator = CreateInstanceRank("core.TransformedGridPropagator", rank)
		self.TransformRank = transformRank

	def ApplyConfigSection(self, configSection):
		configSection.Apply(self.Propagator)

	def SetupStep(self, dt):
		param = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param
		
		self.Propagator.Setup(param, dt, self.psi, self.TransformRank)

	def AdvanceStep(self, t, dt):
		self.Propagator.AdvanceStep(self.psi)

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		self.Propagator.ApplyDifferentiationMatrix(self.psi, dstPsi)
	
	def SetupStepConjugate(self, dt):
		pass

	def AdvanceStepConjugate(self, t, dt):
		self.AdvanceStep(t, dt)
	
	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass

	def GetBasisFunction(self, rank, basisIndex):
		return self.Propagator.GetEigenvectors()[basisIndex,:]


