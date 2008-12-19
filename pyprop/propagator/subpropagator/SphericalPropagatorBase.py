
class SphericalPropagatorBase(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		if transformRank != psi.GetRank() - 1:
			raise "SphericalTransform can only be used on the last rank"

	def ApplyConfigSection(self, config):
		self.Mass = config.mass
		self.RadialRank = config.radial_rank
		self.Config = config

	def SetupStep(self, dt):
		raise "Implement this in a inheriting class"

	def SetupStepConjugate(self, dt):
		#Transform from Grid to Spherical Harmonics
		self.ForwardTransform(self.psi)

	def AdvanceStep(self, t, dt):
		if self.Potential != None:
			self.Potential.AdvanceStep(t, dt)

		self.InverseTransform(self.psi)

	def AdvanceStepConjugate(self, t, dt):
		self.ForwardTransform(self.psi)

		if self.Potential != None:
			self.Potential.AdvanceStep(t, dt)


	def MultiplyHamiltonian(self, srcPsi, dstPsi, t, dt):
		if self.Potential != None:
			self.Potential.MultiplyPotential(srcPsi, dstPsi, t, dt)

		self.InverseTransform(srcPsi)
		self.InverseTransform(dstPsi)

	def MultiplyHamiltonianConjugate(self, srcPsi, dstPsi, t, dt):
		self.ForwardTransform(srcPsi)
		self.ForwardTransform(dstPsi)

	def InverseTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.InverseTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)

	def ForwardTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.ForwardTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)


