
class SphericalPropagatorBase:

	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

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


	def MultiplyHamiltonian(self, dstPsi, t, dt):
		if self.Potential != None:
			self.Potential.MultiplyPotential(self.psi, dstPsi, t, dt)

		self.InverseTransform(self.psi)
		self.InverseTransform(dstPsi)

	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		self.ForwardTransform(self.psi)
		self.ForwardTransform(dstPsi)

	def InverseTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.InverseTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)

	def ForwardTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.ForwardTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)


