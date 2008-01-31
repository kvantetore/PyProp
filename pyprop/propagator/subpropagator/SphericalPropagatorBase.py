
class SphericalPropagatorBase:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		if transformRank != psi.GetRank() - 1:
			raise "SphericalTransform can only be used on the last rank"

	def ApplyConfigSection(self, config):
		self.Mass = config.mass
		self.RadialRank = config.radial_rank

	def SetupStep(self, dt):
		raise "Implement this in a inheriting class"

	def SetupStepConjugate(self, dt):
		#Transform from Grid to Spherical Harmonics
		self.Transform.ForwardTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)

	def AdvanceStep(self, t, dt):
		self.Potential.AdvanceStep(t, dt)
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)

	def AdvanceStepConjugate(self, t, dt):
		self.Transform.ForwardTransform(self.psi)
		self.Potential.AdvanceStep(t, dt)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		self.Potential.MultiplyPotential(dstPsi, t, dt)
		self.Transform.InverseTransform(self.psi)
		self.Transform.InverseTransform(dstPsi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)
		dstPsi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)

	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		self.Transform.ForwardTransform(self.psi)
		self.Transform.ForwardTransform(dstPsi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)
		dstPsi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationSphericalHarmonic)



