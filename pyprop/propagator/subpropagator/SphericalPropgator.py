
class SphericalPropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		if transformRank != psi.GetRank() - 1:
			raise "SphericalTransform can only be used on the last rank"

	def ApplyConfigSection(self, config):
		self.Mass = config.mass
		self.RadialRank = config.radial_rank

	def SetupStep(self, dt):
		#Create config section for angular kinetic energy potential
		cfg = Section("AngularPotential")
		cfg.type = PotentialType.Static
		cfg.classname = "core.AngularKineticEnergyPotential"
		cfg.mass = self.Mass
		cfg.l_rank = self.TransformRank
		cfg.radial_rank = self.RadialRank
		#Create potential
		self.Potential = CreatePotentialFromSection(cfg, "AngularKineticEnergy", self.psi)
		self.Potential.SetupStep(dt)

		#Create transform
		self.Transform = CreateInstanceRank("core.SphericalTransform", self.psi.GetRank())
		self.Transform.SetupStep(self.psi, self.TransformRank)
		self.RepresentationSphericalHarmonic = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RepresentationTheta = self.Transform.CreateAngularRepresentation()

		#Transform from Spherical Harmonics to Grid
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)


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



