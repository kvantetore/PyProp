
class SphericalPropagator(SphericalPropagatorBase):
	BASE = SphericalPropagatorBase

	def __init__(self, psi, transformRank):
		self.BASE.__init__(self, psi, transformRank)

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
		self.Transform.SetupStep(self.psi)
		self.RepresentationSphericalHarmonic = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RepresentationTheta = self.Transform.CreateAngularRepresentation()
		self.RepresentationTheta.SetDistributedModel(self.RepresentationSphericalHarmonic.GetDistributedModel())

		#Transform from Spherical Harmonics to Grid
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)




