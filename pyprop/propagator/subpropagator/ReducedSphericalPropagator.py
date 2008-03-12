
class ReducedSphericalPropagator(SphericalPropagatorBase):
	BASE = SphericalPropagatorBase

	def __init__(self, psi, transformRank):
		self.BASE.__init__(self, psi, transformRank)

		repr = psi.GetRepresentation().GetRepresentation(transformRank)
		self.MaxL = len(repr.Range.GetGrid()) - 1

	
	def SetupStep(self, dt):
		#Create config section for angular kinetic energy potential
		cfg = Section("AngularPotential")
		cfg.type = PotentialType.Static
		cfg.classname = "core.ReducedAngularKineticEnergyPotential"
		cfg.mass = self.Mass
		cfg.l_rank = self.TransformRank
		cfg.radial_rank = self.RadialRank

		#Create potential
		if hasattr(self.Config, 'no_centrifugal_potential'):
			self.Potential = None
		else:
			self.Potential = CreatePotentialFromSection(cfg, "ReducedAngularKineticEnergy", self.psi)
			self.Potential.SetupStep(dt)

		#Create transform
		self.Transform = CreateInstanceRank("core.ReducedSphericalTransform", self.psi.GetRank())
		self.Transform.SetupStep(self.psi, self.TransformRank)
		self.RepresentationSphericalHarmonic = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RepresentationTheta = self.Transform.CreateAngularRepresentation()
		self.RepresentationTheta.SetDistributedModel(self.RepresentationSphericalHarmonic.GetDistributedModel())

		#Transform from Spherical Harmonics to Grid
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)


	def GetBasisFunction(self, rank, basisIndex):
		basisFunction = zeros(self.MaxL+1, dtype=double)
		basisFunction[basisIndex] = 1.0
		return basisFunction

