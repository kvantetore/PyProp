
class PolarPropagator(FourierPropagatorBase):
	__BASE = FourierPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

	def ApplyConfigSection(self, configSection):
		self.RadialRank = configSection.radial_rank

		self.__BASE.ApplyConfigSection(self, configSection)

	def CreateDefaultFourierPotential(self):
		fourierPotentialConf = self.StaticEnergyConf(PotentialType.Static, "core.PolarKineticPotential")
		fourierPotentialConf.mass = self.Mass
		fourierPotentialConf.angular_rank = self.TransformRank
		fourierPotentialConf.radial_rank = self.RadialRank
		fourierPot = CreatePotentialFromSection(fourierPotentialConf, "PolarKineticEnergy", self.psi)

		return fourierPot



