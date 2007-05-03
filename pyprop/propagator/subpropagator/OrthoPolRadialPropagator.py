#Author: Mads Lundeland Apr 2007
#Merged into pyprop main branch May 2. 2007 by torebi

class OrthoPolRadialPropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		rank = psi.GetRank()
		self.Propagator = CreateInstanceRank("core.OrthoPolPropagator", rank)

	def ApplyConfigSection(self, configSection):
		configSection.Apply(self.Propagator)
		self.HypersphericalDimensions = configSection.hyperspherical_dimensions

	# Helper class for effective potential
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname

	def GetPolynomialType(self):
		return self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param.Type

	def SetupAddedPotential(self, dt):
		if self.GetPolynomialType() == core.OrthoPolType.HermiteTransform:
			potentialName = "AddedHermitePotential"
		elif self.GetPolynomialType() == core.OrthoPolType.LaguerreTransform:
			potentialName = "AddedLaguerrePotential"
		else:
			raise "Invalid potential type"
		
		conf = self.StaticEnergyConf(PotentialType.Static, potentialName)
		conf.radial_rank = self.TransformRank
		conf.hyperspherical_dimensions = self.HypersphericalDimensions
		conf.alpha = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Alpha
		potential = CreatePotentialFromSection(conf, "AddedPotential", self.psi)
		potential.SetupStep(dt)
		self.AddedPotential = potential

	def ApplyAddedPotential(self, t, dt):
		self.AddedPotential.AdvanceStep(t, dt)

	def SetupStep(self, dt):
		print "	Setup Potential"
		self.SetupAddedPotential(dt)

		print "	Setup Propagator"
		param = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param
		print "	param.Type = ", param.Type
		print "	param.Cutoff = ", param.Cutoff
		print "	param.HypersphericalRank = ", param.HypersphericalRank
		print "	TransformRank = ", self.TransformRank
		self.Propagator.Setup(param, dt, self.psi, self.TransformRank)

	def SetupStepConjugate(self, dt):
		pass

	def AdvanceStep(self, t, dt):
		self.AddedPotential.AdvanceStep(t, dt)
		self.Propagator.AdvanceStep(self.psi)

	def AdvanceStepConjugate(self, t, dt):
		self.Propagator.AdvanceStep(self.psi)
		self.AddedPotential.AdvanceStep(t, dt)





