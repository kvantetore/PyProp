#Author: Mads Lundeland Apr 2007
#Merged into pyprop main branch May 2. 2007 by torebi

class OrthoPolRadialPropagator(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		rank = psi.GetRank()
		self.Propagator = CreateInstanceRank("core.OrthoPolPropagator", rank)

	def ApplyConfigSection(self, configSection):
		configSection.Apply(self.Propagator)

	# Helper class for effective potential
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname

	def GetSubRepresentation(self, psi):
		return psi.GetRepresentation().GetRepresentation(self.TransformRank)

	def GetParam(self, psi):
		return self.GetSubRepresentation(psi).Range.Param

	def SetupAddedPotential(self, dt):
		if self.GetParam(self.psi).Type == PolynomialType.HermitePolynomial:
			potentialName = "AddedHermitePotential"
		elif self.GetParam(self.psi).Type == core.OrthoPolType.LaguerrePolynomial:
			potentialName = "AddedLaguerrePotential"
		else:
			raise "Invalid potential type"
		
		conf = self.StaticEnergyConf(PotentialType.RankOne, potentialName)
		conf.radial_rank = self.TransformRank
		conf.potential_rank = self.TransformRank
		conf.hyperspherical_dimensions = self.GetParam(self.psi).HypersphericalRank
		conf.scaling = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param.Scaling
		potential = CreatePotentialFromSection(conf, "AddedPotential", self.psi)
		potential.SetupStep(dt)
		self.AddedPotential = potential

	def SetupStep(self, dt):
		print "	Setup Potential"
		self.SetupAddedPotential(dt)

		print "	Setup Propagator"
		param = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param
		print "	param.Type = ", param.Type
		print "	param.Scaling = ", param.Scaling
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

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		self.Propagator.ApplyDifferentiationMatrix(self.psi, dstPsi)
		self.AddedPotential.MultiplyPotential(self.psi, dstPsi, t, dt)

	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass

	def GetBasisFunction(self, rank, basisIndex):
		return self.Propagator.GetEigenvectors()[basisIndex,:]






