from numpy import zeros, double

from pyprop.propagator.subpropagator import SubPropagatorBase
from pyprop.config import Section
from pyprop.createinstance import CreateInstanceRank
import pyprop.potential as potential

import libreducedspherical


class ReducedSphericalPropagator(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		if transformRank != psi.GetRank() - 1:
			raise "SphericalTransform can only be used on the last rank"
		
		repr = psi.GetRepresentation().GetRepresentation(transformRank)
		self.MaxL = len(repr.Range.GetGrid()) - 1
		

	def ApplyConfigSection(self, config):
		self.Mass = config.mass
		self.RadialRank = config.radial_rank
		self.Config = config


	def SetupStep(self, dt):
		#Create config section for angular kinetic energy potential
		cfg = Section("AngularPotential")
		cfg.type = potential.PotentialType.Static
		cfg.classname = "libreducedspherical.ReducedAngularKineticEnergyPotential"
		cfg.mass = self.Mass
		cfg.l_rank = self.TransformRank
		cfg.radial_rank = self.RadialRank
		cfg.storage_model = potential.StaticStorageModel.StorageValue

		#Create potential
		if hasattr(self.Config, 'no_centrifugal_potential'):
			self.Potential = None
		else:
			self.Potential = potential.CreatePotentialFromSection(cfg, "ReducedAngularKineticEnergy", self.psi)
			self.Potential.SetupStep(dt)

		#Create transform
		self.Transform = CreateInstanceRank("libreducedspherical.ReducedSphericalTransform", self.psi.GetRank())
		self.Transform.SetupStep(self.psi, self.TransformRank)
		self.RepresentationSphericalHarmonic = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RepresentationTheta = self.Transform.CreateAngularRepresentation()
		self.RepresentationTheta.SetDistributedModel(self.RepresentationSphericalHarmonic.GetDistributedModel())

		#Transform from Spherical Harmonics to Grid
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationTheta)


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


	def GetBasisFunction(self, rank, basisIndex):
		basisFunction = zeros(self.MaxL+1, dtype=double)
		basisFunction[basisIndex] = 1.0
		return basisFunction



