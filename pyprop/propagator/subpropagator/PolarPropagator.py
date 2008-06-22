
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


	#HACK: Testing 
	def TransformInverse(self, psi):
		self.__BASE.TransformInverse(self, psi)


		"""
		A quick hack to antisymmetrize the wavefunction around 
		phi = pi/4 (z1 == z2) for 2 x 1D electrons
		"""
		"""
		phi = psi.GetRepresentation().GetLocalGrid(self.TransformRank)
		#psi(phi) = - psi(pi/2 - phi)
		mirrorAxis = where(phi == pi/2)[0][0]
		startIndex = where(phi == pi/4)[0][0]

		phiCount = psi.GetData().shape[self.TransformRank]
		curSlice = psi.GetRank() * [slice(None, None, None)]
		mirrorSlice = psi.GetRank() * [slice(None, None, None)]
		for i in range(phiCount/2+1):
			index = i + startIndex
			mirrorIndex = ((mirrorAxis - index) + phiCount) % phiCount
			if mirrorIndex == index:
				curSlice[self.TransformRank] = index
				psi.GetData()[tuple(curSlice)] = 0
			else:
				curSlice[self.TransformRank] = 	index
				mirrorSlice[self.TransformRank] = mirrorIndex
				psi.GetData()[tuple(mirrorSlice)] = - psi.GetData()[tuple(curSlice)]

		"""

