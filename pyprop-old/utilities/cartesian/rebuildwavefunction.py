from numpy import meshgrid

class RankIndex:
	min = 0
	max = 1
	size = 2

class WavefunctionRebuilderCartesian:
	"""
	"""

	def __init__(self, prop):
		self.Propagator = prop

	def ApplyConfigSection(self, configSection):
		
		self.RankInfo = []
		self.RankInfo += [configSection.rank0]
		self.RankInfo += [configSection.rank1]

		self.GridExtent = (self.RankInfo[0][RankIndex.min], self.RankInfo[0][RankIndex.max], \
			self.RankInfo[1][RankIndex.min], self.RankInfo[1][RankIndex.max])

	def Setup(self):
		self.x = self.Propagator.psi.GetRepresentation().GetLocalGrid(0)
		self.y = self.Propagator.psi.GetRepresentation().GetLocalGrid(1)

		self.X, self.Y = meshgrid(self.x, self.y)

	def BuildWavefunction(self):
		"""
		"""
		psi = self.Propagator.psi.GetData()
		return self.X, self.Y, psi

