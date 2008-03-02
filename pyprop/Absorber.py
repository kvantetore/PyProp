

"""

Example config.ini sections

[AbsorbingPotential]
type = PotentialType.Static
storage_model = StaticStorageModel.StorageExpValue
classname = "CombinedAbsorber"
absorbers = ["RadialAbsorber"]

[RadialAbsorber]
type = AbsorbingBoundary
rank = 0
width = 10
absorb_left = True
absorb_right = True

"""


class AbsorbingBoundary(core.AbsorberModel):

	def ApplyConfigSection(self, config):
		self.Width = config.width
		self.AbsorbLeft = config.absorb_left
		self.AbsorbRight = config.absorb_right

	def SetupStep(self, grid):
		minPos = min(grid)
		maxPos = max(grid)

		self.Scaling = ones(len(grid), dtype=double)
		for i in range(len(grid)):	
			pos = grid[i]
			width = self.Width
			if pos < minPos + width and self.AbsorbLeft:
				self.Scaling[i] *= cos(pi/2.0 * (pos - (minPos + width)) / width)**(1./8.);
			if pos > maxPos - width and self.AbsorbRight:
				self.Scaling[i] *= cos(pi/2.0 * (- maxPos + pos + width) / width)**(1./8.);

	def	GetScaling(self):
		return self.Scaling


class CombinedAbsorber:
	"""
	Combined absorber is an absorber that can contain a number of sub-absorbers,
	where each perform in a specific direction.

	"""
	
	def ApplyConfigSection(self, configSection):
		config = configSection.Config

		self.AbsorberList = []
		absorberList = configSection.absorbers
		for absorberName in absorberList:
			#Create absorber
			absorberConfig = config.GetSection(absorberName)
			absorber = absorberConfig.type()
			#Apply configsection
			absorberConfig.Apply(absorber)
			#Add to absorberlist with rank
			self.AbsorberList.append( (absorber, absorberConfig.rank) )
			
	def SetupStep(self, psi, dt):
		rank = psi.GetRank()
		self.AbsorberPotential = CreateInstanceRank("core.CombinedAbsorberPotential", rank)

		repr = psi.GetRepresentation()
		for absorber, rank in self.AbsorberList:
			grid = repr.GetLocalGrid(rank)
			absorber.SetupStep(grid)
			self.AbsorberPotential.AddAbsorber(absorber, rank)


	def MultiplyPotential(self, psi, destPsi, dt, t):
		"""
		MultiplyPotential for absorber does not do anything.
		"""
		pass
	
	def CalculateExpectationValue(self, psi, dt, t):
		return 0

	def ApplyPotential(self, psi, dt, t):
		self.AbsorberPotential.ApplyPotential(psi, dt, t)
		
			

