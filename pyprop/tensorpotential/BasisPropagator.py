class BasisPropagator(PropagatorBase):
	"""
	Propagator that does not transform the wavefunction to the grid basis, but rather applies the potentials
	(TensorPotentials) directly in the basis.

	BasisPropagator assumes that the wavefunction has a CombinedRepresentation.

	BasisPropagator has facilities to use TensorPotentialGenerator to generate basis-representations of potentials,
	and is mostly used to propagate TensorPotential-based problems

	AdvanceStep is not supported by BasisPropagator. Rather, one should use a propagator utilizing MultplyHamiltonian
	such as the ODE and Krylov propagator.
	"""

	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

		#Create a tensor potential generator 		
		#We will use this to construct tensor potentials from potentials specified on the grid
		self.TensorPotentialGenerator = TensorPotentialGenerator(representation = psi.GetRepresentation())

		#We need the temp array for solving overlap matrix eqns
		self.TempPsi2 = self.psi.Copy()

	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)

		self.GeneratePotentials(config)
		self.ConsolidatePotentials()

	def GeneratePotential(self, configSection):
		generator = self.TensorPotentialGenerator

		#Use TensorPotentialGenerator to construct potential in basis
		geometryList = generator.GetGeometryList(configSection)
		potentialData = generator.GeneratePotential(configSection)
		originalPotential = getattr(generator, "OriginalPotential", None)

		#Create PotentialWrapper for TensorPotential
		potential = TensorPotential(self.psi)
		configSection.Apply(potential)
		potential.GeometryList = geometryList
		potential.PotentialData = potentialData
		potential.OriginalPotential = originalPotential
		potential.Name = configSection.name

		return potential

	def GeneratePotentials(self, config):
		"""
		Genereate TensorPotentials from potentials specified on the grid in the
		configuration file
		"""

		#Potentials we should create on the fly
		if hasattr(config.Propagation, "grid_potential_list"):
			potentials = config.Propagation.grid_potential_list
			generator = self.TensorPotentialGenerator

			for potentialName in potentials:
				#Find the corresponding config section
				configSection = config.GetSection(potentialName)
				#generate potential 
				pot = self.GeneratePotential(configSection)
				#add to potential list
				self.PotentialList.append(pot)

	def ConsolidatePotentials(self):
		"""
		Try to consolidate potentials having the same geometry into
		one potential to save evaulations. Potentials containing a time
		dependent part needs to be treated separately, and thus is not considered
		"""
		
		PrintOut( "Consolidating similar potentials: " )
		PrintOut( "Starting with potentials:" )
		for pot in self.PotentialList:
			PrintOut( "    %s" % pot.Name )

		#only non timedependent potentials are considered
		potentials = list(self.PotentialList) #[pot for pot in self.PotentialList if not pot.IsTimeDependent]
		removePotentials = []

		#Consolidate til we're at the last potential
		i = 0
		while(i < len(potentials)-1):
			curPot = potentials[i]

			#Loop over all potentials after curPot
			for otherPot in list(potentials[i+1:]):
				#We can consolidate curPot and otherPot if all the index pairs are the same,
				#And curPot and otherPot has the same time dependency
				canConsolidate = True

				#If one of the potentials is time dependent the other must also be
				if canConsolidate and curPot.IsTimeDependent != otherPot.IsTimeDependent:
					canConsolidate = False

				#If they are both time dependent, they must have the same time function
				if canConsolidate and curPot.IsTimeDependent and otherPot.IsTimeDependent:
					if otherPot.OriginalTimeFunction != curPot.OriginalTimeFunction:
						canConsolidate = False

				#don't consolidate debug potentials
				if canConsolidate and (curPot.DebugPotential or otherPot.DebugPotential):
					canConsolidate = False
					
				#Both potentials must have the same storage 
				if canConsolidate:
					for rank in range(self.Rank):
						if curPot.GeometryList[rank].GetStorageId() != otherPot.GeometryList[rank].GetStorageId():
							canConsolidate = False
							break

				if canConsolidate:
					for rank in range(self.Rank):
						if any(curPot.GeometryList[rank].GetBasisPairs() != otherPot.GeometryList[rank].GetBasisPairs()):
							canConsolidate = False
							break
					
				#Add otherPot to curPot
				if canConsolidate:
					curPot.PotentialData[:] += otherPot.PotentialData
					curPot.Name += "+" + otherPot.Name
					potentials.remove(otherPot)
					removePotentials.append(otherPot)

			#Next potential
			i += 1

		#Remove consolidated potentials from the main potential list
		for pot in removePotentials:
			self.PotentialList.remove(pot)

		PrintOut( "Ended up with potentials:" )
		for pot in self.PotentialList:
			PrintOut( "    %s" % pot.Name )


	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		self.__Base.SetupStep(self, dt)

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		#Multiply potentials
		destPsi.GetData()[:] = 0
		self.MultiplyPotential(srcPsi, destPsi, t, dt)

		#Solve for all overlap matrices
		repr = srcPsi.GetRepresentation()
		if not all([repr.IsOrthogonalBasis(i) for i in range(self.Rank)]):
			repr.SolveOverlap(destPsi)

	def MultiplyHamiltonianNoOverlap(self, srcPsi, destPsi, t, dt):
		#Multiply potentials
		destPsi.GetData()[:] = 0
		self.MultiplyPotential(srcPsi, destPsi, t, dt)

	def MultiplyHamiltonianBalancedOverlap(self, srcPsi, destPsi, t, dt):
		#Store input psi
		def StorePsi():
			self.TempPsi2.GetData()[:] = srcPsi.GetData()

		#Solve for all overlap matrices
		def SolveOverlap1():
			repr = srcPsi.GetRepresentation()
			repr.SolveSqrtOverlap(False, srcPsi)
		
		#Multiply potentials
		def MultiplyPotential():
			destPsi.GetData()[:] = 0
			self.MultiplyPotential(srcPsi, destPsi, t, dt)

		#Solve for all overlap matrices
		def SolveOverlap2():
			repr = destPsi.GetRepresentation()
			repr.SolveSqrtOverlap(True, destPsi)

		#Restore input psi back to its original state
		def RestorePsi():
			srcPsi.GetData()[:] = self.TempPsi2.GetData()

		StorePsi()
		SolveOverlap1()
		MultiplyPotential()
		SolveOverlap2()
		RestorePsi()

	def AdvanceStep(self, t, dt):
		raise NotImplementedException("BasisPropagator does not support AdvanceStep. Use it as a base for explicit propagators")

	def GetBasisFunction(self, rank, basisIndex):
		raise NotImplementedException("Implement GetBasisFunction...")



