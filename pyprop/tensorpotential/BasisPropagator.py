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

		#We need the temp array for solving overlap matrix eqns
		self.TempPsi2 = self.psi.Copy()

	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)

		self.GeneratePotentials(config)
		self.ConsolidatePotentials()

	def GeneratePotential(self, configSection):
		#Create TensorPotential
		potential = TensorPotential(self.psi)
		configSection.Apply(potential)

		return potential

	def GeneratePotentials(self, config):
		"""
		Genereate TensorPotentials from potentials specified on the grid in the
		configuration file
		"""

		#Potentials we should create on the fly
		if hasattr(config.Propagation, "grid_potential_list"):
			potentials = config.Propagation.grid_potential_list

			for potentialName in potentials:
				#Find the corresponding config section
				configSection = config.GetSection(potentialName)
				#generate potential 
				pot = self.GeneratePotential(configSection)
				#check if this potential can be consolidated with an existing one
				for existingPot in self.PotentialList:
					if existingPot.CanConsolidate(pot):
						existingPot.PotentialData[:] += pot.PotentialData
						existingPot.Name += "+" + pot.Name
						pot = None
						break
						
				#add to potential list
				if pot != None:
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
				#Add otherPot to curPot if they can be consolidate
				if curPot.CanConsolidate(otherPot):
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



