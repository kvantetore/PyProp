class BasisPropagatorEpetra(BasisPropagator):
	"""
	EpetraPropagator uses Trilinos::Epetra matrices to store potential matrices in
	a basis representation, and perform matrix-vector (potential-wavefunction) 
	multiplications with these. 
	
	Apart from this it is very similar to BasisPropagator, and uses TensorPotential
	to set up the potential matrices.
	"""

	__Base = PropagatorBase

	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()
		self.PotentialList = []
		self.CutOff = 0.0

		#Need a wavefunction for later
		self.Psi = self.psi.Copy()


	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)
		
		#Generate potentials
		self.GeneratePotentials(config)


	def GeneratePotentials(self, config):
		"""
		Generate TensorPotentials from potentials specified on the grid in the
		configuration file. Then consolidate/copy into Epetra matrices.
		"""

		#Potentials we should create on the fly
		if hasattr(config.Propagation, "grid_potential_list"):
			potentials = config.Propagation.grid_potential_list

			for potentialName in potentials:
				#Find the corresponding config section
				configSection = config.GetSection(potentialName)

				#generate potential 
				tensorPot = self.GeneratePotential(configSection)
				localBasisPairs = [geom.GetBasisPairs() for geom in tensorPot.GeometryList]

				#Setup a new Epetra potential
				epetraPot = EpetraPotential(self.Psi)
				epetraPot.ApplyConfigSection(configSection)

				#check if current potential can be consolidated with an existing one
				for existingPot in self.PotentialList:
					if existingPot.CanConsolidate(epetraPot):
						existingPot.AddTensorPotentialData(tensorPot.PotentialData, localBasisPairs, self.CutOff)
						existingPot.Name += "+" + tensorPot.Name
						tensorPot = None
						epetraPot = None
						break
						
				#add new potential to list
				if epetraPot != None:
					epetraPot.AddTensorPotentialData(tensorPot.PotentialData, localBasisPairs, self.CutOff)
					self.PotentialList.append(epetraPot)
					tensorPot = None

			#Global assemble all potentials
			for pot in self.PotentialList:
				pot.GlobalAssemble()

	
	def MultiplyHamiltonianBalancedOverlap(self, srcPsi, destPsi, t, dt):
		raise NotImplementedException("BasisPropagatorEpetra does not support MultiplyHamiltonianBalancedOverlap!")


