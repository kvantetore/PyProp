
#----------------------------------------------------------------------------------------------------
# Basic Propagator
#----------------------------------------------------------------------------------------------------
class PropagatorBase:
	"""
	This is the most basic propagator, and should probably be used as a basis for creating new
	propagators, but not used directly. It has functionality for setting up potential and applying 
	potentials, but no functionality for handling difference operators.

	"""
	
	def __init__(self, psi):
		#we'll need this later
		self.psi = psi

	def ApplyConfig(self, config): 
		#Create list of potentials.
		print "Creating Potentials..."
		potentialNames = config.Propagation.potential_evaluation
		self.PotentialList = [CreatePotential(config, name, self.psi) for name in potentialNames]
		
	def ApplyConfigSection(self, configSection):
		self.RenormalizeActive = configSection.renormalization
				
	#Propagator interface. Should most likely be overridden by inheriting classes
	def SetupStep(self, dt):
		print "    Setting up Potentials."
		self.SetupPotential(dt)
	
	def RestartPropagation(self, timestep, startTime, propagationTime):
		"""
		Implement on inheriting propagator if any particular restarting
		needs to be done.
		"""
		pass

	def AdvanceStep(self, t, dt):
		ApplyPotential(self.TimeStep)

		#check if we should renormalize
		if self.RenormalizeActive:
			self.psi.Normalize()

	def MultiplyHamiltonian(self, destPsi, t, dt):
		raise "MultiplyHamiltonian should be implemented by inheriting class %s, but it's not" % (self.__class__) 

	#Basic functionality that inheriting classes should use.
	def SetupPotential(self, dt):
		for potential in self.PotentialList:
			print "    Setting up potential ", potential.Name
			potential.SetupStep(dt)

	def ApplyPotential(self, t, dt):
		if dt == None:
			dt = self.TimeStep

		for potential in self.PotentialList:
			potential.AdvanceStep(t, dt)

	def MultiplyPotential(self, srcPsi, destPsi, t, dt):
		if dt == None:
			dt = self.TimeStep

		for potential in self.PotentialList:
			potential.MultiplyPotential(srcPsi, destPsi, t, dt)

	def GetPotentialList(self):
		raise NotImplementedError

	
	def PerformGridOperation(self, gridFunction):
		"""
		Perform a grid operation, as defined by the function 'gridFunction'.
		The wavefunction is assumed to be in the grid representation.
		"""
		gridFunction()

	def CalculatePotentialExpectationValue(self, tmpPsi, potential, t, dt):
		raise Exception("Must be implemented on subpropagator")
