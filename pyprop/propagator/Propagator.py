
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

	def AdvanceStep(self, t, dt):
		ApplyPotential(self.TimeStep)

		#check if we should renormalize
		if self.RenormalizeActive:
			self.psi.Normalize()

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
			
			
	
