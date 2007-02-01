#Create an instance of a DynamicPotentialEvaluator. Try to guess
#the python name of the class from a descriptive classname.
#className can be either the full name of the potential except the rank, 
#as specified in the corresponding pyste file.
#Optionally, classname can be just the name of the classname of the dynamic potential
#(without the "DynamicPotentialEvaluator_..."), if it was not renamed in the pyste
#file. 
def CreateDynamicPotentialInstance(className, rank):
	#Create globals
	glob = dict(sys.modules['__main__'].__dict__)
	glob.update(globals())
	
	potential = None
	longClassName = str("core.DynamicPotentialEvaluator") + "_" + className + "_" + str(rank) + "_" + str(rank)
	try:
		potential = eval(longClassName + "()", glob)
	except: pass
	
	if potential == None:
		shortClassName = className + "_" + str(rank)
		try:
			potential = eval(shortClassName + "()", glob)
		except: pass

	if potential == None:
		raise "Unknown potential", className 

	return potential


def CreateFiniteDifferencePotentialInstance(className, rank):
	#Create globals
	glob = dict(sys.modules['__main__'].__dict__)
	glob.update(globals())
	
	potential = None
	longClassName = "core.ExponentialFiniteDifferenceEvaluator" + "_" + className + "_" + str(rank) + "_" + str(rank)
	try:
		potential = eval(longClassName + "()", glob)
	except: pass
	
	if potential == None:
		shortClassName = "ExponentialFiniteDifferenceEvaluator" + "_" + className + "_" + str(rank)
		try:
			potential = eval(shortClassName + "()", glob)
		except: pass

	if potential == None:
		raise "Unknown potential", className 

	return potential

#Potential Wrapper interface
class PotentialWrapper:
	def ApplyConfigSection(self, configSection):
		self.ConfigSection = configSection
		self.Type = configSection.type
	
	def SetupStep(self, timestep):
		None
		
	def AdvanceStep(self, t, timestep):
		None

#Wrapper for a static potential. The potential is set up during
#SetupStep, and never changed. Timestep must be fixed during 
#propagation.
class StaticPotentialWrapper(PotentialWrapper):

	def __init__(self, psi):
		self.Potential = CreateInstanceRank("core.StaticPotential", psi.GetRank())
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
	
	def SetupStep(self, timeStep):
		#allocate memory
		self.Potential.InitializePotential(self.psi)
	
		if hasattr(self.ConfigSection, "function"):
			func = self.ConfigSection.function
			updateFunc = eval("core.SetPotentialFromGridFunction_" + str(self.psi.GetRank()))
			updateFunc(self.Potential, timeStep, self.psi, self.psi.GetRepresentation(), func, self.ConfigSection)

		elif hasattr(self.ConfigSection, "classname"):
			potentialEvaluator = CreateDynamicPotentialInstance(self.ConfigSection.classname, self.psi.GetRank())
			self.ConfigSection.Apply(potentialEvaluator)
			potentialEvaluator.UpdateStaticPotential(self.Potential, self.psi, timeStep, 0.0)
			self.PotentialEvaluator = potentialEvaluator
			print "dt = ", timeStep

		else:
			raise "Invalid potential config. Must specify either 'classname' or 'function'"
		
		
	def AdvanceStep(self, t, dt):
		self.Potential.ApplyPotential(self.psi)

	def GetPotential(self, dt):
		return real(log(self.Potential.GetPotentialData()) / (- 1.0j * dt))
		

#Wrapper for dynamic potentials. Dynamic potentials are classes implemented in 
#C++, and wrapped by DynamicPotentialEvaluator to make evaluation simple.
#Dynamic potentials have no large memory footprint, and can change each timestep,
#and also timeStep may change any time.
class DynamicPotentialWrapper(PotentialWrapper):

	def __init__(self, psi):
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()
	
		self.Potential = CreateDynamicPotentialInstance(configSection.classname, rank)
		configSection.Apply(self.Potential)
	
	def AdvanceStep(self, t, dt):
		self.Potential.ApplyPotential(self.psi, dt, t)

class FiniteDifferencePotentialWrapper(PotentialWrapper):
	def __init__(self, psi):
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()
	
		self.Evaluator = CreateFiniteDifferencePotentialInstance(configSection.classname, rank)
		configSection.Apply(self.Evaluator)
	
	def AdvanceStep(self, t, dt):
		self.Evaluator.UpdateWavefunction(self.psi, t, dt/2, 0)
		self.Evaluator.UpdateWavefunction(self.psi, t, dt, 1)
		self.Evaluator.UpdateWavefunction(self.psi, t, dt/2, 0)

