 
def CreateInstanceRank(className, rank):
	try:
		return eval("%s_%i()" % (className, rank) )
	except Exception, ex:
		print "Could not create instance of class ", className, "with rank ", rank
		raise ex
		
def CreateDistribution(config, rank=None):
	#Currently we only support one distribution model
	if config.Distribution.model != "LargestStride":
		raise Exception, "Invalid DistributionModel '%s'" % config.Distribution.model

	#Instance Distribution class which is templated over rank
	if rank == None:
		rank = config.Representation.rank
	distrib = CreateInstanceRank("core.DistributedModel", rank)
	
	#Apply configuration section to distribution object 
	#(if it supports applying config data)
	config.Apply(distrib)
	config.Distribution.Apply(distrib)
	
	return distrib
	
	
def CreateRepresentation(config, distribution):
	#Create instance
	representation = config.Representation.type()

	#Set distribution model
	print "setting distributed model"
	representation.SetDistributedModel(distribution)
	
	#Apply configuration section
	config.Representation.Apply(representation)

	#If the representation is a base class of CombinedRepresentation,
	#we must set up the 1d sub-representations.
	combinedRepr = None
	try:
		combinedRepr = eval("core.CombinedRepresentation_" + str(config.Representation.rank))
	except:
		pass
	if combinedRepr != None and combinedRepr in config.Representation.type.__mro__:
		CreateSubRepresentations(representation, config)
	
	return representation

def CreateSubRepresentations(combinedRepr, config):
	rank = config.Representation.rank
	for i in range(rank):
		sectionName = config.Representation.Get("representation" + str(i))
		print "ConfigSection for rank %i is %s" % (i, sectionName)
		section = config.GetSection(sectionName)

		#create instance
		repr = section.type()
		repr.SetBaseRank(i)
		print "Representation for rank %i is %s" % (i, repr)

		#set distributed model
		#TODO: Fix sub-distributed models to reflect the same distribution
		#as the full distributed model.
		distrib = CreateDistribution(config, 1)
		repr.SetDistributedModel(distrib)

		#apply configuration
		section.Apply(repr)

		#Attach this repr to the main representation
		combinedRepr.SetRepresentation(i, repr)

	
def CreateWavefunction(config, representation):
	#Create instance
	print "    Creating instance"
	rank = config.Representation.rank
	psi = CreateInstanceRank("core.Wavefunction", rank)
	
	#Set reresentation
	print "    Setting representation"
	psi.SetRepresentation(representation)
	
	#Allocate data
	print "    Allocating data"
	psi.AllocateData()
	
	return psi
	
	
def CreatePropagator(config, psi):
	#Create instance
	propagator = None
	if hasattr(config.Propagation, "propagator"):
		propagator = config.Propagation.propagator(psi)
		config.Apply(propagator)
		config.Propagation.Apply(propagator)
	else:
		print "WARNING: No propagator specified in config file. Make sure your potential evaluator includes kinetic energy."
	
	return propagator
	
	
def CreatePotential(config, potentialName, psi):
	potentialConfig = config.GetSection(potentialName)
	return CreatePotentialFromSection(potentialConfig, potentialName, psi)

def CreatePotentialFromSection(potentialConfig, potentialName, psi):
	potential = None
	if potentialConfig.type == PotentialType.Static:
		potential = StaticPotentialWrapper(psi)

	elif potentialConfig.type == PotentialType.Dynamic:
		potential = DynamicPotentialWrapper(psi)
	
	elif potentialConfig.type == PotentialType.FiniteDifference:
		potential = FiniteDifferencePotentialWrapper(psi)
	
	else:
		raise "Unknown potential type", potentialConfig.type

	potentialConfig.Apply(potential)
	potential.Name = potentialName
	return potential

	
