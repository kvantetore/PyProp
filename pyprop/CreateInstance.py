 
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
	
	if config.Representation.type == core.SphericalRepresentation3D:
		CreateSubRepresentations(representation, config)
	
	#Set distribution model
	representation.SetDistributedModel(distribution)
	
	#Apply configuration section
	#If representation is SphericalRepresentation3D object don't need
	# to apply config since it was applied in CreateSubRepresentations
	if config.Representation.type != core.SphericalRepresentation3D:
		config.Representation.Apply(representation)

	return representation

def CreateSubRepresentations(representation, config):
	#1D distributed models
	distrib1d1 = CreateDistribution(config, 1)
	distrib1d2 = CreateDistribution(config, 1)

	print "Creating RadialRepresentation..."
	#radial part
	radialRepr = config.Representation.radialtype()
	radialRepr.SetBaseRank(0)
	print "Radial Repr: %s" % radialRepr
	radialRepr.SetDistributedModel(distrib1d1)
	config.RadialRepresentation.Apply(radialRepr)
	representation.SetRepresentation(0,radialRepr)
	
	print "Creating AngularRepresentation..."
	#angular part
	angularRepr = config.Representation.angulartype()
	angularRepr.SetBaseRank(1)
	angularRepr.SetDistributedModel(distrib1d2)
	config.AngularRepresentation.Apply(angularRepr)
	representation.SetRepresentation(1,angularRepr)

	print "done creating AngularRepresentation"
	
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
		#why was this here
		#config.Representation.Apply(propagator)
	else:
		print "WARNING: No momentum evaluator specified in config file. Make sure your potential evaluator includes kinetic energy."
	
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

	
