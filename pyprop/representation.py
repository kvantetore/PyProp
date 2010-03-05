
	
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
		fullDistrib = combinedRepr.GetDistributedModel()
		distrib = fullDistrib.CreateSubDistributedModel()
		repr.SetDistributedModel(distrib)

		#apply configuration
		section.Apply(repr)

		#Attach this repr to the main representation
		combinedRepr.SetRepresentation(i, repr)


