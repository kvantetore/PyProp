import core
from createinstance import FindObjectStack
from distribution import PrintOut
	
def CreateRepresentation(config, distribution):
	rank = config.Representation.rank

	#Create instance of representation, if it is a string, find the
	#corresponding type
	reprType = config.Representation.type
	if type(reprType) == str:
		repr = None
		try:
			repr = FindObjectStack(reprType)
		except:
			repr = FindObjectStack("%s_%i" % (reprType, rank))
		reprType = repr
	representation = reprType()

	#Set distribution model
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
	if combinedRepr != None and combinedRepr in reprType.__mro__:
		CreateSubRepresentations(representation, config)
	
	return representation

def CreateSubRepresentations(combinedRepr, config):
	rank = config.Representation.rank
	for i in range(rank):
		sectionName = config.Representation.Get("representation" + str(i))
		section = config.GetSection(sectionName)

		subRank = getattr(section, "rank", 1)
		if subRank != 1:
			raise Exception("Only representations of rank 1 is supported as sub-representations")
		section.rank = subRank

		#create instance
		reprType = section.type
		if type(reprType) == str:
			r = None
			try:
				r = FindObjectStack(reprType)
			except:
				r = FindObjectStack("%s_%i" % (reprType, subRank))
			reprType = r
		repr = reprType()
		repr.SetBaseRank(i)
		PrintOut("Representation for rank %i is %s" % (i, repr))

		#set distributed model
		fullDistrib = combinedRepr.GetDistributedModel()
		distrib = fullDistrib.CreateSubDistributedModel()
		repr.SetDistributedModel(distrib)

		#apply configuration
		section.Apply(repr)

		#Attach this repr to the main representation
		combinedRepr.SetRepresentation(i, repr)


