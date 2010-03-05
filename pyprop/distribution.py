import core

def CreateDistribution(config, rank=None):
	#Instance Distribution class which is templated over rank
	if rank == None:
		rank = config.Representation.rank
	distrib = CreateInstanceRank("core.DistributedModel", rank)
	
	#hack in initial distribution into configuration
	#TODO: Get initial distribution from propagator
	class sec:
		pass

	if hasattr(config, "Distribution"):
		distrSection = config.Distribution
	else:
		distrSection = sec()

	if not hasattr(distrSection, "proc_array_rank"):	
		distrSection.proc_array_rank = 1
	
	#DO NOT CHANGE THIS UNLESS YOU ARE ABSOLUTELY SURE. IF THE LAST RANK IS USED
	#SaveWavefunctionHDF WILL HAVE ABSOLUTELY HORRIBLE PERFORMANCE!!!
	if not hasattr(distrSection, "initial_distribution"):
		distrSection.initial_distribution = array([0], dtype=int)

	if distrSection.initial_distribution[0] != 0:
		print "WARNING: Not distributing first rank of wavefunction. Consider removing [Distribution] from the config file"

	#apply configuration
	distrib.ApplyConfigSection(distrSection)
	
	return distrib
	

def GetAnotherDistribution(distrib, rank):
	if rank <= len(distrib):
		raise Exception("Can not have distribution with length (%i) >= rank (%i)" % (rank, len(distrib)))

	ranks = [j for j in r_[0:rank] if j not in distrib]
	return array([min(ranks)], dtype=int)


def GetAnotherDistribution2(distrib, curRank, rank):
	if rank <= len(distrib):
		raise Exception("Can not have distribution with length (%i) >= rank (%i)" % (rank, len(distrib)))

	#find out which procRank curRank is distributed on
	curProcRank = find(array(distrib) == curRank)[0]

	#Ranks available for distribution
	availableRanks = [j for j in r_[0:rank] if j not in distrib]

	#The new distrib is similar to the old distrib, only with the distribution 
	#curProcRank changed
	newDistrib = array(distrib)
	newDistrib[curProcRank] = min(availableRanks)

	return newDistrib



