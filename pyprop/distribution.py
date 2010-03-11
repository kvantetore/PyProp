import os

import core
from createinstance import CreateInstanceRank
from numpy import where, array, r_, ndarray


#Load mpi unless it is disabled
__DisableMPI = False
if 'PYPROP_SINGLEPROC' in os.environ:
	__DisableMPI = bool(os.environ['PYPROP_SINGLEPROC'])

ProcId = 0
ProcCount = 1
if __DisableMPI:
	print "MPI Disabled"
else:
	try:
		import pypar
		ProcId = pypar.rank()
		ProcCount = pypar.size()
	except:
		print "Warning: unable to load mpi."

if __DisableMPI:
	for name, obj in core.__dict__.iteritems():
		if name.startswith("DistributedModel_"):
			obj.ForceSingleProc()


# ------------- Convenience functions for printing out in parallel: -------------------

def Linearize(printProcId=False):
	"""
	Generator for performing a task on all processors serially, 
	usually for printing out stuff from all processors
	
	if run on 4 processors:
	>>> for i in Linearize():
	        print "Hello from proc %i" % ProcId
	Hello from proc 0
	Hello from proc 1
	Hello from proc 2
	Hello from proc 3
	"""
	for i in range(ProcCount):
		if ProcId == i:
			if printProcId:
				print "Process ", i, ": "
			yield i
		pypar.barrier()

def IsMaster():
	return ProcId == 0
	
def IsSingleProc():
	return ProcCount == 1
	
def ReshapeArray(array, newShape):
	return ndarray.__new__(array.__class__, dtype=array.dtype, shape=newShape, buffer=array.data)

def PrintOut(str=""):
	if IsMaster():
		print str


# ------------- Redistribution -------------------


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
	curProcRank = where(array(distrib) == curRank)[0][0]

	#Ranks available for distribution
	availableRanks = [j for j in r_[0:rank] if j not in distrib]

	#The new distrib is similar to the old distrib, only with the distribution 
	#curProcRank changed
	newDistrib = array(distrib)
	newDistrib[curProcRank] = min(availableRanks)

	return newDistrib



