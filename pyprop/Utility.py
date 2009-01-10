
def GetFullWavefunctionData(prop):
	#allocate data
	fullData = ndarray(prop.psi.GetRepresentation().GetFullShape(), complex)
	fullDataFlat = ReshapeArray(fullData, [fullData.size])

	#get local data
	localData = prop.psi.GetData()
	if IsSingleProc():
		return localData

	localDataFlat = localData.flatten() 
	localSize = localData.size

	for curProc in range(0, ProcCount):
		#fullDataFlat[localSize*curProc:localSize*(curProc+1)] = pympi.bcast(localData, curProc)
		remoteData = pympi.bcast(localData, curProc)
		distribRank = prop.psi.GetRepresentation().GetDistributedModel().GetDistributedRank(prop.psi)
		start = curProc * localData.shape[distribRank]
		end = (curProc+1) * localData.shape[distribRank]
		if distribRank == 0:
			fullData[start:end, :] = remoteData
		elif distribRank == 1:
			fullData[:, start:end] = remoteData
		
	return fullData

def Linearize():
	for i in range(ProcCount):
		pympi.barrier()
		if pympi.rank == i:
			print "Process ", i, ": "
			yield i


def IsMaster():
	return ProcId == 0
	
def IsSingleProc():
	return ProcCount == 1
	
def ReshapeArray(array, newShape):
	return ndarray.__new__(array.__class__, dtype=array.dtype, shape=newShape, buffer=array.data)

def PrintOut(str=""):
	if IsMaster():
		print str
