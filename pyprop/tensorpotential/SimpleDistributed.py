"""
Global Matrix is N-by-N sparse matrix 

The matrix divided on ProcCount processors

Matrix is stored in packed format, where the current 
processors has stored all rows i of the matrix, and only 
cols j, where startCol_p <= j < startCol_p+colCount_p

Assume N % ProcCount = 0

colCount_p = N / ProcCount
startCol_p = procRank * (N / ProcCount)

Variables without prefix are the same on every processor
Variables prefixed with local are indexed relative to the starting point of the current processor
Variables prefixed with global are indexed absolutely in the global arrays

A      = [ a00 a01 a02     a04     ]
         [ a10 a11 a12 a13 a14 a15 ]
		 [ a20 a21 a22 a23 a24     ]
		 [     a31 a32 a33 a34 a35 ]
		 [ a40 a41 a42 a43 a44 a45 ]
		 [     a51     a53 a54 a55 ]


Matrix storage is packed, with an additiona index array specifying the row and col of each element.
There is made no attempt to balance the matrix, and siginificant load imbalancing may occur. this 
should be checked for, by assuring that the length of each localMatrix is more or less the same

localMatrix0 = [ a00 a01 a10 a11 a20 a20 a31 a40 a41 a51]
localMatrix1 = [ a02 a12 a13 a22 a23 a32 a33 a42 a43 a53]
localMatrix2 = [ a04 a14 a15 a24 a34 a35 a44 a45 a54 a55]

localIndices0 = [ (0,0), (0,1), ..., (5,1) ]
localIndices1 = [ (0,2), (1,2), ..., (5,3) ]
localIndices3 = [ (0,4), (1,4), ..., (5,5) ]

(normally the distribution will not be completely balanced like here)
and the different localMatrices will be of different length.

both input and output psi is distributed the same way
            proc0      proc1     proc2
psi =    [  c0  c1 |  c2  c3 |  c4  c5 ]

The algorithm is then that each processor iterates over 
all the rows it has, calculates the product for that row
and sends it to the correct processor.

In order to ensure better load balancing, each processor should
sort the indices by row then col, and then start with the first element
that does not require sending.

the length of the largest localMatrix is the number of "steps" needed to complete 
the matvec multiply.

each processor maintains a list for each step which proc to send data to,
which procs (may be more than one) to recieve data from, and which matrix element to
process.


indexPairs            global list of all index pairs
distribIndexList      indices into indexPairs for each processor
globalSize            size of the matrix (size*size is the number of elements in a dense matrix)
rank                  which tensor rank of the wavefunction this matrix is working on (should be 0)


"""

def GetGlobalStartIndex(size, procCount, procId):
	rest = size % procCount
	if rest != 0:
		size += procCount - rest
	return (size / procCount) * procId

def GetDistributedShape(size, procCount, procId):
	rest = size % procCount
	if rest == 0:
		distrShape = size / procCount;
	else:
		paddedDistrShape = (size + procCount - rest) / procCount 
		shape = size - paddedDistrShape * procId
		shape = max(shape, 0)
		shape = min(shape, paddedDistrShape)
		distrShape = shape
	return distrShape;


class MultiplyStepInfo(object):
	"""
	LocalMatrixIndex	index into local matrix for this step is step-number or -1 to signal empty
	GlobalRow			row to send this result to
	GlobalCol		 	col to get source data from

	SendProc			procId to send the result of this operation to
	RecvProcList		procs to recieve from this
	RecvLocalRowList	local indices into the dest vector where to store the result from the recvs
	"""

	def __init__(self, localMatrixIndex, globalRow, globalCol, sendProc, recvProcList, recvLocalRowList):
		self.LocalMatrixIndex = localMatrixIndex
		self.GlobalRow = globalRow
		self.GlobalCol = globalCol
		self.SendProc = sendProc
		self.RecvProcList = recvProcList
		self.RecvLocalRowList = recvLocalRowList

	def __repr__(self):
		return "idx pair = [%i, %i], recvProcList = %s" % (self.GlobalRow, self.GlobalCol, self.RecvProcList)

def SetupDistributedIndexList(globalSize, indexPairs, distrib, rank):
	"""
	Distributes index pairs on the different processors. 

	indexPairs is a two dimensional array containing a list 
	of [row, col] pairs. the row indices are in indexPairs[:,0]
	and the col indices are in indexPairs[:,1]

	distrib is the distributed model knowing how the wavefunction
	is distributed
	
	A list of arrays is returned, each array contains indices
	into indexPairs for that processor.
	"""
	procCount = distrib.ProcCount
	procId = distrib.ProcId

	distribIndexList = [[] for i in range(procCount)]

	colStart = array([distrib.GetLocalStartIndex(int(globalSize), int(rank), int(curProc)) for curProc in range(procCount)], dtype=int)
	colEnd = array(list(colStart[1:]) + [globalSize])

	#Iterate through all index pairs and figure out on which processor
	#each index pair belongs
	for i, (row, col) in enumerate(zip(indexPairs[:,0], indexPairs[:,1])):
		foundProc = None
		for curProc, (curStart, curEnd) in enumerate(zip(colStart, colEnd)):
			if curStart <= col < curEnd:
				foundProc = curProc
				break
		if foundProc==None:
			raise Exception("what!?")

		distribIndexList[foundProc].append(i)

	#For each processor, start with local part of the matrix. 
	#this gives a better load balancing
	for curProc in range(procCount):
		for listIndex, pairIndex in enumerate(distribIndexList[curProc]):
			row, col = indexPairs[pairIndex, 0], indexPairs[pairIndex, 1]
			rowStart = colStart[curProc]
			rowEnd = colEnd[curProc]
			if rowStart <= row < rowEnd:
				#this is our start index
				break

		#reorganize distribIndex with pairIndex as starting point
		distribIndexList[curProc] = distribIndexList[curProc][listIndex:] + distribIndexList[curProc][:listIndex]

	return distribIndexList
		

def SetupStepList(globalSize, indexPairs, distribIndexList, distrib, rank):
	"""
	Sets up a list of MultiplyStepInfo for the current processor
	"""
	procCount = distrib.ProcCount
	procId = distrib.ProcId

	#Stepcount is the largest number of matrix elements any processor has
	stepCount = max([len(distribIndexList[i]) for i in range(procCount)])

	rowStart = array([distrib.GetLocalStartIndex(int(globalSize), rank, curProc) for curProc in range(procCount)], dtype=int)
	rowEnd = array(list(rowStart[1:]) + [globalSize])

	stepList = []
	for step in range(stepCount):
		matrixIndex = -1
		GlobalRow = -1
		GlobalCol = -1
		sendProc = -1
		recvProcList = []
		recvLocalRowList = []

		#find out which processors are sending to this
		#processor at this step
		for curProc in range(procCount):
			if len(distribIndexList[curProc]) > step:
				row, col = indexPairs[distribIndexList[curProc][step]]
				if rowStart[procId] <= row < rowEnd[procId] and procId != curProc:
					recvProcList.append(curProc)
					localRow = row - rowStart[procId]
					recvLocalRowList.append(localRow)

		#find out wether we have a matrix element this step
		if len(distribIndexList[procId]) > step:
			matrixIndex = step

			indexPairIndex = distribIndexList[procId][step]
			GlobalRow = indexPairs[indexPairIndex, 0]
			GlobalCol = indexPairs[indexPairIndex, 1]

			#and where to send it
			for curProc, (curStart, curEnd) in enumerate(zip(rowStart, rowEnd)):
				if curStart <= GlobalRow < curEnd:
					sendProc = curProc
					break
			if sendProc == -1:
				raise Exception("what!?")

		#add step
		stepList.append(MultiplyStepInfo(matrixIndex, GlobalRow, GlobalCol, sendProc, recvProcList, recvLocalRowList))
	
	return stepList

def StepListToArray(stepList):
	"""
	Convert a stepList to a set of arrays suitable for passing to fortran
	"""
	#find max number of recieves in a single step
	stepCount = len(stepList)
	maxRecvCount = max([len(step.RecvProcList) for step in stepList])

	localMatrixIndex = asarray([step.LocalMatrixIndex for step in stepList], dtype=int32)
	globalRow = asarray([step.GlobalRow for step in stepList], dtype=int32)
	globalCol = asarray([step.GlobalCol for step in stepList], dtype=int32)
	globalSendProc = asarray([step.SendProc for step in stepList], dtype=int32)
	recvCount = asarray([len(step.RecvProcList) for step in stepList], dtype=int32)
	
	recvProcList = asarray(-1 * ones((stepCount, maxRecvCount), dtype=int32), dtype=int32)
	recvLocalRowList = asarray(-1 * ones((stepCount, maxRecvCount), dtype=int32), dtype=int32)
	for i, step in enumerate(stepList):
		curRecvCount = len(step.RecvProcList)
		recvProcList[i,:curRecvCount] = step.RecvProcList
		recvLocalRowList[i,:curRecvCount] = step.RecvLocalRowList

	return localMatrixIndex, globalRow, globalCol, globalSendProc, recvProcList, recvLocalRowList, recvCount

