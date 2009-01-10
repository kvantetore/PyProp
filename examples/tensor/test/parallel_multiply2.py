#execfile("example.py")
import boost.mpi as mpi

"""

IMPORTANT: This file is not being maintained. It was implemented
to serve as a test for the parallelization scheme for matrix-vector multiplication
tensorpotentialmultiply_generator.py can now generate routines for parallelization
even for multidimensional arrays. It might be easier to see how the algorithm
works in this file rather than in the generator file...

Global Matrix is N-by-N sparse matrix 

"""

from pylab import *
from numpy import *

ProcId = mpi.world.rank
ProcCount = mpi.world.size

"""
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

ArrayTranspose = pyprop.core.ArrayTranspose_1(1)
if ArrayTranspose != None:
	def GetGlobalStartIndex(size, procCount, procId):
		return ArrayTranspose.GetLocalStartIndex(size, 0, procId)
	
	def GetDistributedShape(size, procCount, procId):
		return ArrayTranspose.CreateDistributedShape(size, 0, procId)

else:
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

	colStart = array([distrib.GetLocalStartIndex(globalSize, rank, curProc) for curProc in range(procCount)], dtype=int)
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

	rowStart = array([distrib.GetLocalStartIndex(globalSize, rank, curProc) for curProc in range(procCount)], dtype=int)
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

def SetupSimpleDistributedMatrix(matrix, globalSize, indexPairs, distrib, rank):
	procCount = distrib.ProcCount
	procId = distrib.ProcId

	distribIndexList = SetupDistributedIndexList(globalSize, indexPairs, distrib, rank)
	stepList = SetupStepList(globalSize, indexPairs, distribIndexList, distrib, rank)

	localMatrix = zeros(len(distribIndexList[procId]), dtype=complex)

	for step in stepList:
		if step.LocalMatrixIndex != -1:
			localMatrix[step.LocalMatrixIndex] = matrix(step.GlobalRow, step.GlobalCol)

	return localMatrix, stepList

class DistributedModelTest:
	"""
	test implementation of DistributedModel 
	"""
	def __init__(self):
		self.ProcId = ProcId      	
		self.ProcCount = ProcCount	

	def GetLocalStartIndex(self, globalSize, rank, procId):
		return GetGlobalStartIndex(globalSize, self.ProcCount, procId)




def SimpleDistributedMatrixVectorMultiply(localMatrix, globalSize, stepList, localSource, localDest, comm):

	localSize = localDest.shape[0]
	assert(localDest.shape[0] == localSource.shape[0])
	globalStartIndex = GetGlobalStartIndex(globalSize, comm.size, comm.rank)

	"""
	For all matrix elements stored on this processor, we have the source data locally.
	We perform the multiplication for the matrix element, and send the data to the correct proc.
	"""
	for step in stepList:
		temp = 0
		localCol = step.GlobalCol - globalStartIndex
		localRow = step.GlobalRow - globalStartIndex

		#Post recvs for the data to be recieved this step
		requestList = [comm.irecv(recvProc, 0) for recvProc in step.RecvProcList]

		#If we have work to do this step
		sendReq = None
		if step.LocalMatrixIndex != -1:
			#Calculate dot product
			temp = localMatrix[step.LocalMatrixIndex] * localSource[localCol]

			#decide which procs to send to
			destProc = step.SendProc
			deltaProc = destProc - comm.rank
				
			#Store data if we have the row locally, or send it to the correct proc
			if deltaProc == 0:
				localDest[localRow] += temp
			else:
				sendReq = comm.isend(destProc, 0, temp)
		
		#Wait for data from other procs
		for localRow, req in zip(step.RecvLocalRowList, requestList):
			value, status = req.wait()
			localDest[localRow] += value

		if sendReq:
			sendReq.wait()

			
def RandomA(N, zeroRatio=0.9):
	A = rand(N*N).reshape(N,N)
	for row in range(N):
		for col in range(N):
			if A[row,col] < zeroRatio and row!=col:
				A[row, col] = 0
	A = mpi.broadcast(mpi.world, A, 0)
	if mpi.world.rank == 0:
		print "A = ", A

	return A


def test():

	A =[[ 1,   2,   3,  0,  4,  0 ], \
  	    [ 2,   5,   6,  7,  8,  9 ], \
	    [ 3,   6,  10, 11, 12,  0 ], \
	    [ 0,   7,  11, 13, 14, 15 ], \
	    [ 4,   8,  12, 14, 16, 17 ], \
	    [ 0,   9,   0, 16, 17, 18 ]]

	#A = RandomA(100, zeroRatio=0.95)
	#A = A + A.transpose()

	A = array(A)
	N = A.shape[0]

	indexPairs = []
	for row in range(N):
		for col in range(N):
			if A[row, col] != 0:
				indexPairs.append([row, col])
	indexPairs = array(indexPairs)

	def matrix(row, col):
		return A[row, col]

	#Create in-vector
	psi = rand(N) + 0.0j
	#psi = r_[:N] + 0.0j
	#send psi from proc0 to everyone
	psi = mpi.broadcast(mpi.world, psi, 0)

	#output
	refOutput = dot(A, psi)

	#Create local vectors and matrices
	localSize = GetDistributedShape(N, ProcCount, ProcId)
	globalStartIndex = GetGlobalStartIndex(N, ProcCount, ProcId)
	globalEndIndex = globalStartIndex+localSize

	localPsi = psi[globalStartIndex:globalEndIndex]
	localRefOutput = refOutput[globalStartIndex:globalEndIndex]

	distrib = DistributedModelTest()
	rank = 0
	localMatrix, stepList = SetupSimpleDistributedMatrix(matrix, N, indexPairs, distrib, rank)
	localMatrixIndex, globalRow, globalCol, globalSendProc, recvProcList, recvLocalRowList, recvCount = StepListToArray(stepList)

	print "SHAPE= ", recvLocalRowList.shape

	if ProcId == 10:
		print recvProcList
	
	#if ProcId == 0:
	#	print indexPairs

	
	localTestOutput = zeros(localSize, dtype=complex)
	
	for i in range(ProcCount):
		mpi.world.barrier()
		if i == ProcId:
			#if i == 0:
			#	print "psi =", psi
			print "ProcId == %i" % (i)
			print localSize
			print globalStartIndex, " -> ", globalEndIndex
			print "%i steps, %i nonzero steps" % (len(stepList), len([1 for step in stepList if step.LocalMatrixIndex!=-1]))
			print "max recv length = %i" % (max([len(step.RecvProcList) for step in stepList]))
			#for step in stepList:
			#	print step
			print ""
		mpi.world.barrier()

	#SimpleDistributedMatrixVectorMultiply(localMatrix, N, stepList, localPsi, localTestOutput, mpi.world)
	pyprop.core.TensorPotentialMultiply_SimpD(localMatrix, 1.0, localPsi, localTestOutput, N, localMatrixIndex, globalRow, globalCol, globalSendProc, recvProcList, recvLocalRowList, recvCount)
	
	
	#the verdict
	for i in range(ProcCount):
		if i == ProcId:
			#if i == 0:
			#	print ""
			#	print refOutput
			#	print ""
			print "ProcId == %i" % (i)
			print sqrt(sum(abs(localRefOutput)**2))
			print sqrt(sum(abs(localTestOutput)**2))
			print sqrt(sum(abs(localRefOutput - localTestOutput)**2))
			#print localTestOutput
			#print localRefOutput
			print ""
		mpi.world.barrier()

