execfile("example.py")

"""

Global Matrix is N-by-N banded Hermitian matrix with
k offdiagonal bands

The matrix divided on ProcCount processors

Matrix is stored in BLAS format, where the current 
processors has stored all rows i of the matrix, and only 
cols j, where startCol_p <= j < startCol_p+colCount_p

Assume N % ProcCount = 0

colCount_p = N / ProcCount
startCol_p = procRank * (N / ProcCount)

"""

from pylab import *
from numpy import *

import pypar

ProcId = pypar.rank()
ProcCount = pypar.size()

"""
Variables without prefix are the same on every processor
Variables prefixed with local are indexed relative to the starting point of the current processor
Variables prefixed with global are indexed absolutely in the global arrays

A      = [ a00 a01 a02             ]
         [ a10 a11 a12 a13         ]
		 [ a20 a21 a22 a23 a24     ]
		 [     a31 a32 a33 a34 a35 ]
		 [         a42 a43 a44 a45 ]
		 [             a53 a54 a55 ]

Packed = [ a20 a10 a00   0   0 ] 
         [ a31 a21 a11 a01   0 ]
		 [ a42 a32 a22 a12 a02 ]
		 [ a53 a43 a33 a23 a13 ]
		 [   0 a54 a44 a34 a24 ]
		 [   0   0 a55 a45 a35 ]

Packed = [   0   0 a00 a10  a20] 
         [   0 a01 a11 a21  a31]
		 [ a02 a12 a22 a32  a42]
		 [ a13 a23 a33 a43  a53]
		 [ a24 a34 a44 a54    0]
		 [ a35 a45 a55   0    0]



k is the number of super diagonals (2 in the above ex)
A(i,j) = Packed(j, k+1 - j + i )

Example:
N = 6
P = 3

this corresponds to the following distribution
of the original matrix

            proc0     proc1     proc2
             0
             0   0
A      = [ a00 a01 | a02     |         ]
         [ a10 a11 | a12 a13 |         ]
		 [ a20 a21 | a22 a23 | a24     ]
		 [     a31 | a32 a33 | a34 a35 ]
		 [         | a42 a43 | a44 a45 ]
		 [         |     a53 | a54 a55 ]
		                         0   0
								     0

localMatrix0 = [ a20 a10 a00   0   0 ] 
               [ a31 a21 a11 a01   0 ]

localMatrix1 = [ a42 a32 a22 a12 a02 ]
		       [ a53 a43 a33 a23 a13 ]

localMatrix2 = [   0 a54 a44 a34 a24 ]
  		       [   0   0 a55 a45 a35 ]


	         
             0
             0   0
A      = [ a00 a01 | a02     |     | ]
         [ a10 a11 | a12 a13 |     | ]
		 [ a20 a21 | a22 a23 | a24 | ]
		 [     a31 | a32 a33 | a34 | ]
		 [         | a42 a43 | a44 | ]
		                   0     0   
								 0    
		


both input and output psi is distributed the same way
            proc0      proc1     proc2
psi =    [  c0  c1 |  c2  c3 |  c4  c5 ]

The algorithm is then that each processor iterates over 
all the rows it has, calculates the product for that row
and sends it to the correct processor.
In order to ensure proper load balancing, proc0 should not do
anything on the first two iterations other than recieve row0 and row1
from proc1

"""

def MapRowColToPacked(row, col, size, bands):
	packedRow = col
	packedCol = bands - col + row
	return packedRow, packedCol

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

def GetOwnerProcId(size, procCount, globalIndex):
	rest = size % procCount
	if rest != 0:
		size += procCount - rest
	
	return globalIndex / (size / procCount)
	
	


def BandedMatrixVectorMultiply(localMatrix, size, bands, localSource, localDest, procCount, procId):

	localSize = localDest.shape[0]
	#assert(localSize == size / procCount)
	assert(localDest.shape[0] == localSource.shape[0])
	assert(localMatrix.shape[0] == localSize)
	assert(localMatrix.shape[1] == 2*bands+1)
	globalStartIndex = GetGlobalStartIndex(size, procCount, procId)

	"""
	For all rows which has nonzero elements on this processor
	"""
	for localRow in range(-bands, localSize + bands + 1):
		tempValue = 0
		globalRow = localRow + globalStartIndex
		if 0 <= globalRow < size:
			#Calculate dot product
			localStartIndex = max(0, localRow - bands)
			localEndIndex = min(localRow + bands + 1, localSize)
			for localIndex in range(localStartIndex, localEndIndex):
				globalCol = localIndex + globalStartIndex	

				globalPackedRow, globalPackedCol = MapRowColToPacked(globalRow, globalCol, size, bands)
				
				localPackedRow = globalPackedRow - globalStartIndex
				localPackedCol = globalPackedCol

				tempValue += localMatrix[localPackedRow, localPackedCol] * localSource[localIndex]

		#decide which procs to send or recieve from 
		destProc = GetOwnerProcId(size, procCount, globalRow)
		deltaProc = destProc - procId
		sourceProc = procId - deltaProc
	
		#for i in range(
		#print procId, ", ", deltaProc, ", ", destProc

		if deltaProc == 0:
			if 0 <= globalRow < size:
				localDest[localRow] += tempValue

		else:
			recvTemp = zeros(1, dtype=localSource.dtype)
			sendTemp = array(tempValue, dtype=localSource.dtype)

			#NB! This is NOT the way to do it, one should post and irecv, and then 
			#use send. however, pypar does not have isend or irecv,  so we'll have
			#to settle with this
			#
			#this actual code will not scale as expected, as it will effectively run in serial
			if 0 <= destProc < procCount and 0 <= globalRow < size:
				#if we have data to send
				pypar.send(sendTemp, destProc, use_buffer=True)

			sourceGlobalStartIndex = GetGlobalStartIndex(size, procCount, sourceProc)
			sourceGlobalRow = localRow + sourceGlobalStartIndex
			if 0 <= sourceProc < procCount and 0 <= sourceGlobalRow < size:
				#if we have data to recv
				pypar.receive(sourceProc, recvTemp)

				#decide where to put this value
				sourceLocalRow = sourceGlobalRow - globalStartIndex

				#add at the correct row
				localDest[sourceLocalRow] += recvTemp[0]





def test():
	xmax = 2.
	N = 5
	bands = 3 
	
	dx = 2 * xmax / N

	#Create original matrix
	A = zeros((N, N), dtype=complex)
	for i in range(N):
		x = -xmax + i*dx
		for j in range(-bands, bands+1):
			if 0 <= j+i < N:
				A[i,j+i] = x**2 / (abs(j)+1)
	

	#Create packed matrix
	PackedA = zeros((N, 2*bands+1), complex)
	for i in range(N):
		for j in range(-bands, bands+1):
			if 0 <= j+i < N:
				row = i
				col = i+j
				packedRow, packedCol = MapRowColToPacked(row, col, N, bands)
				PackedA[packedRow, packedCol] = A[row, col]

	"""
	figure()
	imshow(A, interpolation="nearest")
	figure()
	imshow(PackedA, interpolation="nearest")
	"""

	#Create in-vector
	psi = rand(N) + 0.0j
	#send psi from proc0 to everyone
	pypar.broadcast(psi, 0)

	#output
	refOutput = dot(A, psi)

	#Create local vectors and matrices
	localSize = GetDistributedShape(N, ProcCount, ProcId)
	globalStartIndex = GetGlobalStartIndex(N, ProcCount, ProcId)
	globalEndIndex = globalStartIndex+localSize

	localPackedA = PackedA[globalStartIndex:globalEndIndex, :]
	localPsi = psi[globalStartIndex:globalEndIndex]
	localRefOutput = refOutput[globalStartIndex:globalEndIndex]

	localTestOutput = zeros(localSize, dtype=complex)
	
	for i in range(ProcCount):
		if i == ProcId:
			print "ProcId == %i" % (i)
			print localSize
			print globalStartIndex, " -> ", globalEndIndex
			print ""
		pypar.barrier()
	
	#BandedMatrixVectorMultiply(localPackedA, N, bands, localPsi, localTestOutput, ProcCount, ProcId)
	#BandedMatrixMultiply_Wrapper(localPackedA.reshape(localPackedA.size), 1.0, localPsi, localTestOutput, N, bands)
	TensorPotentialMultiply_BandedDistributed(localPackedA.reshape(localPackedA.size), 1.0, localPsi, localTestOutput, N, bands)

	#the verdict
	for i in range(ProcCount):
		if i == ProcId:
			if i == 0:
				print ""
				print refOutput
				print ""
			print "ProcId == %i" % (i)
			print sqrt(sum(abs(localRefOutput)**2))
			print sqrt(sum(abs(localTestOutput)**2))
			print sqrt(sum(abs(localRefOutput - localTestOutput)**2))
			#print localTestOutput
			#print localRefOutput
			print ""
		pypar.barrier()
