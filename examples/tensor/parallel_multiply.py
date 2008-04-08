
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

k is the number of super diagonals (2 in the above ex)
A(i,j) = Packed(j, j-i+k)

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
def BandedMatrixVectorMultiply(localMatrix, size, bands, localSource, localDest, procCount, procId):

	localSize = localDest.shape[0]
	assert(localSize == size / procCount)
	assert(localDest.shape[0] == localSource.shape[0])
	assert(localMatrix.shape[0] == localSize)
	assert(localMatrix.shape[1] == 2*bands+1)
	globalStartIndex = procId * localSize

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

				globalPackedRow = globalCol
				globalPackedCol = globalCol - globalRow + bands

				#globalPackedCol = bands + 1 - globalCol + globalRow

				localPackedRow = globalPackedRow - globalStartIndex
				localPackedCol = globalPackedCol

				#print localPackedRow, localPackedCol
				#localPackedRow = localIndex
				#localPackedCol = localIndex - localRow - bands
				
				tempValue += localMatrix[localPackedRow, localPackedCol] * localSource[localIndex]

		#decide which procs to send or recieve from 
		deltaProc = localRow/localSize
		destProc = procId + deltaProc
		sourceProc = procId - deltaProc

		if deltaProc == 0:
			localDest[localRow] += tempValue
		else:
			recvTemp = zeros(1, dtype=localSource.dtype)
			sendTemp = array(tempValue, dtype=localSource.dtype)

			#NB! This is NOT the way to do it, one should post and irecv, and then 
			#use send. however, pypar does not have isend or irecv,  so we'll have
			#to settle with this
			#
			#this actual code will not scale as expected, as it will effectively run in serial
			if 0 <= destProc < procCount:
				#if we have data to send
				pypar.send(sendTemp, destProc, use_buffer=True)

			if 0 <= sourceProc < procCount:
				#if we have data to recv
				pypar.receive(sourceProc, recvTemp)

				#decide where to put this value
				sourceGlobalStartIndex = sourceProc * localSize
				sourceGlobalRow = localRow + sourceGlobalStartIndex
				sourceLocalRow = sourceGlobalRow - globalStartIndex

				#add at the correct row
				localDest[sourceLocalRow] += recvTemp[0]





def test():
	xmax = 2.
	N = 10
	bands = 4 
	
	dx = 2 * xmax / N

	#Create original matrix
	A = zeros((N, N), dtype=double)
	for i in range(N):
		x = -xmax + i*dx
		for j in range(-bands, bands+1):
			if 0 <= j+i < N:
				A[i,j+i] = x**2 / (abs(j)+1)
	

	#Create packed matrix
	PackedA = zeros((N, 2*bands+1), double)
	for i in range(N):
		for j in range(-bands, bands+1):
			if 0 <= j+i < N:
				row = i
				col = i+j
				PackedA[col, col-row+bands] = A[row, col]

	"""
	figure()
	imshow(A, interpolation="nearest")
	figure()
	imshow(PackedA, interpolation="nearest")
	"""

	#Create in-vector
	psi = rand(N)
	#send psi from proc0 to everyone
	pypar.broadcast(psi, 0)

	#output
	refOutput = dot(A, psi)

	#Create local vectors and matrices
	localSize = N / ProcCount
	localPackedA = PackedA[localSize*ProcId:localSize*(ProcId+1), :]
	localPsi = psi[localSize*ProcId:localSize*(ProcId+1)]
	localRefOutput = refOutput[localSize*ProcId:localSize*(ProcId+1)]

	localTestOutput = zeros(localSize, dtype=double)
	BandedMatrixVectorMultiply(localPackedA, N, bands, localPsi, localTestOutput, ProcCount, ProcId)

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
			print localTestOutput
			print localRefOutput
			print ""
		pypar.barrier()
