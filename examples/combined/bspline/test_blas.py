import pyprop

def TestMatrixVectorHermitianBanded():
	matrixSize = 5
	bands = 1
	lda = bands + 1

	A = zeros((lda, matrixSize), dtype=complex)
	x = ones(matrixSize, dtype=complex)
	y = zeros(matrixSize, dtype=complex)

	A[0,:] = [i+1 for i in range(matrixSize)]

	pyprop.core.MatrixVectorMultiplyHermitianBanded(A, x, y)


def TestMatrixVectorBanded():

	rows = 4
	cols = 5
	bands = 1
	lda = 2 * bands + 1
	diagIndex = (lda - 1) / 2
	diagLength = max(rows, cols)

	x = ones(rows, dtype=complex)
	y = zeros(cols, dtype=complex)
	A = zeros((rows, cols), dtype=complex)
	Ablas = zeros((cols, lda), dtype=complex)

	A += diag([i+1 for i in range(0,diagLength)])[:rows,:cols]
	A += diag([.5*(i+1) for i in range(0,diagLength)],1)[:rows,:cols]
	#A += diag([i+1 for i in range(0,diagLength)],-1)
	#A += diag([i+1 for i in range(0,diagLength)],1)

	for j in range(cols):
		iMin = max(j - bands, 0);
		iMax = min(j + bands, rows);
		for i in range(iMin, iMax):
			I = bands - j + i;
			J = j;
			Ablas[J,I] = A[i,j]

	pyprop.core.MatrixVectorMultiplyBanded(Ablas, x, y, 1+0j, 0j ,0)

	return A, Ablas,  x, y

