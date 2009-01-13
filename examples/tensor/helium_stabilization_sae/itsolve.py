import sys

class IterCount:
	def __init__(self):
		self.TotalCount = 0

	def Callback(self, x):
		#sys.stdout.write("iter\n")
		self.TotalCount += 1


def JacobiPreconditioner(A, b, **args):	
	print "Applying Jacobi preconditioner..."
	invM = diag(1.0 / diag(A))
	preA = dot(invM, A)
	preb = dot(invM, b)
	return preA, preb


def BlockPreconditioner(A, b, **args):
	print "Applying Block preconditioner..."
	probSize = A.shape[0]
	blockSize = args["blockSize"]
	invM = zeros(shape(A))
	for k in range(probSize / blockSize):
		curSlice = [slice(k*(blockSize), (k+1)*blockSize)] * 2
		matrixBlock = A[curSlice]
		invM[curSlice] = inv(matrixBlock)
	preA = dot(invM, A)
	preb = dot(invM, b)
	return preA, preb


def SymmetricBlockPreconditioner(A, b, **args):
	print "Applying Symmetric Block preconditioner..."
	probSize = A.shape[0]
	blockSize = args["blockSize"]
	if mod(probSize, blockSize) != 0:
		raise Exception("Preconditioner block size must divide matrix size!")
	C = zeros(shape(A))
	invC = zeros(shape(A))
	for k in range(probSize / blockSize):
		curSlice = [slice(k*(blockSize), (k+1)*blockSize)] * 2
		matrixBlock = A[curSlice]
		L = linalg.cholesky(matrixBlock)
		C[curSlice] = L[:]
		invC[curSlice] = inv(L)

	preA = dot(invC, A)
	preA = dot(preA, transpose(invC))
	preb = dot(invC, b)
	return preA, preb, invC, C



def TestPreconditionedCG(blockSize = 11):
	#Setup overlap matrix
	prop = SetupProblem()
	S = SetupOverlapMatrix(prop).real

	#Random vector
	b = random.random(S.shape[0])

	#Iteration counter callback
	iterObj = IterCount()
	iterObjPC = IterCount()

	#Preconditioner
	preS, preb, invC, C = SymmetricBlockPreconditioner(S, b, blockSize=blockSize)

	#CG without preconditioning
	t1_nopc = time.time()
	y,info = scipy.linalg.cg(S, b, tol=1e-14, callback=iterObj.Callback)
	t2_nopc = time.time()

	#CG with preconditioning
	t1_pc = time.time()
	y,info = scipy.linalg.cg(preS, preb, tol=1e-14, callback=iterObjPC.Callback)
	t2_pc = time.time()

	#Exact solution
	x = linalg.solve(S, b)

	print "CG without preconditioner"
	print "   runtime = %3.2f s, iterations = %s" % (t2_nopc - t1_nopc, iterObj.TotalCount)
	print
	print "CG with preconditioner"
	print "   runtime = %3.2f s, iterations = %s" % (t2_pc - t1_pc, iterObjPC.TotalCount)
	print
	print "||x_exact - x_cg|| = %s" % linalg.norm(x - dot(transpose(invC), y))
