from pyprop.utilities.folderwatch import FolderWatch
import os

def FindEigenvaluesNearShift(shift, **args):
	messageId = args["messageId"]

	#Setup output files
	messagesFolder = "output/messages/"
	dataFolder = "output/data/"
	if not os.path.exists(messagesFolder):
		os.makedirs(messagesFolder)
	if not os.path.exists(dataFolder):
		os.makedirs(dataFolder)
	messageFile = os.path.join([messagesFolder, "message_%s.h5" % (messageId)])
	eigenvectorFile = os.path.join([dataFolder, "eigenvectors_%s.h5" % (messageId)])

	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, **args)

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#save eigenvectors
	SaveEigenvectors(eigenvectorFile, solver, shift)

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues()
	estimates = solver.Solver.GetConvergenceEstimates()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift

	sortIdx = argsort(real(eigenvalues))
	eigenvalues = eigenvalues[sortIdx]

	#Save message file
	f = table.openFile(messageFile, "w"):
	f.createArray(f.root, "eigenvalues", eigenvalues)
	f.close()


def RunSubmit(function, procCount=1, procPerNode=4, *arglist, **argdict):
	if isinstance(function, str):
		arg1 = function
	else
		arg1 = function.func_name

	arg2 = commands.mkarg(repr(arglist))
	arg3 = commands.mkarg(repr(argdict))

	if INSTALLATION == "hexagon":
		submit = submitpbs.SubmitScript()
		submit.proc_count = procCount
		submit.ppn = min(procPerNode, procCount)
		submit.executable = "./python-exec run-function.py"
		submit.parameters = arg1 + arg2 + arg3
		submit.WriteScript("test.job")

	elif INSTALLATION == "stallo":
		raise Exception("please to implement")
	
	elif INSTALLATION == "local":
		raise Exception("please to implement")
	
	else:
		raise Exception("Unknown installation '%s'" % INSTALLATION)

def FindEigenvalueSpectrum(startEigenvalue, endEigenvalue, **args):
	evStart = FindEigenvaluesNearShift(startEigenvalue, **args)	
	evEnd = FindEigenvaluesNearShift(endEigenvalue, **args)

	eigenvalues = {startEigenvalue: evStart, endEigenvalue:evEnd}
	activeIntervals = [(startEigenvalue, endEigenvalue)]

	while len(activeIntervals) > 0:
		newIntervals = []
		for i, (start, end) in enumerate(list(activeIntervals)):
			evStart = eigenvalues[start]
			evEnd = eigenvalues[end]
		
			#if evStart end evEnd is not overlapping, we should split the interval
			if real(evStart[-1]) < real(evEnd[0]):
				midpoint = (start + end) / 2.
			
				#Find eigenvalues at midpoint
				evMidpoint = FindEigenvaluesNearShift(midpoint, **args)
				eigenvalues[midpoint] = evMidpoint

				#mark subintervals as active
				newIntervals.append((start, midpoint))
				newIntervals.append((midpoint, end))

		#the new intervals are now the active ones
		activeIntervals = newIntervals

	eigenvalueList = []
	overlappingEigenvalues = 0
	for ev in sorted(eigenvalues.keys()):
		for curE in eigenvalues[ev]:
			if len(eigenvalueList) == 0 or curE > eigenvalueList[-1]:
				eigenvalueList.append(curE)
			else:
				overlappingEigenvalues += 1

	print "Eiganvalues found       = %i" % (len(eigenvalueList))
	print "Overlapping eigenvalues = %i" % (overlappingEigenvalues)

	return eigenvalueList
