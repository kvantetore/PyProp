

def FindEigenvaluesNearShift(shift, **args):
	prop = SetupProblem(silent=True, eigenvalueShift=shift, **args)

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues()
	estimates = solver.Solver.GetConvergenceEstimates()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift

	sortIdx = argsort(real(eigenvalues))
	eigenvalues = eigenvalues[sortIdx]

	return eigenvalues


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
