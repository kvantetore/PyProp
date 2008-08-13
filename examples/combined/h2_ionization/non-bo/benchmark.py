ExactValue0 = -1.845261353672151
ExactValue1 = -1.682283647669635

def BenchmarkImtime(**args):
	"""
	Measures the time it takes to run imtime to completion
	with the supplied args, as well as the calculated ground
	state energy
	"""
	
	print args

	args["imtime"] = True
	args["silent"] = True

	duration = -time.time()
	prop = SetupProblem(**args)
	for t in prop.Advance(False): pass
	duration += time.time()

	energy = prop.GetEnergy()

	return duration, energy


def BenchmarkEigenstate(**args):
	"""
	Measures the time it takes to run imtime to completion
	with the supplied args, as well as the calculated ground
	state energy
	"""
	
	print args

	args["imtime"] = True
	args["silent"] = True

	duration = -time.time()
	solver = FindEigenstates(**args)
	duration += time.time()

	energy = solver.GetEigenvalues()[1]

	return duration, energy



def RunBenchmarkTimestep(**args):
	dtList = 0.5**r_[0:12]
	bench = array([BenchmarkImtime(dt=dt, **args) for dt in dtList])

	exactValue = bench[-1, 1]
	error = abs((bench[:,1] - ExactValue1) / ExactValue1)
	time = bench[:,0]

	semilogy(dtList, error)

	return dtList, time, error


def RunBenchmarkBSplineOrder(**args):
	orderList = r_[3:10]
	bench = array([BenchmarkImtime(bsplineOrder=order, **args) for order in orderList])

	error = abs((bench[:,1] - ExactValue1) / ExactValue1)
	time = bench[:,0]

	semilogy(time, error)

	return orderList, time, error


def RunBenchmarkDifferenceOrder(**args):
	orderList = r_[3:19:2]
	bench = array([BenchmarkEigenstate(differenceOrder=order, **args) for order in orderList])

	exactValue = bench[-1,1]
	error = abs((bench[:,1] - exactValue) / exactValue)
	time = bench[:,0]

	print bench[:,1]

	semilogy(orderList, error)

	return orderList, time, error

def RunBenchmarkPolarCount(**args):
	countList = array([16,32,48,64,96,128])
	bench = array([BenchmarkImtime(polarCount=polarCount, **args) for polarCount in countList])

	error = abs((bench[:,1] - ExactValue1) / ExactValue1)
	time = bench[:,0]

	semilogy(countList, error)

	return coutList, time, error
