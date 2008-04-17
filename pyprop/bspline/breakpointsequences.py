from numpy import *
		

def ExponentialBreakpointSequence(rmin, rmax, n, gamma):
	"""
	Compute exponential breakpoint sequence on the interval [rmin, rmax] with n points.
	The parameter gamma determines point spacing:

			gamma -> 0    =>  linear sequence
			gamma -> \inf =>  all points exponentially close to rmin
	
	
	Breakpoints are computed according to

										  exp(gamma*(i-1)/(n-1)) - 1
			xi_i = rmin + (rmax - rmin) * ---------------------------
												exp(gamma) - 1
	 
	"""

	h = rmax - rmin
	xi = [ rmin + h * ( exp(gamma * float(i - 1) / float(n -1)) - 1 ) / ( exp(gamma) - 1 ) \
	       for i in range(1, n + 1) ] 

	return xi


def LinearBreakpointSequence(rmin, rmax, n):
	"""
	Compute linear breakpoint sequence on the interval [rmin, rmax] with n points.
	"""

	xi = [rmin + (rmax - rmin) * i / float(n - 1) for i in range(n)]

	return xi


def QuadraticLinearBreakpointSequence(rmin, rmax, n, joinPoint):
	"""
	Compute quadratic/linear breakpoint sequence on the interval [rmin, rmax] with n points.
	The sequence consists of quadratic and linear sub-sequences, joined at	'joinPoint'. 
	The quadratic sequence extends from 'rmin' to 'joinPoint'-1, while the linear
	sequence extends from 'joinPoint' to 'rmax'. This sequence is particularly suited for 
	problems where both bound and continuum states needs to be resolved.
	"""
	
	#joinPoint = joinPoint + 1
	jointPoint = float(joinPoint)
	rmin = float(rmin)
	rmax = float(rmax)

	# Compute first point of sequence
	r0 = (rmax * (joinPoint - 1) + rmin * (n - joinPoint)) / (2.0 * n - joinPoint - 1.0)

	# Scaling parameters for the two parts
	alpha = (r0 - rmin) / float(joinPoint - 1)**2
	beta = (rmax - r0) / float(n - joinPoint)
	
	# Parabolic part of sequence
	xi = [rmin + alpha * (i - 1)**2 for i in range(1, joinPoint)]
	
	# Linear part of sequence
	xi += [r0 + beta * (i - joinPoint) for i in range(joinPoint, n+1)]

	#print "r0 = ", r0, " alpha = ", alpha, " beta = ", beta

	return xi

def ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gamma):
	"""
	Compute exponential/linear breakpoint sequence. An inner region defined by
	[rmin, rpartition] have exponentially spaced points, whereas the outer region
	defined by (rpartition, rmax] is linearly spaced. The parameter 'n' specifies
	the number of points in the inner region. The outer region is constructed such
	that the grid point spacing equals that of the two last outermost points in the
	inner region.
	"""

	#Setup inner/exponential region
	h = rpartition - rmin
	xi = [ rmin + h * ( exp(gamma * float(i - 1) / float(n -1)) - 1 ) / ( exp(gamma) - 1 ) \
	       for i in range(1, n + 1) ]

	#Setup outer/linear region
	h = xi[-1] - xi[-2]
	xi += [x for x in frange(rpartition + h, rmax, h)]

	return xi

