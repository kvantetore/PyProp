from numpy import *
from scipy.optimize import fsolve

DEBUG = False
		

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

def CenterExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gamma):
	"""
	Similar to 'ExponentialLinearBreakpointSequence', but densest region of points
	is in the middle. More precicely, the right half of the grid will be exponential-
	linear, and this will then be mirrored across the leftmost point to create the
	full grid. Note: Number of grid points must be even.
	"""

	#
	# Set up right part of grid
	# 
	centerPoint = (rmax + rmin) / 2.0

	#Setup inner/exponential region
	h = rpartition - centerPoint
	xi = [ centerPoint + h * ( exp(gamma * float(i - 1) / float(n -1)) - 1 ) / ( exp(gamma) - 1 ) \
	       for i in range(1, n + 1) ]

	#Setup outer/linear region
	h = xi[-1] - xi[-2]
	xi += [x for x in frange(rpartition + h, rmax, h)]

	#
	# Mirror to get left part of grid
	#
	fullGrid = [-x for x in reversed(xi[1:])] + xi

	return fullGrid


def VarExponentialLinearBreakpointSequence(rmin, rmax, minSep, maxSep, N):
	"""
	Same as VarExponentialLinearBreakpointSequence2, but solves a different set of
	equations.
	"""
	
	def pointsInner(h):
		return int(round(N - (rmax - rmin) / maxSep - h / maxSep - 1))

	def Equation1(gam, h, n):
		#n = pointsInner(h)
		return h * (exp(gam/(n-1.0)) - 1.0) / (exp(gam) - 1.0) - minSep

	def Equation2(gam, h, n):
		#n = pointsInner(h)
		return h * (exp(gam) - exp(gam * (n - 2.0)/(n - 1.0))) / (exp(gam) - 1.0) - maxSep

	def Equation3(gam, h, n):
		return N - (rmax - rmin) / maxSep - h / maxSep - 1 - n

	def func(x):
		out = [Equation1(x[0], x[1], x[2])]
		out.append(Equation2(x[0], x[1], x[2]))
		out.append(Equation3(x[0], x[1], x[2]))

		return out

	#Solve system of nonlinear equations to obtain gamma and h
	gam, h, n = fsolve(func, [2.0, 10.0, 20])

	n = ceil(n)

	#Calculate partition point
	rpartition = rmin + h

	if DEBUG:
		print "gamma = %s" % gam
		print "rpartition = %s" % rpartition
		print "number of inner points = %s" % n
		print "Equation 1 = %s" % Equation1(gam, h, n)
		print "Equation 2 = %s" % Equation2(gam, h, n)

	return ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gam)


def VarExponentialLinearBreakpointSequence2(rmin, rmax, minSep, maxSep, N):
	"""
	Compute exponential-linear sequence from a different set of parameters: grid extension,
	smallest and largest point seperatation, and total number of points. A non-linear system
	of equations is solved to find the parameters that best fulfill these conditions (gamma,
	number of inner points and extension of the exponential part). 
	"""
	
	def pointsInner(h):
		return int(round(N - (rmax - rmin) / maxSep - h / maxSep - 1))

	def Equation1(gam, rpartition, n):
		n = round(n)
		return abs(len(ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gam)) - N)

	def Equation2(gam, rpartition, n):
		n = round(n)
		xi = ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gam)
		return abs(xi[1] - xi[0] - minSep)

	def Equation3(gam, rpartition, n):
		n = round(n)
		xi = ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gam)
		return abs(xi[-2] - xi[-1] - maxSep)

	def func(x):
		out = [Equation1(x[0], x[1], x[2])]
		out.append(Equation2(x[0], x[1], x[2]))
		out.append(Equation3(x[0], x[1], x[2]))

		return out

	#Solve system of nonlinear equations to obtain gamma and h
	gam, rpartition, n = fsolve(func, [3.0, 10.0, 20])

	n = ceil(n)

	if DEBUG:
		print "gamma = %s" % gam
		print "rpartition = %s" % rpartition
		print "number of inner points = %s" % n
		print "Equation 1 = %s" % Equation1(gam, rpartition, n)
		print "Equation 2 = %s" % Equation2(gam, rpartition, n)
		print "Equation 3 = %s" % Equation2(gam, rpartition, n)

	return ExponentialLinearBreakpointSequence(rmin, rpartition, rmax, n, gam)
