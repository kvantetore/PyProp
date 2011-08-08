from numpy import array, double, exp, log, r_, sign, sqrt

"""
Grid generation for CustomGridRepresentation

"""
from numpy import array, r_, log, exp, sign, double, sqrt
import copy


def GetGridLinear(conf, xmax=None, xmin=None):
	if xmin == None: xmin = conf.xmin
	if xmax == None: xmax = conf.xmax
	count = conf.count

	start = xmin
	end = xmax

	if not conf.include_left_boundary:
		count += 1
	if not conf.include_right_boundary:
		count += 1

	dx = (xmax - xmin) / float(count-1)
	if not conf.include_left_boundary:
		start += dx
	if conf.include_right_boundary:
		end += dx
	grid = r_[start:end:dx]

	return array(grid, dtype=double)

def GetGridQuadratic(conf):
	xmin = sign(conf.xmin) * sqrt(abs(conf.xmin))
	xmax = sign(conf.xmax) * sqrt(abs(conf.xmax))
	linearGrid = GetGridLinear(conf,xmin=xmin, xmax=xmax)
	grid = sign(linearGrid) * linearGrid**2
	return array(grid, dtype=double)

def GetBidirectionalGridExponential(conf):
	xmax = float(conf.xmax)
	#count = conf.count
	gamma = float(conf.gamma)

	xmin = - log(xmax+1) / gamma
	xmax = log(xmax+1) / gamma
	linearGrid = GetGridLinear(conf, xmin=xmin, xmax=xmax)

	grid = sign(linearGrid) * (exp(gamma * abs(linearGrid)) - 1)

	return grid

def GetBidirectionalGridExponentialLinear(conf):
	innerBoundary = float(conf.inner_boundary)
	outerBoundary = float(conf.outer_boundary)
	count = conf.inner_count
	gamma = conf.gamma

	i = array(r_[0:count], dtype=double)
	n = float(count)

	#Setup inner (exponential) region
	innerGrid = innerBoundary * (exp(gamma*i/n) - 1 ) / ( exp(gamma) - 1 ) 

	#Setup outer (linear) region
	h = innerGrid[-1] - innerGrid[-2]
	positiveGrid = array(list(innerGrid) + list(r_[h+innerGrid[-1]:outerBoundary:h]))

	#Create negative grid values
	grid = array(list(reversed(-1*positiveGrid))[:-1] + list(positiveGrid))
	
	return grid


def GetGridBilinear(conf):
	"""Return a "bilinear" grid, i.e. two joined grids with different linear spacings.

	Config options
	--------------
	xmin:         leftmost grid point
	xmax:         rightmost grid point
	inner_count:  number of grid points of inner (leftmost) grid portion
	outer_count:  number of grid points of outer (rightmost) grid portion
	partition:    divions point between inner and outer grid
	include_left_boundary
	include_right_boundary

	"""
	localConf = copy.copy(conf)

	#Setup inner grid
	localConf.Set("count", conf.inner_count)
	localConf.Set("xmin", conf.xmin)
	localConf.Set("xmax", conf.partition)
	localConf.Set("include_right_boundary", True)
	innerGrid = list(GetGridLinear(localConf))

	#Setup outer grid
	localConf.Set("count", conf.outer_count)
	localConf.Set("xmin", conf.partition)
	localConf.Set("xmax", conf.xmax)
	localConf.Set("include_left_boundary", False)
	localConf.Set("include_right_boundary", conf.include_right_boundary)
	outerGrid = list(GetGridLinear(localConf))

	return array(innerGrid + outerGrid)


