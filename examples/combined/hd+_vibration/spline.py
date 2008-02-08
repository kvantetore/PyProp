from numpy import *

"""

if you have a set of function values y[i], evaluated on a set of grid points x[i],
and whish to interpolate these points to a new set of grid points x1[i], the following 
will do the trick

interp = Interpolator(x, y)
y1 = [interp.Evaluate(v) for v in x1]


"""

class Interpolator:
	"""
	A very simple spline interpolator. It is not very fast or memory
	efficient, but it works.

	For a more efficient implementation, see scipy.interpolate
	"""

	def __init__(self, x, y):
		self.GridPoints = array(x, dtype=double)
		self.GridValues = array(y, dtype=double)
		self.Precompute()

	def Precompute(self):
		"""
		Calculates the spline coefficients
		"""

		x = self.GridPoints
		y = self.GridValues
	
		N = len(x) - 1

		matrix = zeros((N+1,N+1), dtype=double)
		h = zeros(N+1, dtype=double)
		rhs = zeros(N+1, dtype=double)

		h[0] = x[1] - x[0]
		for i in range(1,N):
			h[i]   = x[i+1] - x[i]

			rhs[i] = 6 * ( (y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
			matrix[i, i-1] = h[i-1]
			matrix[i, i]   = 2 * (h[i-1] + h[i])
			matrix[i, i+1] = h[i]
			
		matrix = matrix[1:N, 1:N]
		rhs = rhs[1:N]

		x = linalg.solve(matrix, rhs)
		self.SplineCoeffs = zeros(N+1, dtype=double)
		self.SplineCoeffs[1:N] = x


	def FindIndex(self, x):
		minIndex = 0
		maxIndex = len(self.GridPoints)-1
		while maxIndex - minIndex > 1:
			minValue = self.GridPoints[minIndex]
			maxValue = self.GridPoints[maxIndex]

			nextIndex = (maxIndex + minIndex) / 2
			nextValue = self.GridPoints[nextIndex]

			if nextValue < x:
				minIndex = nextIndex
			else:
				maxIndex = nextIndex

		return minIndex
			
	def Evaluate(self, x):
		i = self.FindIndex(x)
		
		h = self.GridPoints[i+1] - self.GridPoints[i]
		x0 = self.GridPoints[i]
		z1 = self.SplineCoeffs[i]
		z2 = self.SplineCoeffs[i+1]
		y1 = self.GridValues[i]
		y2 = self.GridValues[i+1]

		a = (z2 - z1) / (6 * h)
		b = z1 / 2.
		c = (y2 - y1) / h - (2*h*z1 + h * z2) / 6.
		d = y1

		return a*(x-x0)**3 + b*(x-x0)**2 + c*(x-x0) + d


