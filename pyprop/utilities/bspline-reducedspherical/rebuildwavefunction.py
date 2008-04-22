from numpy import array
from numpy import meshgrid
try:
	import scipy.special
except:
	print "Could not load scipy, rebuilder will not work!"

class RankIndex:
	min = 0
	max = 1
	size = 2

class WavefunctionRebuilderBRS:
	"""
	"""

	def __init__(self, prop):
		self.Propagator = prop

	def ApplyConfigSection(self, configSection):
		
		self.RankInfo = []
		self.RankInfo += [configSection.rank0]
		self.RankInfo += [configSection.rank1]
		self.LMax = configSection.lmax

		self.GridExtent = (self.RankInfo[0][RankIndex.min], self.RankInfo[0][RankIndex.max], \
			self.RankInfo[1][RankIndex.min], self.RankInfo[1][RankIndex.max])

		self.InfoLevel = 0

	def Setup(self):
		self.x, self.z = self.SetupCartesianMesh()
		self.X, self.Z = meshgrid(self.x, self.z)
		self.R, self.Theta = self.SetupSphericalCoordinates(self.x, self.z)

		#Get BasePropagator
		if hasattr(self.Propagator.Propagator, "BasePropagator"):
			self.baseprop = self.Propagator.Propagator.BasePropagator
		else:
			self.baseprop = self.Propagator.Propagator

		#Get b-spline object
		self.bsplineObj = self.baseprop.SubPropagators[0].BSplineObject

		self.r_flat = self.R.ravel()

	def BuildWavefunction(self):
		"""
		"""

		try:
			pyprop.Redirect.Enable(True)
		
			xSize = self.x.size
			zSize = self.z.size
			
			psi = zeros((xSize, zSize))

			if self.InfoLevel > 0:
				pyprop.Redirect.redirect_stdout.stdout.write("l = ")

			for l in range(self.LMax + 1):
				if self.InfoLevel > 0:
					pyprop.Redirect.redirect_stdout.stdout.write("%i, " % l)
					pyprop.Redirect.redirect_stdout.stdout.flush()

				legendre = scipy.special.legendre(l)
				legendreGrid = legendre(cos(self.Theta))
				psiSlice = self.Propagator.psi.GetData()[:,l]
				tmpPsi = self.bsplineObj.ConstructFunctionFromBSplineExpansion(psiSlice, self.r_flat)
				tmpPsi = tmpPsi.reshape((xSize, zSize))
				psi += tmpPsi * legendreGrid

		finally:
			pyprop.Redirect.Disable()

		return self.X, self.Z, psi

	def SetupSphericalCoordinates(self, x, z):
		r = zeros((x.size, z.size))
		theta = zeros((x.size, z.size))

		for i, xi in enumerate(x):
			for j, zj in enumerate(z):
				r[i,j] = sqrt(xi**2 + zj**2)
				theta[i,j] = arctan2(abs(xi), zj)
		
		return r, theta
		

	def SetupCartesianMesh(self, **args):
		"""
		"""

		xMin = self.RankInfo[0][RankIndex.min]
		xMax = self.RankInfo[0][RankIndex.max]
		xSize = self.RankInfo[0][RankIndex.size]
		zMin = self.RankInfo[0][RankIndex.min]
		zMax = self.RankInfo[0][RankIndex.max]
		zSize = self.RankInfo[0][RankIndex.size]

		x = linspace(xMin, xMax, xSize)
		z = linspace(zMin, zMax, zSize)

		return x, z
