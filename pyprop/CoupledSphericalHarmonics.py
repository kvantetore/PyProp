
#Import CoupledIndex from core, and add a pretty print to it
from core import CoupledIndex
CoupledIndex.__str__ = lambda self: "l1=%i, l2=%i, L=%i, M=%i" % (self.l1, self.l2, self.L, self.M)

class DefaultCoupledIndexIterator:
	"""
	Creates an iterator for giving all possible values of l1, l2, L, M to
	the CoupledSphericalHarmonicRepresentation

	lmax gives the maximum l value of either l1 or l2

	L can be
		1) None: range(lmax+1) (default)
		2) a single value
		3) an iterator yielding the values of L

	M can be
		1) a single value (default = 0)
		2) an iterator yielding the values of M

	yields a CoupledIndex(l1, l2, L, M) for every permissible combination of 
	l1, l2, L, M such that |l1-l2| <= L <= l1+l2. For combinations not satisfying
	this criterion, y{l1,l2,L,M} is zero and can therefore safely be removed

	This generator is most likely used in a configuration file like

	[AngularRepresentation]
	type = core.CoupledSphericalHarmonicRepresentation
	index_iterator = DefaultCoupledIndexIterator(lmax=4, L=[0,1,2], M=0)

	"""
	def __init__(self, lmax, L=None, M=[0], fullTensor=False):

		if L==None:
			L = range(lmax+1)

		#Make sure L and M are iterable
		if not hasattr(L, "__iter__"):
			L = [L]
		if not hasattr(L, "__iter__"):
			M = [M]

		self.L = L
		self.M = M
		self.lmax = lmax
		self.FullTensor = fullTensor

	def __iter__(self):
		#Iterate through all permutations
		for curM in self.M:
			for curL in self.L:
				for l1 in range(self.lmax+1):
					for l2 in range(self.lmax+1):
						#Only yield if l1, l2, L satisfies clebsch gordan condition
						#or we have explicitly asked fo the full tensor product of states
						if (abs(l1 - l2) <= curL <= l1 + l2) or self.FullTensor:
							yield CoupledIndex(l1, l2, curL, curM)
		

	def __repr__(self):
		return "sys.modules['pyprop'].DefaultCoupledIndexIterator(%s, L=%s, M=%s, fullTensor=%s)" % \
			(self.lmax, self.L, self.M, self.FullTensor)
