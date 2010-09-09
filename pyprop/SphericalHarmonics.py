from numpy import unique
#Import CoupledIndex from core, and add a pretty print to it
def LMIndexIter(self):
	yield self.l
	yield self.m

from core import LmRange
LmRange.__str__ = lambda self: "l=%i, m=%i" % (self.l, self.m)
LmRange.__repr__ = lambda self: "sys.modules['pyprop'].LmRange(%i, %i)" % (self.l, self.m)	
LmRange.__iter__ = LMIndexIter

class DefaultCoupledIndexIterator:
	"""
	Creates an iterator for giving all possible values of l and m to
	the SphericalHarmonicRepresentation

	lmax gives the maximum value of l


	M can be
		1) a single value (default = 0)
		2) an iterator yielding the values of m

	yields a CoupledIndex(l1, l2, L, M) for every permissible combination of 
	l1, l2, L, M such that
	 
	    |l1-l2| <= L <= l1+l2. (Triangle inequality),
	    M = m1+ m2,
	    -L <= M <= L (See Mathworld, Wigner 3j symbol section)
	     
	For combinations not satisfying	these criteria, y{l1,l2,L,M} is zero 
	and can therefore safely be removed.

	This generator is most likely used in a configuration file like

	[AngularRepresentation]
	type = core.CoupledSphericalHarmonicRepresentation
	index_iterator = DefaultCoupledIndexIterator(lmax=4, L=[0,1,2], M=0)

	"""
	def __init__(self, lmax, m=[0], fullTensor=False):

		if L==None:
			L = range(lmax+1)

		#Make sure L and M are iterable
		if not hasattr(L, "__iter__"):
			L = [L]
		if not hasattr(L, "__iter__"):
			M = [M]
			
		#Make sure L and M are unique
		M = unique(M).tolist()
		L = unique(L).tolist()

		self.L = L
		self.M = M
		self.lmax = lmax
		self.FullTensor = fullTensor

	def __iter__(self):
		#Iterate through all permutations
		for curM in self.M:
			for curL in self.L:
				if (abs(curM) > curL) and (not self.FullTensor):
					continue
				for l1 in range(self.lmax+1):
					for l2 in range(self.lmax+1):
						#Only yield if l1, l2, L satisfies clebsch gordan condition
						#or we have explicitly asked fo the full tensor product of states
						if (abs(l1 - l2) <= curL <= l1 + l2) or self.FullTensor:
							yield CoupledIndex(l1, l2, curL, curM)
		

	def __repr__(self):
		return "sys.modules['pyprop'].DefaultCoupledIndexIterator(%s, L=%s, M=%s, fullTensor=%s)" % \
			(self.lmax, self.L, self.M, self.FullTensor)
