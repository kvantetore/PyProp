from numpy import unique
#Import CoupledIndex from core, and add a pretty print to it
def CoupledIndexIter(self):
	yield self.l1
	yield self.l2
	yield self.L
	yield self.M

from libcoupledspherical import CoupledIndex
CoupledIndex.__str__ = lambda self: "l1=%i, l2=%i, L=%i, M=%i" % (self.l1, self.l2, self.L, self.M)
CoupledIndex.__repr__ = lambda self: "sys.modules['pyprop'].CoupledIndex(%i, %i, %i, %i)" % (self.l1, self.l2, self.L, self.M)	
CoupledIndex.__iter__ = CoupledIndexIter

class DefaultCoupledIndexIterator(object):
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
	def __init__(self, lmax, L=None, M=[0], fullTensor=False):

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
			
			
class EssentialStateCoupledIndexIterator(DefaultCoupledIndexIterator):
	"""
	Creates an iterator for yielding values of l1, l2, L, M that fullfill
	l1+l1 < lcutoff, to the CoupledSphericalHarmonicRepresentation

	lcutoff gives the maximum value of l1+l2
	
	Otherwise identical to DefaultCoupledIndexIterator, see that class for more
	information.
	 
	"""
	def __init__(self, lmax, lcutoff, L=None, M=[0], fullTensor=False):
		self.lcutoff = lcutoff
		super(EssentialStateCoupledIndexIterator, self).__init__(lmax, L, M, fullTensor)


	def __iter__(self):
		for cpldIdx in super(EssentialStateCoupledIndexIterator, self).__iter__():
			if (cpldIdx.l1 + cpldIdx.l2 < self.lcutoff):
				yield cpldIdx
	

	def __repr__(self):
		return "sys.modules['pyprop'].EssentialStateCoupledIndexIterator(%s, %s, L=%s, M=%s, fullTensor=%s)" % \
			(self.lmax, self.lcutoff, self.L, self.M, self.FullTensor)
