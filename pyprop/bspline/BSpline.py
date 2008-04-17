from numpy import *
execfile(__path__[0] + "/bspline/gausslegendrequadrature.py")
execfile(__path__[0] + "/bspline/breakpointsequences.py")

def InitBSpline(conf):
	bspline = BSPLINE()
	bspline.ApplyConfigSection(conf)
	bspline.CreateBSplineTable()
	bspline.CreateBSplineDerivative2Table()
	bspline.ComputeOverlapMatrix()
	bspline.SetupOverlapMatrixExpert()
	return bspline


class BSPLINE(core.BSpline):

	def ApplyConfigSection(self, conf):
		"""
		"""
		self.MaxSplineOrder = conf.order

		self.CreateBreakpointSequence(conf)
		self.CreateContinuitySequence(conf)
		self.CreateKnotSequence()

		# Get Gauss-Legendre quadrature nodes&weights
		self.SetupQuadratureRule(conf)

		self.SetupGlobalGridAndWeights()

		# Calculate number of b-splines. Since all b-splines extend over
		# exactly k+1 knotpoints, we must ensure that the last b-spline has
		# enough knotpoints left.
		self.NumberOfBSplines = self.GetKnotSequence().size - self.MaxSplineOrder

		# Set projection algorithm type
		self.ProjectionAlgorithm = conf.Get("projection_algorithm")

		#Set LAPACK solver algorithm
		if hasattr(conf, "lapack_algorithm"):
			self.LapackAlgorithm = conf.Get("lapack_algorithm")
		else:
			#Use the fast one as default
			self.LapackAlgorithm = 1

		
	def CreateBreakpointSequence(self, conf):
		"""
		Create Breakpointsequence of given type
		"""

		if(conf.bpstype == 'linear'):

			bpsSeq = LinearBreakpointSequence(conf.xmin, \
			                                  conf.xmax, \
			                                  conf.xsize) \

		elif(conf.bpstype == 'exponential'):

			bpsSeq = ExponentialBreakpointSequence(conf.xmin, \
			                                       conf.xmax, \
			                                       conf.xsize, \
			                                       conf.gamma)

		elif(conf.bpstype == 'quadraticlinear'):

			bpsSeq = QuadraticLinearBreakpointSequence(conf.xmin, \
			                                           conf.xmax, \
			                                           conf.xsize, \
			                                           conf.joinpoint)
		elif(conf.bpstype == 'exponentiallinear'):

			bpsSeq = ExponentialLinearBreakpointSequence(conf.xmin, \
			                                           conf.xpartition, \
			                                           conf.xmax, \
			                                           conf.xsize, \
			                                           conf.gamma)
		
		else:
			raise NameError, "Sequence not recognized!"
		
		self.NumberOfBreakpoints = len(bpsSeq)
		self.ResizeBreakpointSequence( self.NumberOfBreakpoints )
		self.GetBreakpointSequence()[:] = array(bpsSeq)


	def CreateContinuitySequence(self, conf):
		"""
		Create a sequence [v_i] defining the continuity condition at each breakpoint.
		The continuity condition is C^(v_i-1), where v_i is associated with
		breakpoint i. C^k-continuity implies continuous k'th derivative. 
		C^-1 means that the function (spline) itself is discontinuous.
		"""

		n_ = self.NumberOfBreakpoints
		self.ResizeContinuitySequence(n_)
		continuitySeq = self.GetContinuitySequence()

		# Vanilla: maximum continuity at interior points,
		#          maximum multiplicity at end points.
		if(conf.continuity == 'vanilla'):
			continuitySeq[:] = (self.MaxSplineOrder - 1) * ones(n_, dtype='int')
			continuitySeq[0] = 0
			continuitySeq[-1] = 0

		# Zero: maximum continuity at interior points,
		#       maximum multiplicity - 1 at end points.
		#       This will remove the nonzero spline and
		#       endpoints, giving zero boundary conditions.
		if(conf.continuity == 'zero'):
			continuitySeq[:] = (self.MaxSplineOrder - 1) * ones(n_, dtype='int')
			continuitySeq[0] = 1
			continuitySeq[-1] = 1


	def CreateKnotSequence(self):
		"""
		Create a B-spline knot sequence given a breakpoint sequence, continuity
		conditions at each point (continuitySequence) and the order of the B-splines.
		"""

		knotSequence = []
		n_ = self.NumberOfBreakpoints
		
		#knotBreakpointMap = []
		topKnotMap = []
		
		knotIndex = -1
		for i in range(n_):
			knotPointMultiplicity = self.MaxSplineOrder - self.GetContinuitySequence()[i]
			knotIndex += knotPointMultiplicity
			for j in range(knotPointMultiplicity):
				knotSequence.append(self.GetBreakpointSequence()[i])
				topKnotMap.append(knotIndex)
				#topKnotMap.append(knotIndex)

		self.ResizeKnotSequence( len(knotSequence) )
		self.GetKnotSequence()[:] = array(knotSequence)

		#self.ResizeKnotBreakpointMap( len(knotBreakpointMap) )
		#self.GetKnotBreakpointMap()[:] = array(knotBreakpointMap)

		self.ResizeTopKnotMap( len(topKnotMap) )
		self.GetTopKnotMap()[:] = array(topKnotMap)


	def SetupQuadratureRule(self, conf):
		"""
		Obtain nodes and weights for Gauss-Legendre quadrature on each
		breakpoint interval.
		"""
	
		# Get points&weights
		quadOrder = conf.Get("quad_order_additional") + conf.Get("order")
		quadrule = GaussQuadratureRule(ruleOrder = quadOrder);

		# Set storage vectors to correct size
		self.ResizeWeights( len(quadrule.Weights) );
		self.ResizeNodes( len(quadrule.Nodes) );

		# Store points and weights
		self.GetNodes()[:] = quadrule.Nodes[:]
		self.GetWeights()[:] = quadrule.Weights[:]


	def SetupGlobalGridAndWeights(self):
		"""
		Using Gauss-Legendre nodes and weights for each breakpoint interval,
		we make a global grid covering [x_min, x_max] by concatenating these.
		A linear transformation is applied on the coordinates to map the nodes
		to a given interval [a, b] (a, b in [x_min, x_max]) and the weights
		are scaled accordingly. In addition, we compute a map from knot point
		index to grid point index (one-to-one).
		"""

		nodes = self.GetNodes()
		weights = self.GetWeights()

		# Resize global grid and weights vectors
		numberOfNodes = nodes.shape[0]
		globalGridSize = numberOfNodes * (self.NumberOfBreakpoints - 1)
		self.ResizeQuadratureGridGlobal(globalGridSize)
		self.ResizeQuadratureWeightsGlobal(globalGridSize)
		
		# This is the knot index -> grid index map
		#self.ResizeKnotGridIndexMap(globalGridSize)

		globalGrid = []
		globalWeights = []
		knotToGridIndexMap = []
		gridPointCounter = 0
		for j in range(self.NumberOfBreakpoints - 1):
			a = self.GetBreakpointSequence()[j]
			b = self.GetBreakpointSequence()[j+1]

			# Take knot point multiplicity into account.
			# Degenerate knots map to the same grid index.
			knotPointMultiplicity = self.MaxSplineOrder - self.GetContinuitySequence()[j]
			for k in range(knotPointMultiplicity):
				knotToGridIndexMap.append(gridPointCounter)

			for i in range(numberOfNodes):
				x = nodes[i]
				w = weights[i]
				y = self.ScaleAndTranslate(x, a, b)

				globalGrid += [y]
				globalWeights += [(b - a) / 2.0 * w]

				gridPointCounter += 1;
					
		# Store global grid and weights
		self.GetQuadratureGridGlobal()[:] = array(globalGrid)[:]
		self.GetQuadratureWeightsGlobal()[:] = array(globalWeights)[:]

		# Store knot->grid map
		self.ResizeKnotGridIndexMap(len(knotToGridIndexMap))
		self.GetKnotGridIndexMap()[:] = array(knotToGridIndexMap)[:]
