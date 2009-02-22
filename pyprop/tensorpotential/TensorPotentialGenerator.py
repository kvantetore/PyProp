#------------------------------------------------------------------------------------
#                       TensorPotentialGenerator main
#------------------------------------------------------------------------------------

def CreateBasisFromRepresentation(representation):
	if representation.__class__ == core.BSplineRepresentation:
		basis = BasisfunctionBSpline()
		basis.SetupBasis(representation.GetBSplineObject())
	
	elif representation.__class__ == core.ReducedSphericalHarmonicRepresentation:
		basis = BasisfunctionReducedSphericalHarmonic()
		basis.SetupBasis(representation.Range.MaxL)
	
	elif representation.__class__ == core.CoupledSphericalHarmonicRepresentation:
		basis = BasisfunctionCoupledSphericalHarmonic()
		basis.SetupBasis(representation)

	elif representation.__class__ == core.CustomGridRepresentation:
		basis = BasisfunctionFiniteDifference()
		basis.SetupBasis(representation)
	
	elif representation.__class__ == core.VectorRepresentation:
		basis = BasisfunctionVector()
		basis.SetupBasis(representation)

	else:
		raise NotImplementedException("Unknown representation %s" % representation)

	return basis


class TensorPotentialGenerator(object):

	def __init__(self, **args):
		self.BasisList = []

		if "config" in args:
			config = args["config"]
			self.Rank = config.Representation.rank
			self.SetupBasisFromConfig(config)

		elif "representation" in args:
			representation = args["representation"]	
			self.Rank = len(representation.GetFullShape())
			self.SetupBasisFromRepresentation(representation)

		else:
			raise Exception("Specify either 'config' or 'representation'")

	def SetupBasisFromConfig(self, config):
		self.BasisList = []
		for i in range(self.Rank):
			#Get the config section for this rank
			rankSection = config.GetSection("Rank%i" % i)
			#Create basisfunction object
			basisObject = rankSection.basis()
			rankSection.Apply(basisObject)
			self.BasisList.append(basisObject)

	def SetupBasisFromRepresentation(self, representation):
		self.BasisList = []
		for i in range(self.Rank):
			subRepr = representation.GetRepresentation(i)
			basisObject = CreateBasisFromRepresentation(subRepr)
			self.BasisList.append(basisObject)

	def GeneratePotential(self, configSection):
		"""
		0) Check the validity of the geometries of the potential
		1) Setup CombinedRepresentation
		2) Create wavefunction
		3) Create potential evaluator
		4) Evaluate potential on grid
		5) Represent potential in basis
		"""

		#Check geometry for each basis
		geometryList = self.GetGeometryList(configSection)

		#1) Create Representation:
		repr = self.SetupRepresentation(geometryList)

		#2) Create a wavefunction in the potential representation. 
		#This will be used to hold the and representation in order 
		#to use the existing potential evaluator framework in pyprop
		psi = CreateWavefunctionInstance(repr)

		#3) Create potential evaluator
		potentialEvaluator = self.SetupPotentialEvaluator(configSection, psi)
		
		#4) Initialize potential with basis pairs for custom potentials
		isCustomPotential = hasattr(potentialEvaluator, "SetBasisPairs")
		potentialShape = list(psi.GetRepresentation().GetInitialShape())
		if isCustomPotential:
			for rank, geometryInfo in enumerate(geometryList):
				if not geometryInfo.UseGridRepresentation():
					#We currently can not be distributed in a basis rank
					if repr.GetDistributedModel().IsDistributedRank(rank) and not geometryInfo.HasParallelMultiply():
						raise Exception("Representation %i can not be distributed during potential evaluation")

					#Get basis pairs from geometry info, and pass it to the potential
					potentialShape[rank] = geometryInfo.GetBasisPairCount()
					try:
						basisPairs = geometryInfo.GetBasisPairs()
						potentialEvaluator.SetBasisPairs(rank, basisPairs)
					except:
						print "BASISPAIRSHAPE=", basisPairs.shape
						raise


		#4) Evaluate the potential on the grid
		potentialData = CreateInstanceRank("core.DataBuffer", self.Rank)
		potentialData.ResizeArray(array(potentialShape))
		potentialEvaluator.UpdatePotentialData(potentialData.GetArray(), psi, 0, 0)
		
		# Save a copy of the potential if we're debugging
		debugPotential = getattr(configSection, "debug_potential", False)
		if debugPotential:	
			self.OriginalPotential = potentialData.GetArray().copy()
			
		#for parallelization
		#We only support the first rank initially distributed
		assert(repr.GetDistributedModel().GetDistribution()[0] == 0)
		assert(len(repr.GetDistributedModel().GetDistribution()) == 1)
		fullShape = repr.GetFullShape()
		if repr.GetDistributedModel().IsSingleProc():
			distribution = []
			origDistribution = []
		else:
			origDistribution = list(repr.GetDistributedModel().GetDistribution())
			distribution = list(repr.GetDistributedModel().GetDistribution())
		transpose = repr.GetDistributedModel().GetTranspose() 

		for distribRank in distribution:
			geomInfo = geometryList[distribRank]
			if not geometryList[distribRank].HasParallelMultiply():
				raise Exception("Distributed Rank %i (%s), has no support for parallel multiplication" % (distribRank, geomInfo.GetStorageId()))

		#5) Represent the potential in the bases

		source = potentialData
		dest = None
		del potentialData
		for rank in reversed(range(self.Rank)):
			basis = self.BasisList[rank]
			geometryInfo = geometryList[rank]

			#If we should not use the grid representation for this rank, it means that the potential
			#is already expressed in the basis representation, and thus we dont need to do anything
			#for this rank
			if geometryInfo.UseGridRepresentation():
				differentiation = 0
				if hasattr(configSection, "differentiation%i" % rank):
					differentiation = configSection.Get("differentiation%i" % rank)

				#Check if this rank is distributed. If it is, we must redistribute
				if rank in distribution:
					PrintOut( "Rank %i is distributed" % rank )
					assert(rank != self.Rank-1)
					assert(not rank+1 in distribution)

					#Create new distribution
					distribIndex = distribution.index(rank)
					newDistribution = list(distribution)
					newDistribution[distribIndex] = rank+1

					#Create shape of transposed function and allocate dest buffer
					transposedShape = transpose.CreateDistributedShape(fullShape, array(newDistribution, dtype=int32))
					
					del dest
					#dest = zeros(transposedshape, dtype=complex)
					dest = CreateInstanceRank("core.DataBuffer", self.Rank)
					dest.ResizeArray(array(transposedShape))

					#Transpose
					transpose.Transpose(fullShape, source.GetArray(), array(distribution, int), dest.GetArray(), array(newDistribution, int))

					#Use the new buffer
					del source
					source = dest
					distribution = newDistribution

				#Calculate dest shape
				destShape = array(source.GetArray().shape)
				destShape[rank] = geometryInfo.GetBasisPairCount()

				#Update full shape
				fullShape[rank] = geometryInfo.GetBasisPairCount()
			
				#Allocate the destination array
				del dest
				#dest = zeros(destShape, dtype=complex)
				dest = CreateInstanceRank("core.DataBuffer", self.Rank)
				dest.ResizeArray(array(destShape))

				#Represent this rank in the basis
				basis.RepresentPotentialInBasis(source.GetArray(), dest.GetArray(), rank, geometryInfo, differentiation) 

				#Use the destination from this rank as the source to the next
				del source
				source = dest

		#done!

		#Transpose back to the original distribution
		if distribution != origDistribution:
			#Create shape of transposed function and allocate dest buffer
			transposedShape = transpose.CreateDistributedShape(fullShape, array(origDistribution, dtype=int32))
			del dest
			#dest = zeros(transposedShape, dtype=complex)
			dest = CreateInstanceRank("core.DataBuffer", self.Rank)
			dest.ResizeArray(array(transposedShape))

			#Transpose
			transpose.Transpose(fullShape, source, array(distribution, int), dest, array(origDistribution, int))

			#Use the new buffer
			del source
			source = dest
			distribution = origDistribution

		del dest
		return source


		
	def GetGeometryList(self, configSection):
		"""
		Creates a list of geometry infos for each basos on the
		the given potential (configSection)
		"""
		geometryList = []
		for i, basis in enumerate(self.BasisList):
			geometryName = configSection.Get("geometry%i" % i)
			geometryInfo = basis.GetGeometryInfo(geometryName)
			geometryList.append(geometryInfo)

		return geometryList


	def SetupRepresentation(self, geometryList):
		"""
		Creates a CombinedRepresentation with subrepresentations to match the
		current potential
		"""

		#Create Combined Representation
		repr = CreateInstanceRank("core.CombinedRepresentation", self.Rank)
	
		#Set distributed model
		distrib = CreateDistribution(None, rank=self.Rank)
		repr.SetDistributedModel(distrib)
	
		#Setup sub representations
		for i, basis in enumerate(self.BasisList):
			geometryInfo = geometryList[i]
	
			#Create representation
			if geometryInfo.UseGridRepresentation():
				subRepr = basis.GetGridRepresentation()
			else:
				subRepr = basis.GetBasisRepresentation()
		
			#Set distributed model
			subDistrib = distrib.CreateSubDistributedModel()
			subRepr.SetDistributedModel(subDistrib)
			subRepr.SetBaseRank(i)
			repr.SetRepresentation(i, subRepr)


		return repr

	def SetupPotentialEvaluator(self, configSection, psi):
		evaluatorPrefix = "core.DynamicPotentialEvaluator"
		classname = configSection.classname
		potentialEvaluator = CreatePotentialInstance(classname, self.Rank, evaluatorPrefix)
		configSection.Apply(potentialEvaluator)
		if hasattr(potentialEvaluator, "Setup"):
			potentialEvaluator.Setup(psi)
			
		return potentialEvaluator
		


