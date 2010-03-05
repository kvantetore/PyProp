
class CombinedPropagator(PropagatorBase):
	__Base = PropagatorBase

	TransposeForward = 1
	TransposeBackward = -1
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
	
		#Create all sub propagators
		self.SubPropagators = []
		for i in range(self.Rank):
			sectionName = config.Propagation.Get("propagator" + str(i))
			if sectionName != None:
				section = config.GetSection(sectionName)

				#Create sub-propagator
				prop = section.propagator(self.psi, i)
				print "Propagator for rank %i is %s" % (i, prop)
			
				#Apply config to sub propagator
				config.Apply(prop)
				section.Apply(prop)
				
				self.SubPropagators.append(prop)
			
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		self.DistributionList = []
		self.DistributionList.append( (distrModel.GetDistribution().copy(), array(self.psi.GetData().shape)) )
		#if len(self.DistributionList[0][0]) > 1: 
		#	raise "Does not support more than 1D proc grid"
		self.RankDistributionMap = zeros(len(self.SubPropagators)+1, dtype=int)

		for prop in self.SubPropagators:
			curRank = prop.TransformRank

			#parallelization
			if not IsSingleProc():
				if distrModel.IsDistributedRank(curRank) and not prop.SupportsParallelPropagation():
					#Get the current distribution	
					distribIndex = self.RankDistributionMap[curRank]
					startDistrib, startShape = self.DistributionList[distribIndex] 

					#Find the next distribution
					finalDistrib = GetAnotherDistribution2(startDistrib, curRank, self.Rank)
		
					#Find shape of the next distribution
					transpose = distrModel.GetTranspose()
					fullShape = self.psi.GetRepresentation().GetFullShape()
					finalShape = transpose.CreateDistributedShape(fullShape, finalDistrib)

					#Update DistributionList
					self.DistributionList.append( (finalDistrib, finalShape) )
					self.RankDistributionMap[curRank+1:] = distribIndex+1

					#Transform into the finalDistrib
					self.Transpose(curRank, self.TransposeForward, self.psi)

			#print "(fwd) currank = %i, distr = %s" % (curRank, self.psi.GetRepresentation().GetDistributedModel().GetDistribution())
			prop.SetupStep(dt/2.)

		#print "DistributionList =", self.DistributionList
		#print "RankDistributionMap =", self.RankDistributionMap

		self.SetupPotential(dt)

		for prop in reversed(self.SubPropagators):
			curRank = prop.TransformRank
			#print "(bwd) currank = %i, distr = %s" % (curRank, self.psi.GetRepresentation().GetDistributedModel().GetDistribution())
			prop.SetupStepConjugate(dt/2.)

			#parallelization:
			if not IsSingleProc():
				self.Transpose(curRank, self.TransposeBackward, self.psi)

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		#first halfstep
		for prop in self.SubPropagators:
			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				self.Transpose(curRank, self.TransposeForward, srcPsi)
				self.Transpose(curRank, self.TransposeForward, destPsi)

			#advance step
			prop.MultiplyHamiltonian(srcPsi, destPsi, t, dt/2.)

		#the potential is in the middle, so it does one full step
		self.MultiplyPotential(srcPsi, destPsi, t, dt)

		#second halfstep
		for prop in reversed(self.SubPropagators):
			#advance step
			prop.MultiplyHamiltonianConjugate(srcPsi, destPsi, t, dt/2.)
			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				self.Transpose(curRank, self.TransposeBackward, srcPsi)
				self.Transpose(curRank, self.TransposeBackward, destPsi)

	
	def CalculatePotentialExpectationValue(self, tmpPsi, potential, t, dt):
		"""
		Calculate expectation value of a grid potential
		"""
		#Callback function, multiplies grid potential on psi
		def MultiplyGridPotential():
			tmpPsi.Clear()
			potential.MultiplyPotential(self.psi, tmpPsi, t, dt)

		#Perform the grid potential multiplication
		self.PerformGridOperation(MultiplyGridPotential, [tmpPsi, self.psi])

		#inner product yields expectation value
		return self.psi.InnerProduct(tmpPsi)

	def AdvanceStep(self, t, dt):
		#Advance one step using strang splitting

		#first halfstep
		for prop in self.SubPropagators:
			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				self.Transpose(curRank, self.TransposeForward, self.psi)

			#advance step
			prop.AdvanceStep(t, dt/2.)

		#the potential is in the middle, so it does one full step
		self.ApplyPotential(t, dt)

		#second halfstep
		for prop in reversed(self.SubPropagators):
			#advance step
			prop.AdvanceStepConjugate(t, dt/2.)
			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				self.Transpose(curRank, self.TransposeBackward, self.psi)

	def Transpose(self, curRank, direction, psi):
		distrModel = psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			nextRank = curRank + 1
			#if curRank < self.Rank-1: nextRank = curRank + 1
			#else: nextRank = curRank - 1
			
			"""
                  0     1    2 
			0)	  D     -    -
			1)    -     D    -
			2)    D     -    -

			curRank = 0
				0->1, do rank0

			curRank = 1
				1->2, do rank1

			curRank = 2
				2->1, do rank2

			do potential

			curRank = 2
				do rank2, 1->2

			curRank = 1
			    do rank1, 2->1

			curRank = 0
				do rank0, 1->0
				
			"""

			#If we're going backward, we should do the
			#inverse transpose in order to ensure that the
			#wavefunction distributed the same way as when going 
			#forward
			if direction == self.TransposeBackward:
				curRank, nextRank = nextRank, curRank

			#Look up the next distribution
			curDistribIndex = self.RankDistributionMap[curRank]
			curDistrib, curShape = self.DistributionList[curDistribIndex]
			nextDistribIndex = self.RankDistributionMap[nextRank]
			nextDistrib, nextShape = self.DistributionList[nextDistribIndex]
			fullShape = self.psi.GetRepresentation().GetFullShape()
			nextShape = distrModel.GetTranspose().CreateDistributedShape(fullShape, nextDistrib)
			#print "rank %i -> %i, distrib %s -> %s, shape %s -> %s" % (curRank, nextRank, curDistrib, nextDistrib, curShape, nextShape)

			#If only transpose if the distribution has changed
			if curDistribIndex != nextDistribIndex:
				#Get transpose buffer
				transposeBufferName = psi.GetAvailableDataBufferName(nextShape)
				if transposeBufferName == -1:
					transposeBufferName = psi.AllocateData(nextShape)
		
				#Change Distribution
				#print "ProcId = %i, direction = %i, curRank = %i (in)" % (ProcId, direction, curRank)
				distrModel.ChangeDistribution(psi, nextDistrib, transposeBufferName)
				#distrModel.GlobalBarrier()
				#print "ProcId = %i, direction = %i, curRank = %i (out)" % (ProcId, direction, curRank)
		
				#HACK: Update all sub-representation with the new distribution
				for i in range(self.Rank):
					subRepr = psi.GetRepresentation().GetRepresentation(i)
					subRepr.GetDistributedModel().SetDistribution(nextDistrib)

	def GetBasisFunction(self, rank, basisIndex):
		prop = self.SubPropagators[rank]
		return prop.GetBasisFunction(rank, basisIndex)


		
	def PerformGridOperation(self, gridFunction, wavefunctionList):
		"""
		Perform a grid operation, as defined by the function 'gridFunction'.
		The wavefunction is transformed to the grid representation, then 
		'gridFunction' is applied, and finally the wavefunction is transformed 
		back to the basis representation.
		"""

		#Transform to grid
		for prop in self.SubPropagators:
			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				for psi in wavefunctionList:
					self.Transpose(curRank, self.TransposeForward, psi)

			#transform step
			for psi in wavefunctionList:
				prop.InverseTransform(psi)

		#Perform grid operation
		gridFunction()

		#Transform to basis
		for prop in reversed(self.SubPropagators):
			#transform step
			for psi in wavefunctionList:
				prop.ForwardTransform(psi)

			#parallelization
			curRank = prop.TransformRank
			if not IsSingleProc():
				for psi in wavefunctionList:
					self.Transpose(curRank, self.TransposeBackward, psi)

