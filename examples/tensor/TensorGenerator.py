class UnsupportedGeometryException(Exception):
	def __init__(self, string):
		Exception.__init__(self, string)

#------------------------------------------------------------------------------------
#                       Interfaces
#------------------------------------------------------------------------------------

class GeometryInfoBase(object):
	"""
	base for all Geometry information objects
	a GeometryInfo object gives information about how the geometry behaves
	"""

	def UseGridRepresentation(self):
		raise NotImplementedException()

	def GetBasisPairCount(self):
		raise NotImplementedException()

	def GetBasisPairs(self):
		raise NotImplementedException()

	def GetOverlapMatrix(self):
		raise NotImplementedException()



class BasisfunctionBase(object):
	"""
	Abstract class which bases all basis function classes.
	Its interface provides 
	"""

	def ApplyConfigSection(self, config):
		raise NotImplementedException()

	def GetGridRepresentation(self):
		"""
		Creates a one dimensional representation for this rank
		in the grid representation
		"""
		raise NotImplementedException()

	def GetBasisRepresentation(self):
		"""
		Creates a one dimensional representation for this rank
		in the basis representation
		"""
		raise NotImplementedException()

	def GetGeometryInfo(self, geometryName):
		raise NotImplementedException()

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		raise NotImplementedException


#------------------------------------------------------------------------------------
#                       BSpline
#------------------------------------------------------------------------------------

class GeometryInfoBSplineDense(GeometryInfoBase):
	"""
	Geometry information for BSpline geometries
	"""
	def __init__(self, bsplineObject):
		#Set member variables 
		self.BSplineObject = bsplineObject	

	def UseGridRepresentation(self):
		return True
	
	def GetBasisPairCount(self):
		return self.BSplineObject.NumberOfBSplines**2
		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()
		
		pairs = zeros((count, 2), dtype=int)
		index = 0
		for i in xrange(self.BSplineObject.NumberOfBSplines):
			for j in xrange(self.BSplineObject.NumberOfBSplines):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		return pairs

class GeometryInfoBSplineBanded(GeometryInfoBase):
	"""
	Geometry information for BSpline geometries
	"""
	def __init__(self, bsplineObject):
		#Set member variables 
		self.BSplineObject = bsplineObject	

	def UseGridRepresentation(self):
		return True
	
	def GetBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return (N * k - (k) * (k-1) / 2) * 2  
		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(max(0, i-k), min(i+k, N)):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		return pairs



from pyprop import InitBSpline

class BasisfunctionBSpline(BasisfunctionBase):
	"""
	Basisfunction class for BSplines
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis( InitBSpline(configSection) )

	def SetupBasis(self, bsplineObject):
		self.BSplineObject = bsplineObject

	def GetGridRepresentation(self):
		repr = pyprop.core.BSplineGridRepresentation()
		repr.SetupRepresentation(self.BSplineObject)
		return repr

	def GetBasisRepresentation(self):
		repr = pyprop.core.BSplineRepresentation()
		repr.SetupRepresentation(self.BSplineObject)
		return repr

	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			raise Exception("TODO: Implement Identity!")
		elif geom == "banded":
			return GeometryInfoBSplineBanded(self.BSplineObject)
		elif geom == "dense":
			return GeometryInfoBSplineDense(self.BSplineObject)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionBSpline" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		#TODO: Implement for general rank

		pairs = geometryInfo.GetBasisPairs()
		RepresentPotentialInBasisBSpline(self.BSplineObject, source, dest, pairs, rank, differentiation)



#------------------------------------------------------------------------------------
#                       Reduced Spherical Harmonic
#------------------------------------------------------------------------------------

class GeometryInfoReducedSphHarmDense(GeometryInfoBase):
	"""
	Geometry information for Reduced spherical harmonic geometries
	with no symmetries
	"""
	def __init__(self, sphericalHarmonicObject):
		#Set member variables 
		self.SphericalHarmonicObject = sphericalHarmonicObject	

	def UseGridRepresentation(self):
		return True
	
	def GetBasisPairCount(self):
		return (self.SphericalHarmonicObject .GetLMax()+1)**2
		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()
		basisCount = sqrt(count)
		
		pairs = zeros((count, 2), dtype=int)
		index = 0
		for i in xrange(basisCount):
			for j in xrange(basisCount):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		return pairs

class GeometryInfoReducedSphHarmSelectionRule(GeometryInfoBase):
	"""
	Geometry information for Reduced spherical harmonic geometries
	with delta l = +- 1 symetry
	"""
	def __init__(self, sphericalHarmonicObject):
		#Set member variables 
		self.SphericalHarmonicObject = sphericalHarmonicObject	

	def UseGridRepresentation(self):
		return True
	
	def GetBasisPairCount(self):
		return self.SphericalHarmonicObject.GetLMax()*2
		
	def GetBasisPairs(self):
		count = self.SphericalHarmonicObject.GetLMax()
		
		pairs = zeros((count*2, 2), dtype=int)
		index = 0
		for i in xrange(count):
			pairs[index, 0] = i
			pairs[index, 1] = i-1
			index+=1
			pairs[index, 0] = i
			pairs[index, 1] = i+1
			index+=1

		return pairs


class GeometryInfoReducedSphHarmDiagonal(GeometryInfoBase):
	"""
	Geometry information for Reduced spherical harmonic geometries
	with \delta l = 0 symmetry. The potential is specified in the 
	l-basis
	"""
	def __init__(self, sphericalHarmonicObject):
		#Set member variables 
		self.SphericalHarmonicObject = sphericalHarmonicObject	

	def UseGridRepresentation(self):
		return False
	
	def GetBasisPairCount(self):
		return self.SphericalHarmonicObject.GetLMax()+1
		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()
		
		pairs = zeros((count, 2), dtype=int)
		index = 0
		for i in xrange(count):
			pairs[index, 0] = i
			pairs[index, 1] = i
			index+=1

		return pairs



class BasisfunctionReducedSphericalHarmonic(BasisfunctionBase):
	"""
	Basisfunction class for BSplines
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis(configSection.lmax)

	def SetupBasis(self, lmax):
		self.LMax = lmax
		self.SphericalHarmonicObject = pyprop.core.ReducedSphericalTools()
		self.SphericalHarmonicObject.Initialize(self.LMax)

	def GetGridRepresentation(self):
		repr = pyprop.core.ThetaRepresentation()
		repr.SetupRepresentation(self.LMax)
		return repr

	def GetBasisRepresentation(self):
		repr = pyprop.core.ReducedSphericalHarmonicRepresentation()
		repr.SetupRepresentation(self.LMax)
		return repr

	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			raise Exception("TODO: Implement Identity!")
		if geom == "diagonal":
			return GeometryInfoReducedSphHarmDiagonal(self.SphericalHarmonicObject)
		elif geom == "dense":
			return GeometryInfoReducedSphHarmDense(self.SphericalHarmonicObject)
		elif geom == "dipoleselectionrule":
			return GeometryInfoReducedSphHarmSelectionRule(self.SphericalHarmonicObject)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionReducedSpherical" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		#TODO: Implement for general rank

		assert(differentiation==0, "Differentiation is not supported on Spherical Harmonics yet")

		pairs = geometryInfo.GetBasisPairs()
		RepresentPotentialInBasisReducedSphericalHarmonic(self.SphericalHarmonicObject, source, dest, pairs, rank, differentiation)




#------------------------------------------------------------------------------------
#                       TensorPotentialGenerator main
#------------------------------------------------------------------------------------

from pyprop import CreateInstanceRank
from pyprop import CreateDistribution
from pyprop import CreateWavefunctionInstance
from pyprop import CreatePotentialInstance

def CreateBasisFromRepresentation(representation):
	if representation.__class__ == pyprop.core.BSplineRepresentation:
		basis = BasisfunctionBSpline()
		basis.SetupBasis(representation.GetBSplineObject())
	
	elif representation.__class__ == pyprop.core.ReducedSphericalHarmonicRepresentation:
		basis = BasisfunctionReducedSphericalHarmonic()
		basis.SetupBasis(representation.Range.MaxL)

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
		#This will be used to hold the potential
		#and representation in order to use the existing potential evaluator
		#framework in pyprop
		psi = CreateWavefunctionInstance(repr)

		#3) Create potential evaluator
		evaluatorPrefix = "core.DynamicPotentialEvaluator"
		classname = configSection.classname
		potentialEvaluator = CreatePotentialInstance(classname, self.Rank, evaluatorPrefix)
		configSection.Apply(potentialEvaluator)

		#4) Evaluate the potential on the grid
		data = potentialEvaluator.UpdatePotentialData(psi.GetData(), psi, 0, 0)

		#for parallelization
		#We only support the first rank initially distributed
		assert(repr.GetDistributedModel().GetDistribution()[0] == 0)
		assert(len(repr.GetDistributedModel().GetDistribution()) == 1)
		fullShape = repr.GetFullShape()
		if repr.GetDistributedModel().IsSingleProc():
			distribution = []
		else:
			distribution = list(repr.GetDistributedModel().GetDistribution())
		transpose = repr.GetDistributedModel().GetTranspose() 

		#5) Represent the potential in the bases
		source = psi.GetData()
		for rank in range(self.Rank-1,-1,-1):
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
					print "Rank %i is distributed" % rank
					assert(rank != self.Rank-1)
					assert(not rank+1 in distribution)

					#Create new distribution
					distribIndex = distribution.index(rank)
					newDistribution = list(distribution)
					newDistribution[distribIndex] = rank+1

					#Create shape of transposed function and allocate dest buffer
					transposedShape = transpose.CreateDistributedShape(fullShape, array(newDistribution, dtype=int))
					dest = zeros(transposedShape, dtype=complex)

					#Transpose
					transpose.Transpose(fullShape, source, array(distribution, int), dest, array(newDistribution, int))

					#Use the new buffer
					source = dest
					distribution = newDistribution

				#Calculate dest shape
				destShape = array(source.shape)
				destShape[rank] = geometryInfo.GetBasisPairCount()

				#Update full shape
				fullShape[rank] = geometryInfo.GetBasisPairCount()
			
				#Allocate the destination array
				dest = zeros(destShape, dtype=complex)

				#Represent this rank in the basis
				basis.RepresentPotentialInBasis(source, dest, rank, geometryInfo, differentiation) 

				#Use the destination from this rank as the source to the next
				source = dest

		#done!
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



#------------------------------------------------------------------------------------
#                       Propagation classes
#------------------------------------------------------------------------------------

from pyprop import PotentialWrapper

class TensorPotential(PotentialWrapper):
	"""
	Potential wrapper for TensorPotentials. See PotentialWrapper for more information on the
	PotentialWrapper interface

	A TensorPotential is potential expressed in the basis functions of the wavefunction.
	This makes it possible to apply the potential without transforming the wavefunction to the
	grid. The downside is that potentials are usually diagonal in the grid basis, and not in the
	basisfunction basis. This means that we must do a matrix vector product to apply the potential.

	in 1D
	out(i) = sum_j V(i, j) psi(j) 
	in 2D
	out(i', j') = sum_{i,j} V(i', i, j', j) psi(i', j')

	in the tensor potential V has the same rank as the wavefunction, where i',i and j',j are compressed
	into one rank each. This is done in order to easier be able to exploit symmetries in the potentials and
	bandedness in the basises (such as the dipole selection rule V = 0 for l' != l +/- 1).

	
	"""

	def __init__(self, psi):
		self.GeometryList = None
		self.PotentialData = None
		self.Name = None
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		#Check wheter this is a time dependent potential
		self.IsTimeDependent = False
		if hasattr(configSection, "time_function"):
			self.IsTimeDependent = True
			self.TimeFunction = lambda t: configSection.time_function(configSection, t)

	def SetupStep(self, timestep):
		self.BasisPairs = [self.GeometryList[i].GetBasisPairs() for i in range(self.psi.GetRank())]
		
	def AdvanceStep(self, t, timestep):
		raise NotImplementedException("TensorPotentials can not be exponentiated directly")

	def MultiplyPotential(self, destPsi, t, timestep):
		rank = self.psi.GetRank()
		
		source = self.psi.GetData()
		dest = destPsi.GetData()

		timeScaling = 1.0
		if self.IsTimeDependent:
			timeScaling = self.TimeFunction(t)
		
		#TODO: Implement support for parallelization. 
		if rank == 1:
			pairs = self.BasisPairs[0]
			MultiplyTensorPotential_1(self.PotentialData, timeScaling, pairs, source, dest)

		elif rank == 2:
			pairs0 = self.BasisPairs[0]
			pairs1 = self.BasisPairs[1]

			MultiplyTensorPotential_2(self.PotentialData, timeScaling, pairs0, pairs1, source, dest)

		else:
			raise NotImplementedException("Only rank=1 and rank=2 is currently implemented for TensorPotential")
	
	def GetExpectationValue(self, t, timestep):
		raise NotImplementedException("GetExpectationValue is not implemented for class %s" % (self.__class__))


from pyprop import PropagatorBase

class BasisPropagator(PropagatorBase):
	"""
	Propagator that does not transform the wavefunction to the grid basis, but rather applies the potentials
	(TensorPotentials) directly in the basis.

	BasisPropagator assumes that the wavefunction has a CombinedRepresentation.

	AdvanceStep is not supported by BasisPropagator. Rather, one should use a propagator utilizing MultplyHamiltonian
	such as the ODE and Krylov propagator.
	"""

	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()
	
		#Create a tensor potential generator 		
		#We will use this to construct tensor potentials from potentials specified on the grid
		self.TensorPotentialGenerator = TensorPotentialGenerator(representation = self.psi.GetRepresentation())

		#We need the temp array for solving overlap matrix eqns
		self.TempData = zeros(self.psi.GetData().shape, dtype=complex)

	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)

		self.GeneratePotentials(config)
		self.ConsolidatePotentials()


	def GeneratePotentials(self, config):
		"""
		Genereate TensorPotentials from potentials specified on the grid in the
		configuration file
		"""

		#Potentials we should create on the fly
		if hasattr(config.Propagation, "grid_potential_list"):
			potentials = config.Propagation.grid_potential_list
			generator = self.TensorPotentialGenerator

			for potentialName in potentials:
				#Find the corresponding config section
				configSection = config.GetSection(potentialName)

				#Use TensorPotentialGenerator to construct potential in basis
				geometryList = generator.GetGeometryList(configSection)
				potentialData = generator.GeneratePotential(configSection)

				#Create PotentialWrapper for TensorPotential
				potential = TensorPotential(self.psi)
				configSection.Apply(potential)
				potential.GeometryList = geometryList
				potential.PotentialData = potentialData
				potential.Name = potentialName

				#Add potential to potential list
				self.PotentialList.append(potential)

	def ConsolidatePotentials(self):
		"""
		Try to consolidate potentials having the same geometry into
		one potential to save evaulations. Potentials containing a time
		dependent part needs to be treated separately, and thus is not considered
		"""
		
		print "Consolidating similar potentials: "
		print "Starting with potentials:"
		for pot in self.PotentialList:
			print "    %s" % pot.Name

		#only non timedependent potentials are considered
		potentials = [pot for pot in self.PotentialList if not pot.IsTimeDependent]
		removePotentials = []

		#Consolidate til we're at the last potential
		i = 0
		while(i < len(potentials)-1):
			curPot = potentials[i]

			#Loop over all potentials after curPot
			for otherPot in list(potentials[i+1:]):
				#We can consolidate curPot and otherPot if all the index pairs are the same
				canConsolidate = True
				for rank in range(self.Rank):
					if not numpy.all(curPot.GeometryList[rank].GetBasisPairs() == otherPot.GeometryList[rank].GetBasisPairs()):
						canConsolidate = False
						break

				#Add oterPot to curPot
				curPot.PotentialData[:] += otherPot.PotentialData
				curPot.Name += "+" + otherPot.Name
				potentials.remove(otherPot)
				removePotentials.append(otherPot)

			#Next potential
			i += 1

		#Remove consolidated potentials from the main potential list
		for pot in removePotentials:
			self.PotentialList.remove(pot)

		print "Ended up with potentials:"
		for pot in self.PotentialList:
			print "    %s" % pot.Name


	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		self.__Base.SetupStep(self, dt)

	def MultiplyHamiltonian(self, destPsi, t, dt):
		#Multiply potentials
		self.MultiplyPotential(destPsi, t, dt)

		#Solve for all overlap matrices
		repr = self.psi.GetRepresentation()
		source = destPsi.GetData()
		dest = self.TempData
		hasNonOrthogonalBasis = False
		for i in range(self.Rank):
			isOrthogonal = repr.IsOrthogonalBasis(i)
			if not isOrthogonal:
				#Solve for the overlap matrix for this rank
				overlapMatrix = repr.GetGlobalOverlapMatrix(i)
				self.SolveForOverlapMatrix(source, dest, overlapMatrix, i)
				
				#Use the current dest as the next source, and vice versa
				source, dest = dest, source
				hasNonOrthogonalBasis = True

		#Make sure we end up with the correct array in destPsi
		if hasNonOrthogonalBasis:
			destPsi.GetData()[:] = source
		
	def AdvanceStep(self, t, dt):
		raise NotImplementedException("BasisPropagator does not support AdvanceStep. Use it as a base for explicit propagators")

	def GetBasisFunction(self, rank, basisIndex):
		raise NotImplementedException("Implement GetBasisFunction...")

	def SolveForOverlapMatrix(self, source, dest, overlapMatrix, rank):
		bspline = self.psi.GetRepresentation().GetRepresentation(rank).GetBSplineObject()
		SolveForOverlapMatrix(bspline, source, dest, rank)
	



