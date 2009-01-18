"""
This file contains python functions for all potentials. 
A potential is usually based on a C++ class implementing the GetPotentialValue 
method. GetPotentialValue returns the potential value in a given coordinate.

Then this class is used as a parameter for a PotentialEvaluator class. The 
potential evaluator knows how to loop over all grid points in a wavefunction, 
and also knows how to perform potential exponentiation, mulitply the potential
with a wavefunction, etc. 

The wrappers in this file create a uniform interface to python for all kinds of
potentials. Apart from special types of potentials like the CrankNicholsonPotential,
there are two main potential types: Dynamic and Static.

Dynamic potentials can, as the name suggests, change value during the course of
propagation, while static potentials remain the same for all values of t. 
Moreover, currently the optimizations made in the static potentials prevents
changing the timestep during the course of propagation. To change the timestep, 
the user must call the UpdateStaticPotential method. This may be changed in the
future.
"""

def CreatePotentialInstance(className, rank, evaluatorPrefix, potentialRank=None):
	"""
	Create an instance of a PotentialEvaluator. Tries to guess
	the python name of the class from a descriptive classname.
	className can be either the full name of the potential except the rank, 
	as specified in the corresponding pyste file.
	Optionally, classname can be just the name of the classname of the dynamic potential
	(without the "DynamicPotentialEvaluator_..."), if it was not renamed in the pyste
	file. 
	"""

	#Create globals
	if potentialRank == None:
		potentialRank = rank

	glob = dict(ProjectNamespace)
	glob.update(globals())
	
	potential = None
	longClassName = evaluatorPrefix + "_" + className + "_" + str(potentialRank) + "_" + str(rank)
	try:
		potential = eval(longClassName + "()", glob)
	except: pass
	
	if potential == None:
		shortClassName = evaluatorPrefix + className + "_" + str(rank)
		try:
			potential = eval(shortClassName + "()", glob)
		except: pass

	if potential == None:
		shortClassName = className + "_" + str(rank)
		try:
			potential = eval(shortClassName + "()", glob)
		except: pass

	if potential == None:
		try:
			potential = eval(className + "()", glob)
		except: pass
			
	if potential == None:
		raise "Unknown potential", className 

	return potential




#Potential Wrapper interface
class PotentialWrapper:
	def ApplyConfigSection(self, configSection):
		"""
		ApplyConfigSection gives the potential the possibiltiy to read properties
		from the associated config section
		"""
		self.ConfigSection = configSection
		self.Type = configSection.type
	
	def SetupStep(self, timestep):
		"""
		SetupStep will be called once to set up the potential. All memory allocations
		and similar one-time actions should be performed here
		"""
		pass
		
	def AdvanceStep(self, t, timestep):
		"""
		Advances the solution according to the potential. That is, applies the exponential
		of the potential V. This is used by the split-step based propagators

		psi = exp(- i dt V ) psi
		"""
		pass

	def MultiplyPotential(self, srcPsi, destPsi, t, timestep):
		"""
		Performs "Matrix"-vector multiplication between the potential and the wavefunction. 
		For diagonal potentials this is simply an element-wise multiplication. This
		is used by ODE-based and krylov-based propagators.

		destPsi += V srcPsi
		"""
		raise "MultiplyPotential is not implemented for class %s" % (self.__class__)
	
	def GetExpectationValue(self, t, timestep):
		"""
		Calculates the expectation value of the wavefunction for the given potential.
		This can also be done with MultiplyPotential, but by using GetExpectationValue, 
		it can be done faster and with less memory requirements.
		"""
		raise "GetExpectationValue is not implemented for class %s" % (self.__class__)


class StaticPotentialWrapper(PotentialWrapper):
	"""
	Wrapper for a static potential. The potential is set up during
	SetupStep, and never changed. Timestep must be fixed during 
	propagation.
	"""

	def __init__(self, psi):
		self.Potential = CreateInstanceRank("core.StaticPotential", psi.GetRank())
		self.psi = psi
		self.Storage = -1

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
	
	def SetupStep(self, timeStep):
		#Determine storage model
		print "Setting up static potential %s" % self
		self.Storage = self.Potential.StorageModel.StorageExpValue
		if hasattr(self.ConfigSection, "storage_model"):
			self.Storage = self.Potential.StorageModel(self.ConfigSection.storage_model)
			print "Using static potential storage model %s" % self.Storage

		#Get time function if defined; check for correct storage model
		self.HasTimeFunction = False
		if hasattr(self.ConfigSection, "time_function"):
			if self.Storage == self.Potential.StorageModel.StorageExpValue:
				raise RunTimeException("Cannot have both time_function and StorageExpValue!")
			else:
				self.TimeFunction = self.ConfigSection.time_function
				self.HasTimeFunction = True

		self.Potential.InitializePotential(self.psi, self.Storage)

		if hasattr(self.ConfigSection, "grid_function"):
			func = self.ConfigSection.grid_function
			updateFunc = eval("core.SetPotentialFromGridFunction_" + str(self.psi.GetRank()))
			updateFunc(self.Potential, timeStep, self.psi, self.psi.GetRepresentation(), func, self.ConfigSection)
	
		elif hasattr(self.ConfigSection, "function"):
			func = self.ConfigSection.function
			potentialData = self.Potential.GetPotentialData()
			func(self.psi, self.ConfigSection, potentialData)
			if self.Storage == self.Potential.StorageModel.StorageExpValue:	
				potentialData[:] = exp(- 1.0j * timeStep * potentialData)

		elif hasattr(self.ConfigSection, "classname"):
			evaluatorPrefix = "core.DynamicPotentialEvaluator"
			potentialEvaluator = CreatePotentialInstance(self.ConfigSection.classname, self.psi.GetRank(), evaluatorPrefix)
			self.ConfigSection.Apply(potentialEvaluator)
			self.PotentialEvaluator = potentialEvaluator

			#update the potential
			self.UpdateStaticPotential(0.0, timeStep)

		else:
			raise "Invalid potential config. Must specify either 'classname' or 'function'"

	def UpdateStaticPotential(self, t, dt):
		"""
		Updates the potential using the PotentialEvaluator created during SetupStep, to use
		the given time and timestep.

		If the potentialevaulator has a function called UpdateStaticPotential, this will be used
		to create the potential. Otherwise, the wavefunction will be set to 1, the potential will be
		applied, and the resulting wavefunction will be used as potential-multiplier.
		"""

		pot = self.PotentialEvaluator
		if hasattr(pot, "UpdateStaticPotential"):
			pot.UpdateStaticPotential(self.Potential, self.psi, dt, 0.0, self.Storage)
		else:
			if self.Storage != StaticStorageModel.StorageExpValue:
				raise Exception("PotentialEvaluators without UpdateStaticPotential must use storage model StorageExpValue")
		
			self.psi.GetData()[:] = 1
			if hasattr(pot, "SetupStep"):
				pot.SetupStep(self.psi, dt)
			pot.ApplyPotential(self.psi, dt, t)
			self.Potential.GetPotentialData()[:] = self.psi.GetData()[:]
				
	def AdvanceStep(self, t, dt):
		scaling = 1.0
		if self.HasTimeFunction:
			scaling = complex(self.TimeFunction(self.ConfigSection, t))
		self.Potential.ApplyPotential(self.psi, dt, scaling)

	def MultiplyPotential(self, srcPsi, destPsi, t, dt):
		if self.Storage == StaticStorageModel.StorageExpValue:
			self.PotentialEvaluator.MultiplyPotential(srcPsi, destPsi, dt, t)
		#elif self.Storage == -1:
		#	raise Exception("what! %s" % self)
		else:
			scaling = 1.0
			if self.HasTimeFunction:
				scaling = self.TimeFunction(self.ConfigSection, t)
			self.Potential.MultiplyPotential(srcPsi, destPsi, dt, scaling)

	def GetExpectationValue(self, t, dt):
		return self.PotentialEvaluator.CalculateExpectationValue(self.psi, dt, t)

	def GetPotential(self, dt):
		if self.Potential.UseStorageValue():
			return real(self.Potential.GetPotentialData()).copy()
		else:
			return real(log(self.Potential.GetPotentialDataExp()) / (- 1.0j * dt))

class DynamicPotentialWrapper(PotentialWrapper):
	"""
	Wrapper for dynamic potentials. Dynamic potentials are classes implemented in 
	C++, and wrapped by DynamicPotentialEvaluator to make evaluation simple.
	Dynamic potentials have no large memory footprint, and can change each timestep,
	and also timeStep may change any time.
	"""
	
	def __init__(self, psi):
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()

		evaluatorPrefix = "core.DynamicPotentialEvaluator"
		self.Evaluator = CreatePotentialInstance(configSection.classname, rank, evaluatorPrefix)
		configSection.Apply(self.Evaluator)
	
	def AdvanceStep(self, t, dt):
		self.Evaluator.ApplyPotential(self.psi, dt, t)

	def MultiplyPotential(self, srcPsi, destPsi, t, dt):
		self.Evaluator.MultiplyPotential(srcPsi, destPsi, dt, t)

	def GetExpectationValue(self, t, dt):
		return self.Evaluator.CalculateExpectationValue(self.psi, dt, t)


class RankOnePotentialWrapper(PotentialWrapper):
	"""
	Wrapper for dynamic and static potentials only dependent on the grid in one rank. 
	In this case, a number of optimization techniques are used to speed up the potential
	evaluation
	"""

	def __init__(self, psi):
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()
	
		evaluatorPrefix = "core.RankOnePotentialEvaluator"
		self.Evaluator = CreatePotentialInstance(configSection.classname, rank, evaluatorPrefix, potentialRank=1)
		configSection.Apply(self.Evaluator)
	
	def AdvanceStep(self, t, dt):
		self.Evaluator.ApplyPotential(self.psi, dt, t)

	def MultiplyPotential(self, srcPsi, destPsi, t, dt):
		self.Evaluator.MultiplyPotential(srcPsi, destPsi, dt, t)

	def GetExpectationValue(self, t, dt):
		return self.Evaluator.CalculateExpectationValue(self.psi, dt, t)


class FiniteDifferencePotentialWrapper(PotentialWrapper):
	def __init__(self, psi):
		self.psi = psi

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()

		evaluatorPrefix = "core.ExponentialFiniteDifferenceEvaluator"
		self.Evaluator = CreatePotentialInstance(configSection.classname, rank, evaluatorPrefix)
		configSection.Apply(self.Evaluator)
	
	def AdvanceStep(self, t, dt):
		self.Evaluator.UpdateWavefunction(self.psi, t, dt/2, 0)
		self.Evaluator.UpdateWavefunction(self.psi, t, dt, 1)
		self.Evaluator.UpdateWavefunction(self.psi, t, dt/2, 0)


class CrankNicholsonPotentialWrapper(PotentialWrapper):
	def __init__(self, psi):
		self.psi = psi
		self.EnergyDataName = -1

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()

		evaluatorPrefix = "core.CrankNicholsonEvaluator"
		self.Evaluator = CreatePotentialInstance(configSection.classname, rank, evaluatorPrefix)
		configSection.Apply(self.Evaluator)
	
	def AdvanceStep(self, t, dt):
		self.Evaluator.UpdateWavefunction(self.psi, t, dt)

	def MultiplyPotential(self, srcPsi, dstPsi, t, dt):
		self.Evaluator.MultiplyOperator(srcPsi, dstPsi, t, dt)


class SparseMatrixParticle(tables.IsDescription):
	RowIndex = tables.Int32Col(pos=1)
	ColIndex = tables.Int32Col(pos=2)
	MatrixElement = tables.ComplexCol(itemsize=16, pos=3)

class MatrixType:
	Sparse = "Sparse"
	Dense = "Dense"

class MatrixPotentialWrapper(PotentialWrapper):

	def __init__(self, psi):
		self.psi = psi

	def LoadMatrixPotential(self, conf):
		self.MatrixType = conf.matrix_type 
		
		if hasattr(conf, "matrix_function"):
			self.LoadMatrixPotentialFunction(conf)

		elif hasattr(conf, "filename"):
			self.LoadMatrixPotentialFile(conf)

		else:
			raise NotImplementedException("Missing load-matrix information in Potential")


	def LoadMatrixPotentialFunction(self, conf):
		func = conf.matrix_function

		if self.MatrixType == MatrixType.Sparse:
			row, col, matrixElement = func(self.psi, conf)
			self.Potential.SetMatrixData(row, col, matrixElement)

		if self.MatrixType == MatrixType.Dense:
			data = func(self.psi, conf)
			self.Potential.SetMatrixData(data, self.MatrixRowRank, self.MatrixColRank)

	def LoadMatrixPotentialFile(self, conf):
		filename = conf.filename
		dataset = conf.dataset
		file = tables.openFile(filename, "r")
		try:
			node = file.getNode(dataset)
			if self.MatrixType == MatrixType.Sparse:
				row = node.cols.RowIndex[:]
				col = node.cols.ColIndex[:]
				matrixElement = node.cols.MatrixElement[:]
				self.Potential.SetMatrixData(row, col, matrixElement)

			if self.MatrixType == MatrixType.Dense:
				data = node[:]
				self.Potential.SetMatrixData(data, self.MatrixRowRank, self.MatrixColRank)

		finally:
			file.close()

	def ApplyConfigSection(self, configSection):
		PotentialWrapper.ApplyConfigSection(self, configSection)
		rank = self.psi.GetRank()

		if configSection.matrix_type == MatrixType.Sparse:
			self.Potential = core.SparseMatrixPotentialEvaluator()

		elif configSection.matrix_type == MatrixType.Dense:
			self.Potential = CreateInstanceRank("core.DenseMatrixPotentialEvaluator", rank)
			matrixColRank = None
			matrixRowRank = None
			if rank == 1:
				matrixRowRank = 0
				matrixColRank = 1
			if hasattr(configSection, "matrix_column_rank"):
				matrixColRank = configSection.matrix_column_rank
			if hasattr(configSection, "matrix_row_rank"):
				matrixRowRank = configSection.matrix_row_rank
			if matrixColRank == None or matrixRowRank == None:
				raise Exception("Must specify matrix_column_rank and matrix_row_rank for DenseMatrixPotential")
			self.MatrixColRank = matrixColRank
			self.MatrixRowRank = matrixRowRank

		else:
			raise NotImplementedException("Invalid MatrixType %s" % configSection.MatrixType)
		
		configSection.Apply(self.Potential)
		self.LoadMatrixPotential(configSection)
		self.TimeFunction = configSection.time_function

	def AdvanceStep(self, t, dt):
		self.Potential.ApplyPotential(self.psi, dt, self.GetTimeValue(t))

	def MultiplyPotential(self, srcPsi, destPsi, t, dt):
		self.Potential.MultiplyPotential(srcPsi, destPsi, self.GetTimeValue(t))

	def GetExpectationValue(self, t, dt):
		return self.GetTimeValue(t) * self.Potential.CalculateExpectationValue(self.psi)

	def GetTimeValue(self, t):
		return self.TimeFunction(self.ConfigSection, t)

