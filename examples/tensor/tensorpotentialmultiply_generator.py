
"""

Generator to create matrix multiply functions for TensorPotentials.
This file prints out Fortran90 functions for calculating MultiplyPotential on 
tensor potentials

each rank in a TensorPotential has a storage defined by the GeometryInfo for
that rank. This can typically be dense (the tensor potential is a full N-by-N matrix)
banded (for bspline basises), or some other sparsity definition, such as
|i-j| == 1 for the dipole selection rule in spherical harmonics. 

it is possible to create a general matrix vector product, by listing all the indices
(i,j) for each rank and then do the matrix vector product using these index pairs. in 3D this
would look like

in C++:
/*
 * 2D arrays where row = pair0(i, 0), col = pair0(i, 1) specifies the (row, col) which 
 * potentialData(i, ...) is the matrix element between.
 */
pair0(pairCount0, 2)
pair1(pairCount1, 2)
pair2(pairCount2, 2)

for (int i0=0; i0<pairCount0; i0++)
{
	int row0 = pair0(pairCount0, 0);
	int col0 = pair0(pairCount0, 1);

	for (int i1=0; i1<pairCount1; i1++)
	{
		int row1 = pair1(paircount1, 0);
		int col1 = pair1(paircount1, 1);

		for (int i2=0; i2<pairCount2; i2++)
		{
			int row2 = pair1(paircount2, 0);
			int col2 = pair1(paircount2, 1);

			dest(col0, col1, col2) += potentialData(i0, i1, i2) * source(row0, row1, row2);
		}
	}
}

In Fortran

subroutine func &
    ( &
        potential, potentialExtent0, potentialExtent1, potentialExtent2, &
        scaling, &
        pair0, pair0Extent0, pair0Extent1, &
        pair1, pair1Extent0, pair1Extent1, &
        pair2, pair2Extent0, pair2Extent1, &
        source, sourceExtent0, sourceExtent1, sourceExtent2, &
        dest, destExtent0, destExtent1, destExtent2&
    )
    implicit none
    include "parameters.f"

	complex (kind=dbl), dimension(0:potentialExtent2-1, 0:potentialExtent1-1, 0:potentialExtent0-1), intent(in) :: potential
	integer, intent(in) :: potentialExtent0, potentialExtent1, potentialExtent2

	real (kind=dbl), intent(in) :: scaling

	complex (kind=dbl), dimension(0:pair0Extent1-1, 0:pair0Extent0-1) :: pair0
	integer, intent(in) :: pair0Extent0, pair0Extent1

	complex (kind=dbl), dimension(0:pair1Extent1-1, 0:pair1Extent0-1) :: pair1
	integer, intent(in) :: pair1Extent0, pair1Extent1

	complex (kind=dbl), dimension(0:pair2Extent1-1, 0:pair2Extent0-1) :: pair2
	integer, intent(in) :: pair2Extent0, pair2Extent1

	complex (kind=dbl), dimension(0:sourceExtent2-1, 0:sourceExtent1-1, 0:sourceExtent0-1), intent(in) :: source
	integer, intent(in) :: sourceExtent0, sourceExtent1, sourceExtent2

	complex (kind=dbl), dimension(0:destExtent2-1, 0:destExtent1-1, 0:destExtent0-1) :: dest
	integer, intent(in) :: destExtent0, destExtent1, destExtent2

	integer :: i0, row0, col0
	integer :: i1, row1, col1
	integer :: i2, row2, col2

	do i0=0, potentialExtent0
		row0 = pair0(0, i0)
		col0 = pair0(1, i0)
	
		do i1=0, potentialExtent1
			row1 = pair1(0, i1)
			col1 = pair1(1, i1)

			do i2=0, potentialExtent2
				row2 = pair2(0, i2)
				col2 = pair2(1, i2)

				dest(row2, row1, row0) = dest(row2, row1, row0) + potential(i2, i1, i0) * source(row2, row1, row0)
			enddo
		enddo
	enddo

end subroutine func

Fortran conventions:
-  All arrays are 0-based
- In pararmeter lists, an N-dimensional array is followed by N args named 
  <arrayname>Extent0 through <arrayName>ExtentN. Extent0 refers to the extent
  of LARGEST stride, while ExtentN refers to the extent of LOWEST stride 
  (which will always be 1). 
  This is very intentionally "C-style" ordering even if the arrays will be accessed
  in Col-major ordering

However this require us to store all elements of potentialData, even if most of the potentials we will
be using will be Hermitian. It also will not be able to optimize the inner loop using BLAS if the 
innermost rank is dense or banded. 

By generating a function for every combination of storage type, we can let the storage type define
how it will loop over the indices in potentialData, exploiting all kinds of symmetries and efficient
multiplication methods (ie using BLAS)

A storage type will need to supply the following
1) Initialization code. Run before the loop starts
2) Looping code. Iterating over all rows and cols of this rank
3) Inner loop code. Calling the inner ranks, allowing them to supply their own looping code

/* Begin rant
I would like to use this opportunity to complain over the state of C++ arrays. 
There are currently NO array implementation for C or C++ that gives 
satisfactory performance, while keeping indexing simple. blitz++ does a good
job of creating a simple interface, and for vectorizable operations it can deliver
good performance. However, the performance is STRONGLY dependent on the compiler,
and the reduction operations has troubles on most compilers.

Fortran, on the other hand has a pretty good interface to arrays. Indexing is 
similar to blitz arrays, but as it is a language feature, rather than an external
library, good compilers are able to optimize access in ways not possible for blitz.
The problem with fortran is that anything else that looping over arrays is horrible pain. 
Furthermore, calling fortran functions from other languages are not standarized, so 
in general one cannot create a portable wrapper from fortran to other languages. In 
practice it is not so grim, as fortran functions are post fixed by an underscore on 
most compilers.
End rant */

The solution chosen here, is to use generate fortran kernels doing nothing more than
the computationally intensive loops. Then calling these kernels from C++, where blitz
arrays are used to keep track of array dimensions. The C++ functions are in turned
wrapped in Python for high level low intensity parts of the program.


"""


def PrettyPrintC(str):
	outStr = ""
	lines = [line.strip() for line in str.split("\n")]
	
	indent = 0
	for line in lines:
		if line == "}":
			indent -= 1

		outStr += "\t" * indent	+ line + "\n"

		if line == "{":
			indent += 1

	return outStr



def PrettyPrintFortran(str):
	outStr = ""
	lines = [line.strip() for line in str.split("\n")]
	
	indent = 0
	for line in lines:
		if line.startswith("enddo") or line.startswith("endif") or line.startswith("end subroutine"):
			indent -= 1

		outStr += "\t" * indent	+ line + "\n"

		if line.startswith("if ") or line.startswith("do ") or line.startswith("subroutine"):
			indent += 1

	return outStr


#-----------------------------------------------------------------------------------
#		Fortran code generation routines
#-----------------------------------------------------------------------------------

TypeMapFortranToC = { \
	"integer": "int",\
	"complex (kind=dbl)": "cplx",\
	"real (kind=dbl)": "double",\
	}



def GetFortranArrayParameterList(paramName, rank):
	"""
	Returns a list of parameters required for a rank-dimensional array. This
	is [<paramName>, <paramName>Extent0, ..., <paramName>Extent<rank>]

	i.e. 

	GetFortranArrayParameterList("data", 2)  
	returns
	["data", "dataExtent0", "dataExtent1"]
	"""
	parameterList = [paramName] 
	parameterList += ["%sExtent%i" % (paramName, i) for i in range(rank)]
	return parameterList

def GetFortranArrayDeclaration(paramName, rank, dataType, intent):
	"""
	Returns free form Fortran90 code that contains the declaration of an n-dimensional array

	i.e.
	
	GetFortranArrayDeclaration("data", 2, "complex (kind=dbl)", "inout")
	returns
	'''
	complex, (kind=dbl), dimension(0:dataExtent1-1, 0:dataExtent0,-1), intent(inout) :: data
	integer, intent(in) :: dataExtent1, dataExtent0
	'''

	Note that dataExtent0 is the extent of the rank with the LARGEST stride
	"""
	paramList = GetFortranArrayParameterList(paramName, rank)
	extentList = paramList[1:]
	extentList.reverse()

	dimensionString = ", ".join(["0:%s-1" % extent for extent in extentList])
	extentString = ", ".join(extentList)

	str = """
		%(dataType)s, dimension(%(dimensionString)s), intent(%(intent)s) :: %(paramName)s
		integer, intent(in) :: %(extentString)s
	""" % locals()
	return str


#-----------------------------------------------------------------------------------
#		Generators
#-----------------------------------------------------------------------------------

class StorageGeneratorBase(object):
	def __init__(self, systemRank, curRank, innerGenerator):
		self.SystemRank = systemRank
		self.CurRank = curRank
		self.InnerGenerator = innerGenerator

	def GetParameterList(self):
		"""
		Returns the function parameters that this loop requires
		"""
		raise Exception("Not Implemented")
	
	def GetParameterDeclarationCode(self):
		"""
		"""
		raise Exception("Not Implemented")
	
	def GetInitializationCode(self):
		"""
		Returns the code needed to initialize this loop. All variables
		created should be postfixed by curRank to indicate that the belong to this
		rank
		"""
		raise Exception("NotImplemented")

	def GetLoopingCodeRecursive(self, conjugate):
		"""
		Returns the looping code. All variables created should be postfixed by
		self.CurRank to indicate that they belong to this rank. At the very least, 
		it should define int row<self.CurRank> and int col<self.CurRank>, and an index
		to the potential i<self.CurRank> as these 
		will be referenced by inner loops. 

		This routine should also call self.InnerGenerator.GetLoopingCodeRecursive(conjugate) to get
		the inner loops, and ad insert the code at the appropriate position in it's own 
		looping code
		"""
		raise Exception("NotImplemented")

	def GetIndexString(self, indexPrefix):
		return ", ".join([indexPrefix + str(i) for i in range(self.SystemRank-1,-1,-1)])


#-----------------------------------------------------------------------------------

class StorageGeneratorSimple(StorageGeneratorBase):
	"""
	Generator for simple storage, that is, where the loop over index
	pairs is the way to do matrix multiplication
	"""

	def __init__(self, systemRank, curRank, innerGenerator):
		StorageGeneratorBase.__init__(self, systemRank, curRank, innerGenerator)

	def GetParameterList(self):
		parameterList = [("pair%i" % self.CurRank, "array", 2, "integer")]
		return parameterList

	def GetParameterDeclarationCode(self):
		str = GetFortranArrayDeclaration("pair%i" % self.CurRank, 2, "integer", "in")
		str += """
			integer :: N%(rank)i, i%(rank)i, row%(rank)i, col%(rank)i
		""" % { "rank": self.CurRank }
		return str

	def GetInitializationCode(self):
		str = """
			N%(rank)i = potentialExtent%(rank)i
		""" % { "rank":self.CurRank }
		return str

	def GetLoopingCodeRecursive(self, conjugate):
		str = """
			do i%(rank)i = 0, N%(rank)i-1
				row%(rank)i = pair%(rank)i(0, i%(rank)i);
				col%(rank)i = pair%(rank)i(1, i%(rank)i);
				%(innerLoop)s
			enddo
		""" % { "rank":self.CurRank, "innerLoop": self.GetInnerLoop(conjugate) }
		return str

	def GetInnerLoop(self, conjugate):
		str = ""
		if self.InnerGenerator != None:
			str += self.InnerGenerator.GetLoopingCodeRecursive(conjugate)
		else:
			#We're at the innermost loop
			conjg = ""
			if conjugate:
				conjg = "conjg"
			str += """
				dest(%(rowIndex)s) = dest(%(rowIndex)s) + %(conjg)s(potential(%(potentialIndex)s)) * scaling * source(%(colIndex)s)
			""" % \
				{ \
					"rowIndex": self.GetIndexString("row"), \
					"potentialIndex": self.GetIndexString("i"), \
					"colIndex": self.GetIndexString("col"),  \
					"conjg": conjg, \
				}
		return str



#-----------------------------------------------------------------------------------


class StorageGeneratorHermitian(StorageGeneratorSimple):
	"""
	Generator for simple hermitian storage, that is, where the loop 
	over index pairs is the way to do matrix multiplication, but the 
	index pairs only include the upper part of the matrix (row < col)
	"""

	def __init__(self, systemRank, curRank, innerGenerator):
		StorageGeneratorSimple.__init__(self, systemRank, curRank, innerGenerator)

	def GetLoopingCodeRecursive(self, conjugate):
		str = """
			do i%(rank)i = 0, N%(rank)i-1
				row%(rank)i = pair%(rank)i(0, i%(rank)i);
				col%(rank)i = pair%(rank)i(1, i%(rank)i);
				%(innerLoop)s

				col%(rank)i = pair%(rank)i(0, i%(rank)i);
				row%(rank)i = pair%(rank)i(1, i%(rank)i);
				%(innerLoopConjg)s
			enddo
		""" % { "rank":self.CurRank, "innerLoop": self.GetInnerLoop(conjugate), "innerLoopConjg": self.GetInnerLoop(not conjugate)}
		return str


#-----------------------------------------------------------------------------------

class StorageGeneratorDiagonal(StorageGeneratorSimple):
	"""
	Generator for the case where the potential is diagonal in the 
	current rank, that is the potential is 0 for row != col.
	In this case, we will only loop over row == col
	"""

	def __init__(self, systemRank, curRank, innerGenerator):
		StorageGeneratorSimple.__init__(self, systemRank, curRank, innerGenerator)

	def GetLoopingCodeRecursive(self, conjugate):
		str = """
			do i%(rank)i = 0, N%(rank)i-1
				row%(rank)i = i%(rank)i
				col%(rank)i = i%(rank)i
				%(innerLoop)s
			enddo
		""" % { "rank":self.CurRank, "innerLoop": self.GetInnerLoop(conjugate)}
		return str

#-----------------------------------------------------------------------------------

class StorageGeneratorIdentity(StorageGeneratorSimple):
	"""
	Generator for the case where the potential is
	independent of this rank. In this case, we will use 
	potential index = 0 for, and loop over all row = col
	"""

	def GetLoopingCodeRecursive(self, conjugate):
		str = """
			i%(rank)i = 0
			do row%(rank)i = 0, sourceExtent%(rank)i-1
				col%(rank)i = row%(rank)i
				%(innerLoop)s
			enddo
		""" % { "rank":self.CurRank, "innerLoop": self.GetInnerLoop(conjugate)}
		return str

#-----------------------------------------------------------------------------------


class StorageGeneratorBandedBlas(StorageGeneratorBase):
	"""
	"""

	def __init__(self, systemRank, curRank, innerGenerator):
		StorageGeneratorBase.__init__(self, systemRank, curRank, innerGenerator)

	def GetParameterList(self):
		parameterList = [("pair%i" % self.CurRank, "array", 2, "integer")]
		return parameterList

	def GetParameterDeclarationCode(self):
		str = GetFortranArrayDeclaration("pair%i" % self.CurRank, 2, "integer", "in")
		str += """
			integer :: N%(rank)i, i%(rank)i, row%(rank)i, col%(rank)i, bsplineCount%(rank)i, bandCount%(rank)i
			integer :: hermRow%(rank)i, hermCol%(rank)i 
			integer :: subDiagonals%(rank)i, sourceStride%(rank)i, destStride%(rank)i
			complex (kind=dbl) :: alpha%(rank)i, beta%(rank)i
		""" % { "rank": self.CurRank }
		return str

	def GetInitializationCode(self):
		str = """
			N%(rank)i = potentialExtent%(rank)i
			bsplineCount%(rank)i = sourceExtent%(rank)i
			bandCount%(rank)i = N%(rank)i / bsplineCount%(rank)i

			subDiagonals%(rank)i = bandCount%(rank)i - 1
			sourceStride%(rank)i = 1
			destStride%(rank)i = 1
			alpha%(rank)i = scaling
			beta%(rank)i = 1.0d0
			col%(rank)i = 0
			row%(rank)i = 0
			i%(rank)i = 0

		""" % { "rank":self.CurRank }
		return str

	def GetLoopingCodeRecursive(self, conjugate):
		str = ""
		#If this is the innermost loop, we can optimize it by calling blas
		if self.InnerGenerator != None:
			str += """
				i%(rank)i = 0
				do row%(rank)i = 0, bsplineCount%(rank)i - bandCount%(rank)i - 1
					do col%(rank)i = row%(rank)i, row%(rank)i+bandCount%(rank)i - 1
						%(innerLoop)s	
						i%(rank)i = i%(rank)i + 1
					enddo
				enddo
				
				do row%(rank)i = bsplineCount%(rank)i - bandCount%(rank)i,  bsplineCount%(rank)i - 1
					do col%(rank)i = row%(rank)i, bsplineCount%(rank)i - 1
						%(innerLoop)s	
						i%(rank)i = i%(rank)i + 1
					enddo
					
					i%(rank)i = i%(rank)i + (- bsplineCount%(rank)i + row%(rank)i + bandCount%(rank)i)
				enddo
				
				i%(rank)i = 0
				do hermRow%(rank)i = 0, bsplineCount%(rank)i - bandCount%(rank)i - 1
					i%(rank)i = i%(rank)i + 1
					do hermCol%(rank)i = hermRow%(rank)i+1, hermRow%(rank)i+bandCount%(rank)i - 1
						col%(rank)i = hermRow%(rank)i
						row%(rank)i = hermCol%(rank)i
						%(innerLoopConjg)s	
						i%(rank)i = i%(rank)i + 1
					enddo
				enddo
				
				do hermRow%(rank)i = bsplineCount%(rank)i - bandCount%(rank)i,  bsplineCount%(rank)i - 1
					i%(rank)i = i%(rank)i + 1
					do hermCol%(rank)i = hermRow%(rank)i+1, bsplineCount%(rank)i - 1
						col%(rank)i = hermRow%(rank)i
						row%(rank)i = hermCol%(rank)i
						%(innerLoopConjg)s	
						i%(rank)i = i%(rank)i + 1
					enddo
					
					i%(rank)i = i%(rank)i + (- bsplineCount%(rank)i + hermRow%(rank)i + bandCount%(rank)i)
				enddo
			""" % { "rank":self.CurRank, "innerLoop": self.GetInnerLoop(conjugate), "innerLoopConjg": self.GetInnerLoop(not conjugate)}

		else:
			str += """
				i%(rank)i = 0
				call zhbmv( &
					"L", &
					bsplineCount%(rank)i, &
					subDiagonals%(rank)i, &
					alpha%(rank)i, &
					potential(%(potentialIndex)s), &
					bandCount%(rank)i, &
					source(%(colIndex)s), &
					sourceStride%(rank)i, &
					beta%(rank)i, &
					dest(%(rowIndex)s), &
					destStride%(rank)i &
				)

			""" % \
			{ \
				"rank":self.CurRank, \
				"rowIndex": self.GetIndexString("row"), \
				"potentialIndex": self.GetIndexString("i"), \
				"colIndex": self.GetIndexString("col"),  \
			}
		return str

	def GetInnerLoop(self, conjugate):
		str = ""
		if self.InnerGenerator != None:
			str += self.InnerGenerator.GetLoopingCodeRecursive(conjugate)
		else:
			#We're at the innermost loop
			conjg = ""
			if conjugate:
				conjg = "conjg"
			str += """
				dest(%(rowIndex)s) = dest(%(rowIndex)s) + %(conjg)s(potential(%(potentialIndex)s)) * scaling * source(%(colIndex)s)
			""" % \
				{ \
					"rowIndex": self.GetIndexString("row"), \
					"potentialIndex": self.GetIndexString("i"), \
					"colIndex": self.GetIndexString("col"),  \
					"conjg": conjg, \
				}
		return str


#-----------------------------------------------------------------------------------

storageGeneratorMap = { \
	"Simple": StorageGeneratorSimple, \
	"Hermitian": StorageGeneratorHermitian, \
	"Diagonal": StorageGeneratorDiagonal, \
	"Identity": StorageGeneratorIdentity, \
	"Banded": StorageGeneratorBandedBlas, \
}

class TensorMatrixMultiplyGenerator(object):

	def __init__(self, storageNameList):
		self.SystemRank = len(storageNameList)
		self.GeneratorList = self.GetGeneratorList(storageNameList)
		self.StorageNameList = storageNameList

	def GetGeneratorList(self, storageNameList):
		generatorList = []
		
		#Create recursive list of generators starting on the last rank
		systemRank = self.SystemRank
		innerGenerator = None
		for curRank in range(systemRank-1, -1, -1):
			storageName = storageNameList[curRank]
			storageClass = storageGeneratorMap[storageName]
		
			generator = storageClass(systemRank, curRank, innerGenerator)
			generatorList.insert(0, generator)
			innerGenerator = generator
		
		signature = "_".join(storageNameList)

		return generatorList
	
	def GetParameterList(self):
		generatorList = self.GeneratorList
		systemRank = self.SystemRank

		parameterList = []
		parameterList += [("potential", "array", systemRank, "complex (kind=dbl)")]
		parameterList += [("scaling", "scalar", "real (kind=dbl)")]
		parameterList += [("source", "array", systemRank, "complex (kind=dbl)")]
		parameterList += [("dest", "array", systemRank, "complex (kind=dbl)")]
		for gen in generatorList:
			parameterList += gen.GetParameterList()

		return parameterList
	
	def GetMethodName(self):
		signature = "_".join(self.StorageNameList)
		return "TensorMatrixMultiply_%s" % signature

	def GetFortranCode(self):
		systemRank = self.SystemRank 
		generatorList = self.GeneratorList

		#Get parameter list for Fortran
		parameterList = []
		for param in self.GetParameterList():
			paramName = param[0]
			paramType = param[1]
			if paramType == "array":
				paramRank = param[2]
				parameterList += GetFortranArrayParameterList(paramName, paramRank)
			else:
				parameterList += [paramName]
		parameterString = ", ".join(parameterList)
		
		str = """
		subroutine %(methodName)s(%(parameterString)s)
		
			implicit none
			include "parameters.f"
		""" % { "methodName": self.GetMethodName(), "rank": systemRank, "parameterString": parameterString}
		str += GetFortranArrayDeclaration("potential", systemRank, "complex (kind=dbl)", "in")
		str += "real (kind=dbl) :: scaling\n"
		str += GetFortranArrayDeclaration("source", systemRank, "complex (kind=dbl)", "in")
		str += GetFortranArrayDeclaration("dest", systemRank, "complex (kind=dbl)", "inout")
		str += "".join([gen.GetParameterDeclarationCode() for gen in generatorList])
		str += "".join([gen.GetInitializationCode() for gen in generatorList])
		str += generatorList[0].GetLoopingCodeRecursive(False)
		str += """
		end subroutine %(methodName)s
		""" % { "methodName": self.GetMethodName() }
		
		return PrettyPrintFortran(str)

	def GetWrapperCode(self):
		#Get c signature for fortran routine
		parameterList = []
		for param in self.GetParameterList():
			paramName = param[0]
			paramType = param[1]
			if paramType == "array":
				paramRank = param[2]
				paramDataType = param[3]
				parameterList += ["%s* %s" % (TypeMapFortranToC[paramDataType], paramName)]
				parameterList += ["int* %s" % name for name in GetFortranArrayParameterList(paramName, paramRank)[1:]]
			else:
				paramDataType = param[2]
				parameterList += ["%s* %s" % (TypeMapFortranToC[paramDataType], paramName)]
		fortranParameterString = ", ".join(parameterList)

		#Create extern function declaration
		fortranMethodName = "FORTRAN_NAME(%s)" % self.GetMethodName().lower()
		externDeclaration = """
		extern "C"
		{
			void %(fortranMethodName)s(%(fortranParameterString)s);
		}
		""" % locals()

		#Get parameter list for c++
		parameterList = []
		for param in self.GetParameterList():
			paramName = param[0]
			paramType = param[1]
			if paramType == "array":
				paramRank = param[2]
				paramDataType = param[3]
				parameterList += ["Array< %s, %i > %s" % (TypeMapFortranToC[paramDataType], paramRank, paramName)]
			else:
				paramDataType = param[2]
				parameterList += ["%s %s" % (TypeMapFortranToC[paramDataType], paramName)]
		cParameterString = ", ".join(parameterList)

		#Create temporaries for the extent of all arrays
		wrapperBody = ""
		callList = []
		for param in self.GetParameterList():
			paramName = param[0]
			paramType = param[1]
			if paramType == "array":
				paramRank = param[2]
				paramDataType = param[3]
				wrapperBody += "".join(["int %s = %s.extent(%i);\n" % (name, paramName, i) for i, name in enumerate(GetFortranArrayParameterList(paramName, paramRank)[1:])])

				callList += ["%s.data()" % (paramName)]
				callList += ["&%s" % name for name in GetFortranArrayParameterList(paramName, paramRank)[1:]]
			else:
				paramDataType = param[2]
				callList += ["&%s" % (paramName)]

		callString = ", ".join(callList)

		cMethodName = self.GetMethodName() + "_Wrapper"
		wrapperDeclaration = """
		void %(cMethodName)s(%(cParameterString)s)
		{
			%(wrapperBody)s
			%(fortranMethodName)s(%(callString)s);
		}
		""" % locals()

		return PrettyPrintC(externDeclaration + wrapperDeclaration)

	def GetBoostPythonCode(self):
		return 'def("%(name)s", %(name)s_Wrapper);\n' % {"name": self.GetMethodName()}

def GetAllPermutations(systemRank, curRank):
	if curRank == systemRank:
		yield ()
	else:
		for key in storageGeneratorMap.keys():
			for subperm in GetAllPermutations(systemRank, curRank+1):
				yield (key,) +  subperm


generatorPermutationList = []
for systemRank in range(1,4+1):
	for perm in GetAllPermutations(systemRank, 0):
		generatorPermutationList.append(TensorMatrixMultiplyGenerator(list(perm)))


def PrintFortranCode():
	for generator in generatorPermutationList:
		print generator.GetFortranCode()

def PrintWrapperCode():
	str = """
	#include <core/common.h>
	#include <core/utility/fortran.h>
	#include <boost/python.hpp>

	using namespace boost::python;
	using namespace blitz;
	"""

	for generator in generatorPermutationList :
		str += generator.GetWrapperCode()

	str += """
	void export_tensormatrixmultiply()
	{
		%(exportCode)s
	}
	""" % {"exportCode": "".join([g.GetBoostPythonCode() for g in generatorPermutationList]) }

	print PrettyPrintC(str)
