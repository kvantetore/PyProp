PYPROP_ROOT  := pyprop
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp tensorpotentialmultiply_wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o) tensorpotentialmultiply.o parallel_multiply.o test.o
MODULENAME   := libpotential

-include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

tensorpotentialmultiply_wrapper.cpp: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintWrapperCode()" > $*.cpp
	python -c "execfile('$<'); PrintWrapperHeader()" > $*.h

tensorpotentialmultiply.f: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintFortranCode()" > $*.f

tensorpotentialmultiply.o: tensorpotentialmultiply.f
	$(FC) $(FORTRAN_FLAGS) $< -c -o $*.o 
#	-fbounds-check 
#
parallel_multiply.o: parallel_multiply.f
	mpif90 $(FORTRAN_FLAGS) $< -c -o $*.o 

test.o: test.f
	mpif90 $(FORTRAN_FLAGS) $< -c -o $*.o 


$(MODULENAME).so: $(OBJECTS)
	$(LD) $(LD_SHARED) -o $(MODULENAME).so $(OBJECTS) -lmpi_f90 -lmpi_f77  $(LIBS) -lbspline -lcore -lgfortran -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH))


