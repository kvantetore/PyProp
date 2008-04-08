PYPROP_ROOT  := pyprop
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp tensorpotentialmultiply_wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o) tensorpotentialmultiply.o
MODULENAME   := libpotential

-include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

tensorpotentialmultiply_wrapper.cpp: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintWrapperCode()" > $*.cpp
	
tensorpotentialmultiply.f: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintFortranCode()" > $*.f

tensorpotentialmultiply.o: tensorpotentialmultiply.f
	$(FC) $(FORTRAN_FLAGS) $< -c -o $*.o 
#	-fbounds-check 

$(MODULENAME).so: $(OBJECTS)
	$(LD) $(LD_SHARED) -o $(MODULENAME).so $(OBJECTS) $(LIBS) -lbspline -lcore -lgfortran -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH))


