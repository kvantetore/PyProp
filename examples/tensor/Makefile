PYPROP_ROOT  := pyprop
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp tensorpotentialmultiply_wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o) tensorpotentialmultiply.o
MODULENAME   := libpotential

$(MODULENAME).so: $(OBJECTS)
	$(LD) $(LD_SHARED) -o $(MODULENAME).so $(OBJECTS) $(LIBS) -lbspline -lcore -lgfortran -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH))

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

tensorpotentialmultiply_wrapper.cpp: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintWrapperCode()" > $*.cpp
	
tensorpotentialmultiply.f: tensorpotentialmultiply_generator.py
	python -c "execfile('$<'); PrintFortranCode()" > $*.f

tensorpotentialmultiply.o: tensorpotentialmultiply.f
	gfortran -g -O3 -ffree-form -ffree-line-length-none $< -c -o $*.o -fPIC
#	-fbounds-check 
clean:
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf *.so
	rm -f tensorotentialmultiply.f

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

