PYPROP_PATH = ../pyprop

-include $(PYPROP_PATH)/core/makefiles/Makefile.include
-include $(PYPROP_PATH)/Makefile.platform
-include Makefile.platform

MPICXX ?= mpicxx

LDFLAGS := $(LDFLAGS) -L../pyprop/pyprop/core
CPPFLAGS := $(CPPFLAGS) -I../pyprop

CPPFLAGS := $(CPPFLAGS) 
LDFLAGS := $(LDFLAGS) -lhdf5 -lz -lm -lblitz $(LAPACK_LIBS)

benchmark: benchmark.o tensormultiply.o
	$(MPICXX) benchmark.o tensormultiply.o -o benchmark $(CPPFLAGS) $(LDFLAGS) $(FORTRAN_LIBS)

benchmark.o: benchmark.cpp h5traits.h blitzhdf.h
	$(MPICXX) -c $< -o $*.o $(CPPFLAGS)

tensormultiply.o: tensormultiply.f90
	$(FC) -c $< -o $*.o $(CPPFLAGS) $(FORTRAN_FLAGS) 


clean:
	@rm -f *.o
	@rm -f test
