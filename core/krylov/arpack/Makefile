PYPROP_ROOT  := ../../..

ARPACK_LIB   = $(ARPACK_DIR)/libarpack.a
#PARPACK_LIB   = $(ARPACK_DIR)/libparpack.a

include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := $(call ABS_PATH,$(PYPROP_ROOT)/pyprop/core)

PYSTEFILES   := wrapper.pyste
ARPACK_SOURCE  := \
			arpackpropagator.cpp \
			wrapper.cpp
SOURCEFILES  := $(ARPACK_SOURCE)

OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libarpack
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
PYPROP_LIB_PATH := $(DESTDIR)
LIBS         := -lcore $(PARPACK_LIB) $(ARPACK_LIB) $(FORTRAN_LIBS)
all: $(MODULEFILE)

arpacktest: arpacktest.cpp
	mpicxx  -I/opt/blitz/include -DPYPROP_USE_BLAS   -DBZ_DEBUG -g -I. -I/usr/lib/python2.5/site-packages/numpy/core/include -I/usr/include/python2.5  -I../../../ -ftemplate-depth-255 -DBOOST_PYTHON_DYNAMIC_LIB -fPIC -D__GCCXML__  -MD -MF arpacktest.d -c arpacktest.cpp -o arpacktest.o
	mpicxx  -DBZ_DEBUG -g  -L/opt/blitz/lib -fPIC  -o arpacktest arpacktest.o  /opt/ARPACK/libparpack.a /opt/ARPACK/libarpack.a -lmpi_f77 -lgfortran  -L/home/tore/prog/pyprop/pyprop/core -Wl,-rpath=/home/tore/prog/pyprop/pyprop/core -llapack -lblas -lblitz

$(MODULEFILE): $(OBJECTS)
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH)) 

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

clean:
	rm -f arpacktest
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf $(MODULEFILE)

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

