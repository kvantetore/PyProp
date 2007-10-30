PYPROP_ROOT  := ../..

ODE_LIB   = ode.o

include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := $(call ABS_PATH,$(PYPROP_ROOT)/pyprop/core)

PYSTEFILES   := wrapper.pyste
ODE_SOURCE  := \
			odewrapper.cpp \
			wrapper.cpp
SOURCEFILES  := $(ODE_SOURCE)

OBJECTS      := $(SOURCEFILES:.cpp=.o) $(ODE_LIB) 
MODULENAME   := libode
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
PYPROP_LIB_PATH := $(DESTDIR)
LIBS         := -lcore $(FORTRAN_LIBS)
all: $(MODULEFILE)

$(MODULEFILE): $(OBJECTS)
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH)) 

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

$(ODE_LIB): ode.f
	$(FC) $(FORTRAN_FLAGS) -c ode.f 

clean:
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

