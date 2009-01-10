PYPROP_ROOT  := pyprop
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libpotential

-include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

$(MODULENAME).so: $(OBJECTS)
	$(LD) $(LD_SHARED) -o $(MODULENAME).so $(OBJECTS) $(FORTRAN_LIBS) $(LIBS) -lbspline -lcore -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH))


