PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
CPPFLAGS     := $(CPPFLAGS) 

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libMixedPotential

include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

