PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH = $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
#CPPFLAGS     := $(CPPFLAGS) 

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libhydrogen

include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

