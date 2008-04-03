PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libh2p

include $(PYPROP_ROOT)/core/makefiles/Makefile.extension
