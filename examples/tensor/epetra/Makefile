PYPROP_ROOT  := pyprop
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

PYPROP_LIB_PATH := $(PYPROP_ROOT)/pyprop/core
INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libepetratest

-include $(PYPROP_ROOT)/core/makefiles/Makefile.extension

CPPFLAGS := $(CPPFLAGS)
LAPACK_LIBS  := $(LAPACK_LIBS)
