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

SUPERLU_PATH = /opt/superlu_serial
SUPERLU_LIBS = -L$(SUPERLU_PATH)/lib -lsuperlu_serial
SUPERLU_INCLUDE = -I$(SUPERLU_PATH)/SRC
CPPFLAGS := $(CPPFLAGS) $(SUPERLU_INCLUDE)
LAPACK_LIBS  := $(LAPACK_LIBS) $(SUPERLU_LIBS)


