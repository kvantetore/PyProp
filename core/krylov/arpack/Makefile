PYPROP_ROOT  := ../../..

ARPACK_LIB   = $(ARPACK_DIR)/libarpack.a

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

all: $(MODULEFILE)

$(MODULEFILE): $(OBJECTS)
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) -lcore $(ARPACK_LIB) -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH)) 

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

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

