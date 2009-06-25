PYPROP_ROOT  := ../../..

include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := $(call ABS_PATH,$(PYPROP_ROOT)/pyprop/core)

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := \
			piramsolver.cpp \
			wrapper.cpp

OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libpiram
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
STATICFILE   := $(PYPROP_ROOT)/core/lib/$(MODULENAME).a
PYPROP_LIB_PATH := $(DESTDIR)
LIBS         := $(LIBS) -lcore 

#Make static exec if PYPROP_STATIC is set to 1
ifeq ($(PYPROP_STATIC),1)
STATIC_DEP = $(STATICFILE)
endif

#Make shared lib if PYPROP_DYNAMIC is set to 1
ifeq ($(PYPROP_DYNAMIC),1)
DYNAMIC_DEP = $(MODULEFILE)
endif

all: $(DYNAMIC_DEP) $(STATIC_DEP)

$(MODULEFILE): $(OBJECTS)
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH)) 

$(STATICFILE): $(OBJECTS)
	rm -f $(STATICFILE)
	mkdir -p $(PYPROP_ROOT)/core/lib
	$(AR) cr $(STATICFILE) $(OBJECTS)

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

clean:
	rm -f arpacktest
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf $(MODULEFILE)
	rm -rf $(STATICFILE)

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

