PYPROP_ROOT  := ../../
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := ../../pyprop/core

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp redirect.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libredirect
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
STATICFILE   := $(PYPROP_ROOT)/core/lib/$(MODULENAME).a

#Make static exec if PYPROP_STATIC is set to 1
ifeq ($(PYPROP_STATIC),1)
STATIC_DEP = $(STATICFILE)
endif

#Make shared lib if PYPROP_DYNAMIC is set to 1
ifeq ($(PYPROP_DYNAMIC),1)
DYNAMIC_DEP = $(MODULEFILE)
endif

all: $(DYNAMIC_DEP) $(STATIC_DEP)

$(STATICFILE): $(OBJECTS)
	rm -f $(STATICFILE)
	mkdir -p $(PYPROP_ROOT)/core/lib
	$(AR) cr $(STATICFILE) $(OBJECTS)

$(MODULEFILE): $(OBJECTS) 
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) 

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

clean:
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

