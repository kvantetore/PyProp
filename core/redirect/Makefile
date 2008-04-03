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

all: $(MODULEFILE) $(STATICFILE)

$(STATICFILE): $(OBJECTS)
	rm -f $(STATICFILE)
	mkdir -p $(PYPROP_ROOT)/core/lib
	$(AR) cr $(STATICFILE) $(OBJECTS)

$(MODULEFILE): $(OBJECTS) 
	rm -f $(MODULEFILE)
ifdef PYPROP_SHARED	
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) 
endif	

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

