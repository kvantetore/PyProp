PYPROP_ROOT  := ../..

include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := $(call ABS_PATH,$(PYPROP_ROOT)/pyprop/core)

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp \
                ../transform/bspline/bspline.cpp \
				../transform/bspline/bsplinepropagator.cpp \
				../transform/bspline/bsplinetransform.cpp \
				../representation/bspline/bsplinerepresentation.cpp \
				../representation/bspline/bsplinegridrepresentation.cpp


OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libbspline
STATICFILE   := $(PYPROP_ROOT)/core/lib/$(MODULENAME).a
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
PYPROP_LIB_PATH := $(DESTDIR)
LIBS         := $(LIBS)

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
	$(LD) $(LD_SHARED) -lcore -o $(MODULEFILE) $(OBJECTS) $(LIBS) -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH)) 

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste


$(STATICFILE): $(OBJECTS)
	rm -f $(STATICFILE)
	ar -cr $(STATICFILE) $(OBJECTS)

clean:
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf ../representation/bspline/*.o
	rm -rf ../transform/bspline/*.o
	rm -rf *.o
	rm -rf $(MODULEFILE)

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

