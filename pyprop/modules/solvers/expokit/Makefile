PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
DESTDIR      := $(call ABS_PATH,$(PYPROP_ROOT)/pyprop/core)

PYSTEFILES   := wrapper.pyste
EXPOKIT_SOURCE  := \
			f2c/dspadm.cpp \
			f2c/dgpadm.cpp \
			f2c/dgexpv.cpp \
			f2c/dsexpv.cpp \
			f2c/zhpadm.cpp \
			f2c/zgpadm.cpp \
			f2c/zgexpv.cpp \
			f2c/zhexpv.cpp \
			expokitpropagator.cpp \
			wrapper.cpp
SOURCEFILES  := $(EXPOKIT_SOURCE)

OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libexpokit
MODULEFILE   := $(DESTDIR)/$(MODULENAME).so
PYPROP_LIB_PATH := $(DESTDIR)

all: $(MODULEFILE)

$(MODULEFILE): $(OBJECTS)
	rm -f $(MODULEFILE)
	$(LD) $(LD_SHARED) -o $(MODULEFILE) $(OBJECTS) $(LIBS) -lcore -L$(PYPROP_LIB_PATH) $(call STATIC_LINK_DIR,$(PYPROP_LIB_PATH))

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

clean:
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf f2c/*.o
	rm -rf $(MODULEFILE)

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

