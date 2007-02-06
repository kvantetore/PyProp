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

all: $(MODULEFILE)

$(MODULEFILE): $(OBJECTS)
	$(LD) -shared -o $(MODULEFILE) $(OBJECTS) $(LIBS) 

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

#C++ Compile rule
%.o : %.cpp
	$(CXX) -MD -c $< -o $*.o
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

-include $(DEPENDENCIES)

