PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
CPPFLAGS     := $(CPPFLAGS) 

PYSTEFILES   := wrapper.pyste
SOURCEFILES  := wrapper.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)
MODULENAME   := libMixedPotential

all: $(MODULENAME).so

$(MODULENAME).so: $(OBJECTS)
	$(LD) -shared -o $(MODULENAME).so $(OBJECTS) $(LIBS) -lcore -L$(PYPROP_ROOT)/pyprop/core

wrapper.cpp: wrapper.pyste 
	$(PYSTE) $(INCLUDE) --out=wrapper.cpp --module=$(MODULENAME) wrapper.pyste

clean:
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf *.so

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

