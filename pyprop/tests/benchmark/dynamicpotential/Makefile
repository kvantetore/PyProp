PYPROP_HOME  := ../../../..

include $(PYPROP_HOME)/core/makefiles/Makefile.include
include $(PYPROP_HOME)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_HOME)/
CPPFLAGS     := $(CPPFLAGS) 

SOURCEFILES  := dynpot.cpp
OBJECTS      := $(SOURCEFILES:.cpp=.o)

dynpot: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o dynpot $(OBJECTS) $(LIBS) -lpython2.4 -L$(PYPROP_HOME)/pyprop/core

clean:
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o
	rm -rf dynpot

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

