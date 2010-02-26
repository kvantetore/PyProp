PYPROP_ROOT  := ../../..
include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/
#CPPFLAGS     := $(CPPFLAGS)
OBJECTS      := test.o

blastest: $(OBJECTS)
	$(CXX) -o blastest $(OBJECTS) -lblitz -llapack_pic -lf77blas -lcblas -latlas

clean:
	rm -f blastest
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


