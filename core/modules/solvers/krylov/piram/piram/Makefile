PYPROP_ROOT  := ../../../..

include $(PYPROP_ROOT)/core/makefiles/Makefile.include
include $(PYPROP_ROOT)/Makefile.platform

INCLUDE      := $(INCLUDE) -I$(PYPROP_ROOT)/

all: piramtest

qrtest: qrtest.o
	$(MPICXX) $(LDFLAGS) -o qrtest qrtest.o -lblitz -lblas -llapack

piramtest: piramtest.o
	$(MPICXX) $(LDFLAGS) -o piramtest piramtest.o $(LAPACK_LIBS)

test: test.o
	$(MPICXX) $(LDFLAGS) -o test test.o -lblitz -lblas -llapack

clean:
	rm -f test
	rm -f piramtest
	rm -rf .deps
	mkdir .deps
	rm -rf wrapper.cpp
	rm -rf *.o

#autodependencies
DEPDIR        = .deps
df            = $(DEPDIR)/$(*F)
DEPENDENCIES  = $(addprefix $(DEPDIR)/, $(SOURCEFILES:%.cpp=%.P))

-include $(DEPENDENCIES)

