include makefiles/Makefile.base
include ../Makefile.platform

SUBDIRS := redirect krylov
ifdef PYPROP_USE_ODE
SUBDIRS := $(SUBDIRS) ode
endif

#A small rule to compile only a selected few pyste files
#TESTFILES :=  python/blitzblas.pyste
#TESTFILES :=  python/representation.pyste python/combinedrepresentation.pyste
TESTFILES :=  python/combinedabsorber.pyste
BSPLINEFILES := python/bspline/bsplinerepresentation.pyste
BSPLINEREPFILES := python/bspline/bsplinetransform.pyste

test: test
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES)
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) --generate-main $(PYSTEFILES)

#A small rule to compile only vectorrepresentation
TESTFILES2 :=  python/vectorrepresentation.pyste

test2: test2
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES2)	

bspline: 
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(BSPLINEFILES)
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) --generate-main $(PYSTEFILES)

bsplinerep: 
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(BSPLINEREPFILES)
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) --generate-main $(PYSTEFILES)


