include makefiles/Makefile.base
include ../Makefile.platform

SUBDIRS := redirect krylov
ifdef PYPROP_USE_ODE
SUBDIRS := $(SUBDIRS) ode
endif

#A small rule to compile only a selected few pyste files
TESTFILES :=  python/distributedmodel.pyste

test: test
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES)
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) --generate-main $(PYSTEFILES)

#A small rule to compile only vectorrepresentation
TESTFILES2 :=  python/vectorrepresentation.pyste

test2: test2
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES2)	
