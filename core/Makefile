include makefiles/Makefile.base

SUBDIRS := redirect krylov
ifeq ($(PYPROP_USE_ODE),1)
SUBDIRS := $(SUBDIRS) ode
endif
ifeq ($(PYPROP_USE_BSPLINE),1)
SUBDIRS := $(SUBDIRS) bspline
endif
ifeq ($(PYPROP_USE_RUNGEKUTTA),1)
SUBDIRS := $(SUBDIRS) runge-kutta
endif

#A small rule to compile only a selected few pyste files
#TESTFILES :=  python/blitzblas.pyste
#TESTFILES :=  python/tensorpotential.pyste
TESTFILES :=  python/overlapmatrix.pyste

#TESTFILES :=  python/distributedmodel.pyste
BSPLINEREPFILES := python/bspline/bsplinerepresentation.pyste
BSPLINEFILES := python/bspline/bsplinetransform.pyste

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
