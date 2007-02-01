include makefiles/Makefile.base
include ../Makefile.platform

#A small rule to compile only a selected few pyste files
TESTFILES := python/representation.pyste python/sphericalrepresentation.pyste 
test:
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES)
