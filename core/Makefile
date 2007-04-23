include makefiles/Makefile.base
include ../Makefile.platform

SUBDIRS = redirect

#A small rule to compile only a selected few pyste files
TESTFILES := python/distributedmodel.pyste


test: test
	$(PYSTE) --multiple --out=$(PYSTEOUTPUTDIR) --module=$(MODULE_NAME) $(TESTFILES)
