#include ../../Makefile.platform

SUBDIRS = expokit

all:
ifdef PYPROP_USE_ARPACK
	cd arpack; $(MAKE) $@
endif
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

clean:
ifdef PYPROP_USE_ARPACK
	cd arpack; $(MAKE) $@
endif
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

