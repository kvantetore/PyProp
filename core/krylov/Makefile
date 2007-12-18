include ../../Makefile.platform

SUBDIRS := 
ifdef PYPROP_USE_EXPOKIT
SUBDIRS := $(SUBDIRS) expokit
endif

ifdef PYPROP_USE_ARPACK
SUBDIRS := $(SUBDIRS) arpack
endif

ifdef PYPROP_USE_PIRAM
SUBDIRS := $(SUBDIRS) piram
endif

all:
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

clean:
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

