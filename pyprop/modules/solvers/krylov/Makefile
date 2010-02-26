include ../../Makefile.platform

SUBDIRS := 
ifeq ($(PYPROP_USE_EXPOKIT),1)
SUBDIRS := $(SUBDIRS) expokit
endif

ifeq ($(PYPROP_USE_ARPACK),1)
SUBDIRS := $(SUBDIRS) arpack
endif

ifeq ($(PYPROP_USE_PIRAM),1)
SUBDIRS := $(SUBDIRS) piram
endif

ifeq ($(PYPROP_USE_PAMP),1)
SUBDIRS := $(SUBDIRS) pamp
endif

ifeq ($(PYPROP_USE_GMRES),1)
SUBDIRS := $(SUBDIRS) gmres
endif

all:
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

clean:
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 

