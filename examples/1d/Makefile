
SUBDIRS = pulse transformedgrid

all clean pyclean: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 
