
SUBDIRS = h2+_2d qdot_single

all clean pyclean: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 
