
SUBDIRS = h2+ kulander

all clean pyclean: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) $@); done 
