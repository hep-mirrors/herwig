BUILT_SOURCES  = TMDecayer__all.cc
DISTCLEANFILES = TMDecayer__all.cc

TMDecayer__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
TensorMeson2PScalarDecayer.h\
TensorMesonVectorPScalarDecayer.h \
TensorMesonVectorVectorDecayer.h

ALL_CC_FILES = \
TensorMeson2PScalarDecayer.cc \
TensorMesonVectorPScalarDecayer.cc \
TensorMesonVectorVectorDecayer.cc 
