BUILT_SOURCES  = TMDecayer__all.cc
CLEANFILES = TMDecayer__all.cc

TMDecayer__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
TensorMeson2PScalarDecayer.h\
TensorMesonVectorPScalarDecayer.h \
TensorMesonVectorVectorDecayer.h \
TensorMesonTensorScalarDecayer.h \
TensorMesonVectorScalarDecayer.h \
Spin3Meson2PScalarDecayer.h \
Spin3MesonVectorScalarDecayer.h \
Spin3MesonVectorPScalarDecayer.h \
Spin3MesonTensorPScalarDecayer.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
TensorMeson2PScalarDecayer.cc \
TensorMesonVectorPScalarDecayer.cc \
TensorMesonVectorVectorDecayer.cc \
TensorMesonTensorScalarDecayer.cc \
TensorMesonVectorScalarDecayer.cc \
Spin3Meson2PScalarDecayer.cc \
Spin3MesonVectorScalarDecayer.cc \
Spin3MesonVectorPScalarDecayer.cc \
Spin3MesonTensorPScalarDecayer.cc
