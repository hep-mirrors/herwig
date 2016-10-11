BUILT_SOURCES  = VMDecayer__all.cc
DISTCLEANFILES = VMDecayer__all.cc

VMDecayer__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
a1ThreePionCLEODecayer.h \
a1SimpleDecayer.h\
a1ThreePionDecayer.h \
PVectorMesonVectorPScalarDecayer.h \
VectorMeson2FermionDecayer.h \
VectorMeson2MesonDecayer.h \
VectorMeson3PionDecayer.h \
VectorMesonPScalarFermionsDecayer.h \
VectorMesonPVectorPScalarDecayer.h \
VectorMesonVectorPScalarDecayer.h  \
VectorMesonVectorScalarDecayer.h  \
VectorMesonVectorVectorDecayer.h  \
OniumToOniumPiPiDecayer.h

ALL_CC_FILES = \
a1ThreePionCLEODecayer.cc            \
a1SimpleDecayer.cc                   \
a1ThreePionDecayer.cc                \
PVectorMesonVectorPScalarDecayer.cc  \
VectorMeson2FermionDecayer.cc        \
VectorMeson2MesonDecayer.cc          \
VectorMeson3PionDecayer.cc           \
VectorMesonPScalarFermionsDecayer.cc \
VectorMesonPVectorPScalarDecayer.cc  \
VectorMesonVectorPScalarDecayer.cc   \
VectorMesonVectorScalarDecayer.cc    \
VectorMesonVectorVectorDecayer.cc    \
OniumToOniumPiPiDecayer.cc
