BUILT_SOURCES  = VMDecayer__all.cc
CLEANFILES = VMDecayer__all.cc

VMDecayer__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
a1ThreePionCLEODecayer.h \
a1SimpleDecayer.h\
a1ThreePionDecayer.h \
VectorMeson2FermionDecayer.h \
VectorMeson2SpinHalfBaryonsDecayer.h \
VectorMeson2SpinThreeHalfBaryonsDecayer.h \
VectorMeson2MesonDecayer.h \
VectorMesonPScalarFermionsDecayer.h \
VectorMesonVectorPScalarDecayer.h  \
VectorMesonVectorScalarDecayer.h  \
VectorMesonVectorVectorDecayer.h  \
VectorMesonTensorVectorDecayer.h \
PseudoVectorMesonVectorVectorDecayer.h  \
OniumToOniumPiPiDecayer.h \
f1RhoPiPiDecayer.h \
f1FourPiDecayer.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
a1ThreePionCLEODecayer.cc            \
a1SimpleDecayer.cc                   \
a1ThreePionDecayer.cc                \
VectorMeson2FermionDecayer.cc        \
VectorMeson2SpinHalfBaryonsDecayer.cc \
VectorMeson2SpinThreeHalfBaryonsDecayer.cc \
VectorMeson2MesonDecayer.cc          \
VectorMesonPScalarFermionsDecayer.cc \
VectorMesonVectorPScalarDecayer.cc   \
VectorMesonVectorScalarDecayer.cc    \
VectorMesonVectorVectorDecayer.cc    \
VectorMesonTensorVectorDecayer.cc    \
PseudoVectorMesonVectorVectorDecayer.cc  \
OniumToOniumPiPiDecayer.cc           \
f1RhoPiPiDecayer.cc                  \
f1FourPiDecayer.cc
