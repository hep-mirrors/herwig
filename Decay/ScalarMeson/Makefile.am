BUILT_SOURCES  = SMDecayer__all.cc
CLEANFILES = SMDecayer__all.cc

SMDecayer__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
EtaPiGammaGammaDecayer.h\
EtaPiPiGammaDecayer.h \
EtaPiPiFermionsDecayer.h \
EtaPiPiPiDecayer.h \
PScalar4FermionsDecayer.h\
PScalarLeptonNeutrinoDecayer.h\
PScalarPScalarVectorDecayer.h  \
PScalarVectorFermionsDecayer.h\
PScalarVectorVectorDecayer.h\
ScalarMesonTensorScalarDecayer.h\
ScalarScalarScalarDecayer.h  \
SemiLeptonicScalarDecayer.h  \
ScalarMesonFactorizedDecayer.h \
ScalarVectorVectorDecayer.h \
Scalar2FermionsDecayer.h \
PseudoScalar2FermionsDecayer.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
EtaPiGammaGammaDecayer.cc \
EtaPiPiGammaDecayer.cc \
EtaPiPiFermionsDecayer.cc \
EtaPiPiPiDecayer.cc \
PScalar4FermionsDecayer.cc \
PScalarLeptonNeutrinoDecayer.cc \
PScalarPScalarVectorDecayer.cc \
PScalarVectorFermionsDecayer.cc  \
PScalarVectorVectorDecayer.cc \
ScalarMesonTensorScalarDecayer.cc \
ScalarScalarScalarDecayer.cc \
SemiLeptonicScalarDecayer.cc \
ScalarMesonFactorizedDecayer.cc \
ScalarVectorVectorDecayer.cc \
Scalar2FermionsDecayer.cc  \
PseudoScalar2FermionsDecayer.cc 
