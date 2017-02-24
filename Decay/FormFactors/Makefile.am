BUILT_SOURCES  = Formfactor__all.cc
CLEANFILES = Formfactor__all.cc

Formfactor__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
BallZwickyScalarFormFactor.h \
BallZwickyVectorFormFactor.h \
BaryonSimpleFormFactor.h \
BaryonThreeQuarkModelFormFactor.h \
ChengHeavyBaryonFormFactor.h \
ISGW2FormFactor.h \
ISGWFormFactor.h \
LambdabExcitedLambdacSumRuleFormFactor.h \
LightBaryonQuarkModelFormFactor.h \
SingletonFormFactor.h \
WSBFormFactor.h \
KiselevBcFormFactor.h \
MelikhovFormFactor.h \
MelikhovStechFormFactor.h \
BtoSGammaFlatEnergy.h \
BtoSGammaKagan.h \
HQETFormFactor.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
BallZwickyScalarFormFactor.cc \
BallZwickyVectorFormFactor.cc \
BaryonSimpleFormFactor.cc \
BaryonThreeQuarkModelFormFactor.cc \
ChengHeavyBaryonFormFactor.cc \
ISGW2FormFactor.cc \
ISGWFormFactor.cc \
LambdabExcitedLambdacSumRuleFormFactor.cc \
LightBaryonQuarkModelFormFactor.cc \
SingletonFormFactor.cc \
WSBFormFactor.cc \
KiselevBcFormFactor.cc \
MelikhovFormFactor.cc \
MelikhovStechFormFactor.cc \
BtoSGammaFlatEnergy.cc \
BtoSGammaKagan.cc \
HQETFormFactor.cc
