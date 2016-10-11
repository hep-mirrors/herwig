noinst_LTLIBRARIES = libHwFormFactor.la
libHwFormFactor_la_SOURCES = \
BaryonFormFactor.cc BaryonFormFactor.fh \
BaryonFormFactor.h \
ScalarFormFactor.cc ScalarFormFactor.fh \
ScalarFormFactor.h \
BtoSGammaHadronicMass.h \
BtoSGammaHadronicMass.cc BtoSGammaHadronicMass.fh




BUILT_SOURCES  = Formfactor__all.cc
DISTCLEANFILES = Formfactor__all.cc

Formfactor__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

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
