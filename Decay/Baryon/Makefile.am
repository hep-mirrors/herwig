BUILT_SOURCES  = BaryonDecayer__all.cc
CLEANFILES = BaryonDecayer__all.cc

BaryonDecayer__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
Baryon1MesonDecayerBase.h   \
Baryon1MesonDecayerBase.fh  \
BaryonFactorizedDecayer.h   \
KornerKramerCharmDecayer.h  \
OmegaXiStarPionDecayer.h    \
SemiLeptonicBaryonDecayer.h \
StrongHeavyBaryonDecayer.h  \
NonLeptonicHyperonDecayer.h \
NonLeptonicOmegaDecayer.h   \
RadiativeHyperonDecayer.h   \
RadiativeHeavyBaryonDecayer.h \
SU3BaryonDecupletOctetPhotonDecayer.h\
SU3BaryonDecupletOctetScalarDecayer.h\
SU3BaryonOctetDecupletScalarDecayer.h\
SU3BaryonOctetOctetPhotonDecayer.h   \
SU3BaryonOctetOctetScalarDecayer.h   \
SU3BaryonSingletOctetPhotonDecayer.h \
SU3BaryonSingletOctetScalarDecayer.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
Baryon1MesonDecayerBase.cc   \
BaryonFactorizedDecayer.cc   \
KornerKramerCharmDecayer.cc  \
OmegaXiStarPionDecayer.cc    \
SemiLeptonicBaryonDecayer.cc \
StrongHeavyBaryonDecayer.cc  \
NonLeptonicHyperonDecayer.cc           \
NonLeptonicOmegaDecayer.cc             \
RadiativeHyperonDecayer.cc             \
RadiativeHeavyBaryonDecayer.cc         \
SU3BaryonDecupletOctetPhotonDecayer.cc \
SU3BaryonDecupletOctetScalarDecayer.cc \
SU3BaryonOctetDecupletScalarDecayer.cc \
SU3BaryonOctetOctetPhotonDecayer.cc    \
SU3BaryonOctetOctetScalarDecayer.cc    \
SU3BaryonSingletOctetPhotonDecayer.cc  \
SU3BaryonSingletOctetScalarDecayer.cc
