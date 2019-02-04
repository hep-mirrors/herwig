BUILT_SOURCES  = WeakCurrents__all.cc
CLEANFILES = WeakCurrents__all.cc

WeakCurrents__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
FourPionNovosibirskCurrent.h \
ScalarMesonCurrent.h\
ThreePionDefaultCurrent.h\
OneKaonTwoPionDefaultCurrent.h \
TwoKaonOnePionDefaultCurrent.h \
EtaPiPiDefaultCurrent.h \
ThreePionCLEOCurrent.h\
TwoPionRhoCurrent.h\
KPiKStarCurrent.h\
EtaPhotonCurrent.h\
PionPhotonCurrent.h\
TwoPionPhotonCurrent.h\
TwoPionPhotonSNDCurrent.h\
OmegaPionSNDCurrent.h\
PhiPiCurrent.h\
VectorMesonCurrent.h\
FivePionCurrent.h \
KPiCurrent.h\
OneKaonTwoPionCurrent.h\
TwoKaonOnePionCurrent.h\
TwoPionCzyzCurrent.h\
TwoKaonCzyzCurrent.h\
ThreePionCzyzCurrent.h\
FourPionCzyzCurrent.h\
EtaPiPiCurrent.h\
EtaPrimePiPiCurrent.h\
KKPiCurrent.h\
KStarKCurrent.h\
EtaPhiCurrent.h\
EtaOmegaCurrent.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
FourPionNovosibirskCurrent.cc \
ScalarMesonCurrent.cc \
ThreePionDefaultCurrent.cc \
OneKaonTwoPionDefaultCurrent.cc \
TwoKaonOnePionDefaultCurrent.cc \
EtaPiPiDefaultCurrent.cc \
ThreePionCLEOCurrent.cc \
TwoPionRhoCurrent.cc \
KPiKStarCurrent.cc \
EtaPhotonCurrent.cc \
PionPhotonCurrent.cc \
TwoPionPhotonCurrent.cc \
TwoPionPhotonSNDCurrent.cc \
OmegaPionSNDCurrent.cc \
PhiPiCurrent.cc \
VectorMesonCurrent.cc \
FivePionCurrent.cc \
KPiCurrent.cc \
OneKaonTwoPionCurrent.cc\
TwoKaonOnePionCurrent.cc\
TwoPionCzyzCurrent.cc \
TwoKaonCzyzCurrent.cc \
ThreePionCzyzCurrent.cc \
FourPionCzyzCurrent.cc\
EtaPiPiCurrent.cc\
EtaPrimePiPiCurrent.cc\
KKPiCurrent.cc\
KStarKCurrent.cc\
EtaPhiCurrent.cc\
EtaOmegaCurrent.cc
