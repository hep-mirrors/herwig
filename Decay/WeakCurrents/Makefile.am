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
TwoPionPhotonCurrent.h\
VectorMesonCurrent.h\
FivePionCurrent.h \
KPiCurrent.h\
KaonThreeMesonCurrent.h\
TwoPionCzyzCurrent.h\
TwoKaonCzyzCurrent.h\
ThreePionCzyzCurrent.h\
FourPionCzyzCurrent.h\
EtaPiPiCurrent.h

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
TwoPionPhotonCurrent.cc \
VectorMesonCurrent.cc \
FivePionCurrent.cc \
KPiCurrent.cc \
KaonThreeMesonCurrent.cc \
TwoPionCzyzCurrent.cc \
TwoKaonCzyzCurrent.cc \
ThreePionCzyzCurrent.cc \
FourPionCzyzCurrent.cc\
EtaPiPiCurrent.cc
