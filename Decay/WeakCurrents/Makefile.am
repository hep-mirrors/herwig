BUILT_SOURCES  = WeakCurrents__all.cc
CLEANFILES = WeakCurrents__all.cc

WeakCurrents__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
FourPionNovosibirskCurrent.h \
ScalarMesonCurrent.h\
ThreeMesonCurrentBase.h \
ThreeMesonDefaultCurrent.h\
ThreePionCLEOCurrent.h\
TwoMesonRhoKStarCurrent.h\
TwoPionPhotonCurrent.h\
VectorMesonCurrent.h\
FivePionCurrent.h \
KPiCurrent.h\
KaonThreeMesonCurrent.h \
WeakDecayCurrent.fh WeakDecayCurrent.h \
LeptonNeutrinoCurrent.fh \
LeptonNeutrinoCurrent.h


ALL_CC_FILES = \
FourPionNovosibirskCurrent.cc \
ScalarMesonCurrent.cc \
ThreeMesonCurrentBase.cc \
ThreeMesonDefaultCurrent.cc \
ThreePionCLEOCurrent.cc \
TwoMesonRhoKStarCurrent.cc \
TwoPionPhotonCurrent.cc \
VectorMesonCurrent.cc \
FivePionCurrent.cc \
KPiCurrent.cc \
KaonThreeMesonCurrent.cc 
