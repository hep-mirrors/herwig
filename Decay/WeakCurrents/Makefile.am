pkglib_LTLIBRARIES =
noinst_LTLIBRARIES = libHwWeakCurrent.la 

if WANT_NEWDECAYERS
pkglib_LTLIBRARIES += HwWeakCurrents.la 
endif

libHwWeakCurrent_la_SOURCES = \
WeakDecayCurrent.cc WeakDecayCurrent.fh WeakDecayCurrent.h WeakDecayCurrent.icc

HwWeakCurrents_la_LDFLAGS = -module -version-info 2:0:0
HwWeakCurrents_la_SOURCES = \
FourPionNovosibirskCurrent.cc FourPionNovosibirskCurrent.fh \
FourPionNovosibirskCurrent.h FourPionNovosibirskCurrent.icc \
ScalarMesonCurrent.cc ScalarMesonCurrent.fh ScalarMesonCurrent.h ScalarMesonCurrent.icc \
ThreeMesonCurrentBase.cc ThreeMesonCurrentBase.fh ThreeMesonCurrentBase.h ThreeMesonCurrentBase.icc\
ThreeMesonDefaultCurrent.cc ThreeMesonDefaultCurrent.fh \
ThreeMesonDefaultCurrent.h ThreeMesonDefaultCurrent.icc \
ThreePionCLEOCurrent.cc ThreePionCLEOCurrent.fh ThreePionCLEOCurrent.h ThreePionCLEOCurrent.icc \
TwoMesonRhoKStarCurrent.cc TwoMesonRhoKStarCurrent.fh \
TwoMesonRhoKStarCurrent.h TwoMesonRhoKStarCurrent.icc \
TwoPionPhotonCurrent.cc TwoPionPhotonCurrent.fh TwoPionPhotonCurrent.h TwoPionPhotonCurrent.icc\
VectorMesonCurrent.cc VectorMesonCurrent.fh VectorMesonCurrent.h VectorMesonCurrent.icc\
LeptonNeutrinoCurrent.cc LeptonNeutrinoCurrent.fh LeptonNeutrinoCurrent.h LeptonNeutrinoCurrent.icc \
FivePionCurrent.cc FivePionCurrent.fh FivePionCurrent.h FivePionCurrent.icc

