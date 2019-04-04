noinst_LTLIBRARIES = libHwMatchboxMatching.la
pkglib_LTLIBRARIES = HwQTildeMatching.la HwDipoleMatching.la

libHwMatchboxMatching_la_SOURCES = \
HardScaleProfile.h \
HardScaleProfile.cc \
ShowerApproximation.h \
ShowerApproximation.cc \
ShowerApproximationKernel.h \
ShowerApproximationKernel.cc \
ShowerApproximationGenerator.h \
ShowerApproximationGenerator.cc \
MEMatching.h \
MEMatching.cc

HwQTildeMatching_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 4:0:0
HwQTildeMatching_la_SOURCES = \
QTildeMatching.h \
QTildeMatching.cc

HwDipoleMatching_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 4:0:0
HwDipoleMatching_la_SOURCES = \
DipoleMatching.h \
DipoleMatching.cc
