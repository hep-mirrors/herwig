if WANT_DIPOLE
noinst_LTLIBRARIES = libHwMatchboxMatching.la
pkglib_LTLIBRARIES = HwQTildeMatching.la
endif

libHwMatchboxMatching_la_SOURCES = \
ShowerApproximation.h \
ShowerApproximation.cc \
ShowerApproximationKernel.h \
ShowerApproximationKernel.cc \
ShowerApproximationGenerator.h \
ShowerApproximationGenerator.cc \
DipoleMatching.h \
DipoleMatching.cc \
MEMatching.h \
MEMatching.cc

HwQTildeMatching_la_LDFLAGS = -module -version-info 1:0:0
HwQTildeMatching_la_SOURCES = \
QTildeMatching.h \
QTildeMatching.cc
