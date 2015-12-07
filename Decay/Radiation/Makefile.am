noinst_LTLIBRARIES = libHwDecRad.la
pkglib_LTLIBRARIES = HwSOPHTY.la

libHwDecRad_la_SOURCES = \
DecayRadiationGenerator.cc  DecayRadiationGenerator.h  \
DecayRadiationGenerator.fh \
QEDRadiationHandler.cc  QEDRadiationHandler.fh\
QEDRadiationHandler.h

HwSOPHTY_la_SOURCES = \
SOPHTY.h SOPHTY.cc \
FFDipole.cc  FFDipole.fh  FFDipole.h \
IFDipole.cc  IFDipole.fh  IFDipole.h \
YFSFormFactors.cc  YFSFormFactors.h
HwSOPHTY_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:0:0
