noinst_LTLIBRARIES = libHwDecRad.la

libHwDecRad_la_SOURCES = \
DecayRadiationGenerator.cc  DecayRadiationGenerator.h  \
DecayRadiationGenerator.fh \
QEDRadiationHandler.cc  QEDRadiationHandler.fh\
QEDRadiationHandler.h




BUILT_SOURCES  = Sophty__all.cc
CLEANFILES = Sophty__all.cc

Sophty__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
SOPHTY.h \
FFDipole.fh  FFDipole.h \
IFDipole.fh  IFDipole.h \
YFSFormFactors.h

ALL_CC_FILES = \
SOPHTY.cc \
FFDipole.cc \
IFDipole.cc \
YFSFormFactors.cc
