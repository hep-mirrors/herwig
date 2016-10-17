BUILT_SOURCES  = Sophty__all.cc
CLEANFILES = Sophty__all.cc

Sophty__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
SOPHTY.h \
FFDipole.fh  FFDipole.h \
IFDipole.fh  IFDipole.h \
YFSFormFactors.h \
DecayRadiationGenerator.h \
DecayRadiationGenerator.fh \
QEDRadiationHandler.fh \
QEDRadiationHandler.h


ALL_CC_FILES = \
SOPHTY.cc \
FFDipole.cc \
IFDipole.cc \
YFSFormFactors.cc
