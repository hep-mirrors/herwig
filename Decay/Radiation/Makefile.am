BUILT_SOURCES  = Sophty__all.cc
CLEANFILES = Sophty__all.cc

Sophty__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
SOPHTY.h \
FFDipole.fh  FFDipole.h \
IFDipole.fh  IFDipole.h \
YFSFormFactors.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
SOPHTY.cc \
FFDipole.cc \
IFDipole.cc \
YFSFormFactors.cc
