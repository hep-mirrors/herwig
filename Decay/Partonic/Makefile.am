BUILT_SOURCES  = Partonic__all.cc
CLEANFILES = Partonic__all.cc

Partonic__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
QuarkoniumDecayer.h   \
HeavyDecayer.h        \
WeakPartonicDecayer.h \
BtoSGammaDecayer.h    \
PartonicDecayerBase.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
QuarkoniumDecayer.cc \
HeavyDecayer.cc \
WeakPartonicDecayer.cc\
BtoSGammaDecayer.cc\
PartonicDecayerBase.cc
