BUILT_SOURCES  = Partonic__all.cc
CLEANFILES = Partonic__all.cc

Partonic__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
QuarkoniumDecayer.h   \
HeavyDecayer.h        \
WeakPartonicDecayer.h \
BtoSGammaDecayer.h    \
PartonicDecayerBase.h 

ALL_CC_FILES = \
QuarkoniumDecayer.cc \
HeavyDecayer.cc \
WeakPartonicDecayer.cc\
BtoSGammaDecayer.cc\
PartonicDecayerBase.cc
