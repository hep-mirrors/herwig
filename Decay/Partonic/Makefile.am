BUILT_SOURCES  = Partonic__all.cc
DISTCLEANFILES = Partonic__all.cc

Partonic__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

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
