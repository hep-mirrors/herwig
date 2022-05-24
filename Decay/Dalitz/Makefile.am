BUILT_SOURCES  = DalitzDecay__all.cc
CLEANFILES = DalitzDecay__all.cc

DalitzDecay__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
DalitzResonance.h\
DtoKPiPiE691.h \
FOCUSDptoKmPipPip.h \
DtoKPiPiMarkIII.h \
ScalarTo3ScalarDalitz.h\
CLEODptoKmPipPip.h\
CLEOD0toKmPipPi0.h\
CLEOD0toK0PipPim.h\
BABARDstoKpKmPip.h\
BABAREtactoKpKmEta.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
DtoKPiPiE691.cc\
FOCUSDptoKmPipPip.cc\
DtoKPiPiMarkIII.cc \
ScalarTo3ScalarDalitz.cc\
CLEODptoKmPipPip.cc\
CLEOD0toKmPipPi0.cc\
CLEOD0toK0PipPim.cc\
BABARDstoKpKmPip.cc\
BABAREtactoKpKmEta.cc
