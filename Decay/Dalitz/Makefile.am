BUILT_SOURCES  = DalitzDecay__all.cc
CLEANFILES = DalitzDecay__all.cc

DalitzDecay__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
DalitzResonance.h\
FOCUSDptoKmPipPip.h \
ScalarTo3ScalarDalitz.h\
BABAREtactoKpKmEta.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
FOCUSDptoKmPipPip.cc\
ScalarTo3ScalarDalitz.cc\
BABAREtactoKpKmEta.cc
