BUILT_SOURCES  = Perturbative__all.cc
CLEANFILES = Perturbative__all.cc

Perturbative__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
 SMWDecayer.h\
 SMZDecayer.h\
 SMTopDecayer.h\
 SMTopPOWHEGDecayer.h\
 SMZFermionsPOWHEGDecayer.h\
 SMWFermionsPOWHEGDecayer.h


DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
 SMWDecayer.cc  \
 SMZDecayer.cc  \
 SMTopDecayer.cc \
 SMTopPOWHEGDecayer.cc \
 SMZFermionsPOWHEGDecayer.cc  \
 SMWFermionsPOWHEGDecayer.cc  





BUILT_SOURCES  += PerturbativeHiggs__all.cc
CLEANFILES += PerturbativeHiggs__all.cc

PerturbativeHiggs__all.cc : $(DIR_H_FILES_2) $(DIR_CC_FILES_2) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES_2) > $@

EXTRA_DIST += $(ALL_H_FILES_2) $(ALL_CC_FILES_2)

DIR_H_FILES_2 = $(addprefix $(srcdir)/,$(ALL_H_FILES_2))
ALL_H_FILES_2 = \
 SMHiggsFermionsDecayer.h \
 SMHiggsFermionsPOWHEGDecayer.h \
 SMHiggsGGHiggsPPDecayer.h\
 SMHiggsWWDecayer.h

DIR_CC_FILES_2 = $(addprefix $(srcdir)/,$(ALL_CC_FILES_2))
ALL_CC_FILES_2 = \
 SMHiggsFermionsDecayer.cc \
 SMHiggsFermionsPOWHEGDecayer.cc \
 SMHiggsGGHiggsPPDecayer.cc \
 SMHiggsWWDecayer.cc 
