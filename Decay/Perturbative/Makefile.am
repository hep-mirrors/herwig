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
 SMHiggsFermionsDecayer.h \
 SMHiggsGGHiggsPPDecayer.h\
 SMHiggsWWDecayer.h


DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
 SMWDecayer.cc  \
 SMZDecayer.cc  \
 SMTopDecayer.cc \
 SMHiggsFermionsDecayer.cc \
 SMHiggsGGHiggsPPDecayer.cc \
 SMHiggsWWDecayer.cc 
