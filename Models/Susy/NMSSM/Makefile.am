BUILT_SOURCES  = NMSSM__all.cc
CLEANFILES = NMSSM__all.cc

NMSSM__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
NMSSM.h NMSSM.fh \
NMSSMFFHVertex.h   \
NMSSMWWHVertex.h   \
NMSSMWHHVertex.h   \
NMSSMHSFSFVertex.h \
NMSSMGOGOHVertex.h \
NMSSMHHHVertex.h   \
NMSSMGGHVertex.h   \
NMSSMPPHVertex.h   \
NMSSMWWHHVertex.h  

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
NMSSM.cc  \
NMSSMFFHVertex.cc   \
NMSSMWWHVertex.cc   \
NMSSMWHHVertex.cc   \
NMSSMHSFSFVertex.cc\
NMSSMGOGOHVertex.cc \
NMSSMHHHVertex.cc   \
NMSSMGGHVertex.cc \
NMSSMPPHVertex.cc  \
NMSSMWWHHVertex.cc
