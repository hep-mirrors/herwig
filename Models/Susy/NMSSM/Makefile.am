BUILT_SOURCES  = NMSSM__all.cc
CLEANFILES = NMSSM__all.cc

NMSSM__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

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
