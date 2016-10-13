BUILT_SOURCES  = LH__all.cc
CLEANFILES = LH__all.cc

LH__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
LHModel.h LHModel.fh \
LHFFZVertex.h  \
LHFFPVertex.h  \
LHFFGVertex.h  \
LHFFWVertex.h  \
LHWWWVertex.h  \
LHWWWWVertex.h \
LHFFHVertex.h  \
LHWWHVertex.h  \
LHWHHVertex.h  \
LHWWHHVertex.h 

ALL_CC_FILES = \
LHModel.cc  \
LHFFZVertex.cc  \
LHFFPVertex.cc  \
LHFFGVertex.cc  \
LHFFWVertex.cc  \
LHWWWVertex.cc  \
LHWWWWVertex.cc  \
LHFFHVertex.cc  \
LHWWHVertex.cc  \
LHWHHVertex.cc  \
LHWWHHVertex.cc
