BUILT_SOURCES  = LH__all.cc
CLEANFILES = LH__all.cc

LH__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
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

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
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
