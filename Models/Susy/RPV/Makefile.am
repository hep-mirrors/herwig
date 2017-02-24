BUILT_SOURCES  = RPV__all.cc
CLEANFILES = RPV__all.cc

RPV__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
RPV.h RPV.fh \
RPVhelper.h \
RPVLLEVertex.h \
RPVLQDVertex.h \
RPVUDDVertex.h \
RPVFFZVertex.h \
RPVFFWVertex.h \
RPVWSSVertex.h \
RPVFFSVertex.h \
RPVWWHVertex.h \
RPVSSSVertex.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
RPV.cc  \
RPVLLEVertex.cc\
RPVLQDVertex.cc\
RPVUDDVertex.cc\
RPVFFZVertex.cc\
RPVFFWVertex.cc\
RPVWSSVertex.cc\
RPVFFSVertex.cc\
RPVWWHVertex.cc\
RPVSSSVertex.cc
