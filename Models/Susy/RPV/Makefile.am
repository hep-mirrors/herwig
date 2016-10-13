BUILT_SOURCES  = RPV__all.cc
CLEANFILES = RPV__all.cc

RPV__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES) Makefile

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
