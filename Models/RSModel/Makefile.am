BUILT_SOURCES  = RS__all.cc
CLEANFILES = RS__all.cc

RS__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
RSModel.fh RSModel.h \
RSModelFFGRVertex.h  \
RSModelFFGGRVertex.h \
RSModelFFWGRVertex.h \
RSModelSSGRVertex.h  \
RSModelVVGRVertex.h  \
RSModelWWWGRVertex.h \
RSModelGGGGRVertex.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
RSModel.cc  \
RSModelFFGRVertex.cc  \
RSModelFFGGRVertex.cc \
RSModelFFWGRVertex.cc \
RSModelSSGRVertex.cc  \
RSModelVVGRVertex.cc  \
RSModelWWWGRVertex.cc \
RSModelGGGGRVertex.cc 
