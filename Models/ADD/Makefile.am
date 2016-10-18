BUILT_SOURCES  = ADD__all.cc
CLEANFILES = ADD__all.cc

ADD__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
ADDModel.fh ADDModel.h  \
ADDModelFFGRVertex.h  \
ADDModelFFGGRVertex.h \
ADDModelFFWGRVertex.h \
ADDModelSSGRVertex.h  \
ADDModelVVGRVertex.h \
ADDModelGGGGRVertex.h \
ADDModelWWWGRVertex.h \
GravitonMassGenerator.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
ADDModel.cc  \
ADDModelFFGRVertex.cc  \
ADDModelFFGGRVertex.cc \
ADDModelFFWGRVertex.cc \
ADDModelSSGRVertex.cc  \
ADDModelVVGRVertex.cc  \
ADDModelGGGGRVertex.cc \
ADDModelWWWGRVertex.cc \
GravitonMassGenerator.cc
