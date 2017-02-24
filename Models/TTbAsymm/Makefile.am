BUILT_SOURCES  = TTbA__all.cc
CLEANFILES = TTbA__all.cc

TTbA__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
TTbAModel.h TTbAModel.fh \
TTbAModelWPTDVertex.h \
TTbAModelZPQQVertex.h \
TTbAModelAGQQVertex.h \
TTbAModelSU2XVertex.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
TTbAModel.cc  \
TTbAModelWPTDVertex.cc \
TTbAModelZPQQVertex.cc \
TTbAModelAGQQVertex.cc \
TTbAModelSU2XVertex.cc
