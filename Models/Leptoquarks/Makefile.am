BUILT_SOURCES  = Leptoquark__all.cc
CLEANFILES = Leptoquark__all.cc

Leptoquark__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
LeptoquarkModel.h LeptoquarkModel.fh \
LeptoquarkModelSLQSLQGGVertex.h \
LeptoquarkModelSLQFFVertex.h	\
LeptoquarkModelSLQSLQGVertex.h 	

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
LeptoquarkModel.cc  \
LeptoquarkModelSLQSLQGGVertex.cc \
LeptoquarkModelSLQFFVertex.cc   \
LeptoquarkModelSLQSLQGVertex.cc
