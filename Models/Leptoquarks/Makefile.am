BUILT_SOURCES  = Leptoquark__all.cc
CLEANFILES = Leptoquark__all.cc

Leptoquark__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
LeptoquarkModel.h LeptoquarkModel.fh \
LeptoquarkModelSLQSLQGGVertex.h \
LeptoquarkModelSLQFFVertex.h	\
LeptoquarkModelSLQSLQGVertex.h 	

ALL_CC_FILES = \
LeptoquarkModel.cc  \
LeptoquarkModelSLQSLQGGVertex.cc \
LeptoquarkModelSLQFFVertex.cc   \
LeptoquarkModelSLQSLQGVertex.cc
