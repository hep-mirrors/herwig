BUILT_SOURCES  = Leptoquark__all.cc
DISTCLEANFILES = Leptoquark__all.cc

Leptoquark__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

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
