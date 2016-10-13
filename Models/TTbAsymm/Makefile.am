BUILT_SOURCES  = TTbA__all.cc
CLEANFILES = TTbA__all.cc

TTbA__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
TTbAModel.h TTbAModel.fh \
TTbAModelWPTDVertex.h \
TTbAModelZPQQVertex.h \
TTbAModelAGQQVertex.h \
TTbAModelSU2XVertex.h 

ALL_CC_FILES = \
TTbAModel.cc  \
TTbAModelWPTDVertex.cc \
TTbAModelZPQQVertex.cc \
TTbAModelAGQQVertex.cc \
TTbAModelSU2XVertex.cc
