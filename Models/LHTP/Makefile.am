BUILT_SOURCES  = LHTP__all.cc
DISTCLEANFILES = LHTP__all.cc

LHTP__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
LHTPModel.fh     LHTPModel.h \
LHTPWWHVertex.h \
LHTPFFGVertex.h \
LHTPFFPVertex.h \
LHTPFFWVertex.h \
LHTPFFZVertex.h \
LHTPWWWVertex.h \
LHTPHHHVertex.h \
LHTPWHHVertex.h \
LHTPFFHVertex.h

ALL_CC_FILES = \
LHTPModel.cc      \
LHTPWWHVertex.cc  \
LHTPFFGVertex.cc  \
LHTPFFPVertex.cc  \
LHTPFFWVertex.cc  \
LHTPFFZVertex.cc  \
LHTPWWWVertex.cc  \
LHTPHHHVertex.cc  \
LHTPWHHVertex.cc  \
LHTPFFHVertex.cc  
