BUILT_SOURCES  = Sextet__all.cc
CLEANFILES = Sextet__all.cc

Sextet__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES) Makefile

ALL_H_FILES = \
SextetParticles.h \
SextetModel.h SextetModel.fh \
SextetGSSVertex.h  \
SextetGGSSVertex.h \
SextetGVVVertex.h  \
SextetGGVVVertex.h \
SextetFFSVertex.h \
SextetFFVVertex.h

ALL_CC_FILES = \
SextetModel.cc      \
SextetGSSVertex.cc  \
SextetGGSSVertex.cc \
SextetGVVVertex.cc  \
SextetGGVVVertex.cc \
SextetFFSVertex.cc  \
SextetFFVVertex.cc  
