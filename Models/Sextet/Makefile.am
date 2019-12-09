BUILT_SOURCES  = Sextet__all.cc
CLEANFILES = Sextet__all.cc

Sextet__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
SextetParticles.h \
SextetModel.h SextetModel.fh \
SextetGSSVertex.h  \
SextetPSSVertex.h  \
SextetGGSSVertex.h \
SextetGVVVertex.h  \
SextetPVVVertex.h  \
SextetGGVVVertex.h \
SextetFFSVertex.h \
SextetFFVVertex.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
SextetModel.cc      \
SextetGSSVertex.cc  \
SextetPSSVertex.cc  \
SextetGGSSVertex.cc \
SextetGVVVertex.cc  \
SextetPVVVertex.cc  \
SextetGGVVVertex.cc \
SextetFFSVertex.cc  \
SextetFFVVertex.cc  
