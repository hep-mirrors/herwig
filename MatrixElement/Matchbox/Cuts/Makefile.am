pkglib_LTLIBRARIES = HwMatchboxCuts.la

HwMatchboxCuts_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 4:0:0

nodist_HwMatchboxCuts_la_SOURCES = \
Cuts__all.cc

BUILT_SOURCES  = Cuts__all.cc
CLEANFILES = Cuts__all.cc

Cuts__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
IdentifiedParticleCut.h \
MatchboxDeltaRCut.h \
MissingPtCut.h \
FrixionePhotonSeparationCut.h \
InvariantMassCut.h \
PairPtCut.h \
PairRapidityCut.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
IdentifiedParticleCut.cc \
MatchboxDeltaRCut.cc \
MissingPtCut.cc \
FrixionePhotonSeparationCut.cc \
InvariantMassCut.cc \
PairPtCut.cc \
PairRapidityCut.cc
