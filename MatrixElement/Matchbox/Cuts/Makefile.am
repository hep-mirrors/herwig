pkglib_LTLIBRARIES = HwMatchboxCuts.la

HwMatchboxCuts_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 2:0:0

nodist_HwMatchboxCuts_la_SOURCES = \
Cuts__all.cc

BUILT_SOURCES  = Cuts__all.cc
CLEANFILES = Cuts__all.cc

Cuts__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

ALL_H_FILES = \
IdentifiedParticleCut.h \
MatchboxDeltaRCut.h \
MissingPtCut.h \
FrixionePhotonSeparationCut.h \
InvariantMassCut.h \
PairPtCut.h \
PairRapidityCut.h

ALL_CC_FILES = \
IdentifiedParticleCut.cc \
MatchboxDeltaRCut.cc \
MissingPtCut.cc \
FrixionePhotonSeparationCut.cc \
InvariantMassCut.cc \
PairPtCut.cc \
PairRapidityCut.cc
