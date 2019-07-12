noinst_LTLIBRARIES = libHwMatchboxPhasespace.la

nodist_libHwMatchboxPhasespace_la_SOURCES = \
Phasespace__all.cc

BUILT_SOURCES  = Phasespace__all.cc
CLEANFILES = Phasespace__all.cc

Phasespace__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
FFLightInvertedTildeKinematics.h \
FFLightTildeKinematics.h \
FFMassiveInvertedTildeKinematics.h \
FFMassiveTildeKinematics.h \
FILightInvertedTildeKinematics.h \
FILightTildeKinematics.h \
FIMassiveInvertedTildeKinematics.h \
FIMassiveTildeKinematics.h \
IFLightInvertedTildeKinematics.h \
IFLightTildeKinematics.h \
IFMassiveInvertedTildeKinematics.h \
IFMassiveTildeKinematics.h \
IILightInvertedTildeKinematics.h \
IILightTildeKinematics.h \
InvertedTildeKinematics.fh \
InvertedTildeKinematics.h \
MatchboxPhasespace.h \
MatchboxRambo.h \
PhasespaceHelpers.h \
RandomHelpers.h \
TildeKinematics.fh \
TildeKinematics.h \
TreePhasespace.h \
TreePhasespaceChannels.h \
MatchboxReference.h \
FlatInvertiblePhasespace.h \
FlatInvertibleLabframePhasespace.h \
PhasespaceCouplings.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
FFLightInvertedTildeKinematics.cc \
FFLightTildeKinematics.cc \
FFMassiveInvertedTildeKinematics.cc \
FFMassiveTildeKinematics.cc \
FILightInvertedTildeKinematics.cc \
FILightTildeKinematics.cc \
FIMassiveInvertedTildeKinematics.cc \
FIMassiveTildeKinematics.cc \
IFLightInvertedTildeKinematics.cc \
IFLightTildeKinematics.cc \
IFMassiveInvertedTildeKinematics.cc \
IFMassiveTildeKinematics.cc \
IILightInvertedTildeKinematics.cc \
IILightTildeKinematics.cc \
InvertedTildeKinematics.cc \
MatchboxPhasespace.cc \
MatchboxRambo.cc \
PhasespaceHelpers.cc \
TildeKinematics.cc \
TreePhasespace.cc \
TreePhasespaceChannels.cc \
MatchboxReference.cc \
FlatInvertiblePhasespace.cc \
FlatInvertibleLabframePhasespace.cc \
PhasespaceCouplings.cc
