pkglib_LTLIBRARIES = HwMatchboxBuiltin.la

HwMatchboxBuiltin_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 5:0:0

nodist_HwMatchboxBuiltin_la_SOURCES = \
Amplitudes__all.cc

BUILT_SOURCES  = Amplitudes__all.cc
CLEANFILES = Amplitudes__all.cc

Amplitudes__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
MatchboxCurrents.h \
MatchboxZGammaAmplitude.h \
MatchboxAmplitudellbarqqbargg.h \
MatchboxAmplitudellbarqqbarg.h \
MatchboxAmplitudellbarqqbar.h \
MatchboxAmplitudellbarqqbarqqbar.h \
MatchboxAmplitudelnuqqbargg.h \
MatchboxAmplitudelnuqqbarg.h \
MatchboxAmplitudelnuqqbar.h \
MatchboxAmplitudelnuqqbarqqbar.h \
MatchboxAmplitudehbbbarg.h \
MatchboxAmplitudehbbbar.h \
MatchboxAmplitudehggg.h \
MatchboxAmplitudehgg.h \
MatchboxAmplitudehqqbarg.h \
MatchboxAmplitudeqqbarttbar.h \
MatchboxAmplitudeqqbarttbarg.h \
MatchboxAmplitudeggttbar.h \
MatchboxAmplitudeggttbarg.h \
HelAmps_sm.h \
Parameters_sm.h \
MG_qqx2ttx.h \
MG_qqx2ttxg.h \
MG_gg2ttx.h \
MG_gg2ttxg.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
MatchboxCurrents.cc \
MatchboxZGammaAmplitude.cc \
MatchboxAmplitudellbarqqbar.cc \
MatchboxAmplitudellbarqqbarg.cc \
MatchboxAmplitudellbarqqbargg.cc \
MatchboxAmplitudellbarqqbarqqbar.cc \
MatchboxAmplitudelnuqqbar.cc \
MatchboxAmplitudelnuqqbarg.cc \
MatchboxAmplitudelnuqqbargg.cc \
MatchboxAmplitudelnuqqbarqqbar.cc \
MatchboxAmplitudehbbbar.cc \
MatchboxAmplitudehbbbarg.cc \
MatchboxAmplitudehgg.cc \
MatchboxAmplitudehggg.cc \
MatchboxAmplitudehqqbarg.cc \
MatchboxAmplitudeqqbarttbar.cc \
MatchboxAmplitudeqqbarttbarg.cc \
MatchboxAmplitudeggttbar.cc \
MatchboxAmplitudeggttbarg.cc \
HelAmps_sm.cc \
Parameters_sm.cc \
MG_qqx2ttx.cc \
MG_qqx2ttxg.cc \
MG_gg2ttx.cc \
MG_gg2ttxg.cc
