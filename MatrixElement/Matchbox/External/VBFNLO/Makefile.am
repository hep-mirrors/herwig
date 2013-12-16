SUBDIRS = Diagrams
if WANT_LIBFASTJET
SUBDIRS += VBFAnalysis
endif
pkglib_LTLIBRARIES = HwMatchboxVBFNLO.la HwMatchboxVBFNLO2.la

HwMatchboxVBFNLO_la_LDFLAGS = -module -version-info 1:1:0

HwMatchboxVBFNLO_la_CXXFLAGS = $(AM_CXXFLAGS)
HwMatchboxVBFNLO_la_CXXFLAGS += -Df2cFortran

HwMatchboxVBFNLO_la_CFLAGS = $(AM_CFLAGS)
HwMatchboxVBFNLO_la_CFLAGS += -Df2cFortran

HwMatchboxVBFNLO_la_LIBADD = \
Diagrams/libDiagramContainers.la \
-L$(VBFNLOLIBS) -lVBFNLO -lVBFNLOAmplitudes -lVBFNLOUtilities -lVBFNLOLoops -lHELAS -lVBFNLOPDFsets -lVBFNLOPhasespace -lVBFspecific -lgfortran 

HwMatchboxVBFNLO_la_SOURCES = \
VBFNLOAmplitudeBase.h VBFNLOAmplitudeBase.cc \
VBFNLOAmplitudeVVJJNeutralBase.h VBFNLOAmplitudeVVJJNeutralBase.cc \
VBFNLOAmplitudePP2hJetJet.h VBFNLOAmplitudePP2hJetJet.cc \
VBFNLOAmplitudePP2hJetJetJet.h VBFNLOAmplitudePP2hJetJetJet.cc \
VBFNLOAmplitudeqq2hqqg.h VBFNLOAmplitudeqq2hqqg.cc \
VBFNLOAmplitudeqg2hqqq.h VBFNLOAmplitudeqg2hqqq.cc \
VBFNLOCommonBlocks.h VBFNLOCommonBlocks.cc \
VBFNLOPrototypes.h\
VBFNLOMEBase.h VBFNLOMEBase.cc \
VBFNLOMEVVJJNeutralBase.h VBFNLOMEVVJJNeutralBase.cc \
VBFNLOMEPP2hJetJet.h VBFNLOMEPP2hJetJet.cc \
VBFNLOMEPP2hJetJetJet.h VBFNLOMEPP2hJetJetJet.cc \
VBFNLOMEqq2hqqg.h VBFNLOMEqq2hqqg.cc \
VBFNLOMEqg2hqqq.h VBFNLOMEqg2hqqq.cc 

HwMatchboxVBFNLO2_la_LDFLAGS = -module -version-info 1:0:0

HwMatchboxVBFNLO2_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMatchboxVBFNLO2_la_CPPFLAGS += -I$(VBFNLOINCLUDE)

HwMatchboxVBFNLO2_la_LIBADD = \
-L$(VBFNLOLIBS) \
-lVBFNLO -lVBFNLOAmplitudes -lVBFNLOUtilities -lVBFNLOLoops \
-lHELAS -lVBFNLOPDFsets -lVBFNLOPhasespace -lVBFspecific -lgfortran 

HwMatchboxVBFNLO2_la_SOURCES = \
VBFNLOAmplitude.h VBFNLOAmplitude.cc

