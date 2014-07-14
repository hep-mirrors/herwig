pkglib_LTLIBRARIES = HwMatchboxVBFNLO.la

HwMatchboxVBFNLO_la_LDFLAGS = -module -version-info 1:0:0

HwMatchboxVBFNLO_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMatchboxVBFNLO_la_CPPFLAGS += -I$(VBFNLOINCLUDE)

#HwMatchboxVBFNLO_la_LIBADD = \
#-L$(VBFNLOLIBS) \
#-lVBFNLO -lVBFNLOAmplitudes -lVBFNLOUtilities -lVBFNLOLoops \
#-lHELAS -lVBFNLOPDFsets -lVBFNLOPhasespace -lVBFspecific -lgfortran 

HwMatchboxVBFNLO_la_SOURCES = \
VBFNLOAmplitude.h VBFNLOAmplitude.cc \
VBFNLOPhasespace.h VBFNLOPhasespace.cc
