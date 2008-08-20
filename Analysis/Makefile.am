SUBDIRS = BSM

pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la \
                     HwMPIAnalysis.la
HwAnalysis_la_LDFLAGS = -module -version-info 6:0:0
HwAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwAnalysis_la_LIBADD =
HwAnalysis_la_SOURCES = \
EventShapes.cc EventShapes.fh EventShapes.h EventShapes.icc \
EventShapesMasterAnalysis.cc EventShapesMasterAnalysis.h \
BasicConsistency.cc BasicConsistency.h \
LEPMultiplicityCount.cc LEPMultiplicityCount.h \
MultiplicityInfo.h \
LEPBMultiplicity.cc\
LEPBMultiplicity.h LEPBMultiplicity.icc \
SimpleLHCAnalysis.h SimpleLHCAnalysis.cc\
TTbarAnalysis.h TTbarAnalysis.cc\
LPairAnalysis.h LPairAnalysis.cc\
GammaGammaAnalysis.h GammaGammaAnalysis.cc\
GammaJetAnalysis.h GammaJetAnalysis.cc\
HiggsJetAnalysis.h HiggsJetAnalysis.cc \
DrellYanPT.h DrellYanPT.cc

if HAVE_HEPMC
HwAnalysis_la_SOURCES +=  \
GraphvizPlot.h GraphvizPlot.cc \
HepMCFile.h    HepMCFile.cc
HwAnalysis_la_CPPFLAGS += $(HEPMCINCLUDE)
HwAnalysis_la_LIBADD += $(HEPMCLIBS)
endif


HwLEPAnalysis_la_LDFLAGS = -module -version-info 3:1:0
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.h\
SingleParticleAnalysis.cc SingleParticleAnalysis.h\
LEPEventShapes.cc LEPEventShapes.h\
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.h\
DrellYanDalitzAnalysis.cc DrellYanDalitzAnalysis.h \
BELLECharmAnalysis.h BELLECharmAnalysis.cc\
CLEOCharmAnalysis.h CLEOCharmAnalysis.cc

# analysis code which depends on ktjet
if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module -version-info 4:1:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
HwLEPJetAnalysis_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwLEPJetAnalysis_la_SOURCES = \
LEPFourJetsAnalysis.cc\
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc \
LEPJetAnalysis.cc LEPJetAnalysis.h\
LeptonDalitzAnalysis.cc LeptonDalitzAnalysis.h \
TopDalitzAnalysis.cc TopDalitzAnalysis.h

#code which also depends on ROOT
if WANT_LIBROOT
pkglib_LTLIBRARIES += HwUEAnalysis.la
HwUEAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwUEAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) $(ROOTINCLUDE)
HwUEAnalysis_la_LIBADD = $(KTJETLIBS) $(ROOTLIBS)
HwUEAnalysis_la_SOURCES = \
UEobservables.cc UEobservables.h
endif
endif

HwMPIAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwMPIAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMPIAnalysis_la_LIBADD =
HwMPIAnalysis_la_SOURCES = \
RFieldAnalysis.h RFieldAnalysis.cc \
MPIMultiplicity.h MPIMultiplicity.cc
