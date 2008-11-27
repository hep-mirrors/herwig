SUBDIRS = BSM

pkglib_LTLIBRARIES = HwAnalysis.la
HwAnalysis_la_LDFLAGS = -module -version-info 7:0:0
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

pkglib_LTLIBRARIES +=  HwLEPAnalysis.la
HwLEPAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.h\
SingleParticleAnalysis.cc SingleParticleAnalysis.h\
LEPEventShapes.cc LEPEventShapes.h\
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.h\
BELLECharmAnalysis.h BELLECharmAnalysis.cc\
CLEOCharmAnalysis.h CLEOCharmAnalysis.cc

if HAVE_HEPMC
pkglib_LTLIBRARIES += HwHepMCAnalysis.la
HwHepMCAnalysis_la_LDFLAGS = -module -version-info 1:0:0
HwHepMCAnalysis_la_SOURCES = \
GraphvizPlot.h GraphvizPlot.cc \
HepMCFile.h    HepMCFile.cc
HwHepMCAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(HEPMCINCLUDE)
HwHepMCAnalysis_la_LIBADD = $(HEPMCLIBS)
endif


# analysis code which depends on ktjet
if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module -version-info 5:0:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
HwLEPJetAnalysis_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwLEPJetAnalysis_la_SOURCES = \
LEPFourJetsAnalysis.cc\
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc \
LEPJetAnalysis.cc LEPJetAnalysis.h
endif
