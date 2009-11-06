pkglib_LTLIBRARIES = HwAnalysis.la
HwAnalysis_la_LDFLAGS = -module -version-info 8:0:0
HwAnalysis_la_SOURCES = \
EventShapes.cc EventShapes.fh EventShapes.h  \
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
HwLEPAnalysis_la_LDFLAGS = -module -version-info 5:0:0
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.h\
SingleParticleAnalysis.cc SingleParticleAnalysis.h\
LEPEventShapes.cc LEPEventShapes.h\
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.h\
BELLECharmAnalysis.h BELLECharmAnalysis.cc\
CLEOCharmAnalysis.h CLEOCharmAnalysis.cc

# analysis code which depends on fastjet
if WANT_LIBFASTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module -version-info 6:0:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
HwLEPJetAnalysis_la_LIBADD = $(FASTJETLIBS) 
HwLEPJetAnalysis_la_SOURCES = \
LEPJetAnalysis.cc LEPJetAnalysis.h\
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.h
endif
