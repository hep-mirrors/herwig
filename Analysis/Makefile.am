pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la

HwAnalysis_la_LDFLAGS = -module
HwAnalysis_la_SOURCES = \
EvtShapes.fh EvtShapes.icc EvtShapes.cc EvtShapes.h \
EventShapes.cc EventShapes.fh EventShapes.h EventShapes.icc \
EventShapesMasterAnalysis.cc EventShapesMasterAnalysis.fh \
EventShapesMasterAnalysis.h EventShapesMasterAnalysis.icc \
BasicConsistency.cc BasicConsistency.fh \
BasicConsistency.h BasicConsistency.icc \
MultiplicityCount.cc MultiplicityCount.fh \
MultiplicityCount.h MultiplicityCount.icc \
SimpleLHCAnalysis.h SimpleLHCAnalysis.fh \
SimpleLHCAnalysis.icc SimpleLHCAnalysis.cc\
TTbarAnalysis.h TTbarAnalysis.fh \
TTbarAnalysis.icc TTbarAnalysis.cc\
LPairAnalysis.h LPairAnalysis.fh \
LPairAnalysis.icc LPairAnalysis.cc\
GammaGammaAnalysis.h GammaGammaAnalysis.fh \
GammaGammaAnalysis.icc GammaGammaAnalysis.cc\
GammaJetAnalysis.h GammaJetAnalysis.fh \
GammaJetAnalysis.icc GammaJetAnalysis.cc\
HiggsJetAnalysis.h HiggsJetAnalysis.fh \
HiggsJetAnalysis.icc HiggsJetAnalysis.cc\
GraphvizPlot.h GraphvizPlot.fh GraphvizPlot.icc GraphvizPlot.cc

HwLEPAnalysis_la_LDFLAGS = -module
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.fh \
BFragmentationAnalysisHandler.h BFragmentationAnalysisHandler.icc \
SingleParticleAnalysis.cc SingleParticleAnalysis.fh\
SingleParticleAnalysis.h SingleParticleAnalysis.icc \
LEPEventShapes.cc LEPEventShapes.fh LEPEventShapes.h LEPEventShapes.icc \
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.fh \
IdentifiedParticleAnalysis.h IdentifiedParticleAnalysis.icc \
DrellYanDalitzAnalysis.cc DrellYanDalitzAnalysis.fh \
DrellYanDalitzAnalysis.h DrellYanDalitzAnalysis.icc 

# analysis code which depends on ktjet
if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE)
HwLEPJetAnalysis_la_LIBADD = -L$(KTJETPATH)/lib -lKtEvent
HwLEPJetAnalysis_la_SOURCES = \
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.fh \
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc \
LeptonDalitzAnalysis.cc LeptonDalitzAnalysis.fh \
LeptonDalitzAnalysis.h LeptonDalitzAnalysis.icc \
TopDalitzAnalysis.cc TopDalitzAnalysis.fh \
TopDalitzAnalysis.h TopDalitzAnalysis.icc 
endif
