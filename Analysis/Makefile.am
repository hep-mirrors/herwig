pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la

HwAnalysis_la_LDFLAGS = -module
HwAnalysis_la_SOURCES = \
EvtShapes.fh EvtShapes.icc EvtShapes.cc EvtShapes.h \
EventShapes.cc EventShapes.fh EventShapes.h EventShapes.icc \
EventShapesMasterAnalysis.cc EventShapesMasterAnalysis.fh \
EventShapesMasterAnalysis.h EventShapesMasterAnalysis.icc \
Histogram.cc Histogram.fh Histogram.h Histogram.icc \
Statistic.cc Statistic.fh Statistic.h Statistic.icc \
BasicConsistency.cc BasicConsistency.fh \
BasicConsistency.h BasicConsistency.icc \
MultiplicityCount.cc MultiplicityCount.fh \
MultiplicityCount.h MultiplicityCount.icc \
SimpleLHCAnalysis.h SimpleLHCAnalysis.fh \
SimpleLHCAnalysis.icc SimpleLHCAnalysis.cc

HwLEPAnalysis_la_LDFLAGS = -module
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.fh \
BFragmentationAnalysisHandler.h BFragmentationAnalysisHandler.icc \
SingleParticleAnalysis.cc SingleParticleAnalysis.fh\
SingleParticleAnalysis.h SingleParticleAnalysis.icc \
LEPEventShapes.cc LEPEventShapes.fh LEPEventShapes.h LEPEventShapes.icc \
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.fh \
IdentifiedParticleAnalysis.h IdentifiedParticleAnalysis.icc 

# analysis code which depends on ktjet
if WANT_LIBKTJET
AM_CPPFLAGS+= -I$(KTJETPATH)
HwLEPAnalysis_la_LIBADD = -L$(KTJETPATH)/lib -lKtEvent
HwLEPAnalysis_la_SOURCES +=\
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.fh\
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc 
endif
