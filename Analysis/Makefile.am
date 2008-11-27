SUBDIRS = BSM

pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la \
                     HwMPIAnalysis.la HwDrellYanDalitzAnalysis.la
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
DrellYanPT.h DrellYanPT.cc\
LesHouchesOutput.h LesHouchesOutput.cc

pkglib_LTLIBRARIES += HwBbarAnalysis.la
HwBbarAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwBbarAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwBbarAnalysis_la_SOURCES = \
BbarAnalysis.cc BbarAnalysis.fh \
BbarAnalysis.h BbarAnalysis.icc

pkglib_LTLIBRARIES += HwVV_ME_Analysis.la
HwVV_ME_Analysis_la_LDFLAGS = -module -version-info 4:0:0
HwVV_ME_Analysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwVV_ME_Analysis_la_SOURCES = \
VV_ME_Analysis.cc VV_ME_Analysis.fh \
VV_ME_Analysis.h VV_ME_Analysis.icc

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
BELLECharmAnalysis.h BELLECharmAnalysis.cc\
CLEOCharmAnalysis.h CLEOCharmAnalysis.cc

HwDrellYanDalitzAnalysis_la_LDFLAGS = -module -version-info 1:0:0
HwDrellYanDalitzAnalysis_la_SOURCES = \
DrellYanDalitzAnalysis.cc DrellYanDalitzAnalysis.h 

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
endif

HwMPIAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwMPIAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMPIAnalysis_la_LIBADD =
HwMPIAnalysis_la_SOURCES = \
RFieldAnalysis.h RFieldAnalysis.cc \
MPIMultiplicity.h MPIMultiplicity.cc

pkglib_LTLIBRARIES += HwZpt.la
HwZpt_la_LDFLAGS = -module -version-info 4:0:0
HwZpt_la_CPPFLAGS = $(AM_CPPFLAGS)
HwZpt_la_SOURCES = \
ZpTRun2.cc ZpTRun2.h \
ZpTRun1.cc ZpTRun1.h 

pkglib_LTLIBRARIES += HwWpt.la
HwWpt_la_LDFLAGS = -module -version-info 4:0:0
HwWpt_la_CPPFLAGS = $(AM_CPPFLAGS)
HwWpt_la_SOURCES = \
Wpt.cc Wpt.h 


pkglib_LTLIBRARIES += HwZrap.la
HwZrap_la_LDFLAGS = -module -version-info 4:0:0
HwZrap_la_CPPFLAGS = $(AM_CPPFLAGS)
HwZrap_la_SOURCES = \
Zrapidity.cc Zrapidity.h 
