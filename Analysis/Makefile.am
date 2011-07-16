pkglib_LTLIBRARIES = HwAnalysis.la
HwAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 9:1:0
HwAnalysis_la_SOURCES = \
EventShapes.cc EventShapes.fh EventShapes.h  \
EventShapesMasterAnalysis.cc EventShapesMasterAnalysis.h \
BasicConsistency.cc BasicConsistency.h \
LEPMultiplicityCount.cc LEPMultiplicityCount.h \
MultiplicityInfo.h \
LEPBMultiplicity.cc\
LEPBMultiplicity.h \
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

pkglib_LTLIBRARIES +=  HwLEPAnalysis.la
HwLEPAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 6:0:0
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
HwLEPJetAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 8:0:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
HwLEPJetAnalysis_la_LIBADD = $(FASTJETLIBS) 
HwLEPJetAnalysis_la_SOURCES = \
LEPJetAnalysis.cc LEPJetAnalysis.h\
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.h \
y23.cc y23.h

#pkglib_LTLIBRARIES += HwPP2HAnalysis.la
#HwPP2HAnalysis_la_LDFLAGS = -module -version-info 4:0:0
#HwPP2HAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
#HwPP2HAnalysis_la_LIBADD = $(KTJETLIBS) $(CLHEPLIB)
#HwPP2HAnalysis_la_SOURCES = \
#PP2HAnalysis.cc PP2HAnalysis.h 

pkglib_LTLIBRARIES += HwVV_ME_Analysis.la
HwVV_ME_Analysis_la_LDFLAGS = -module -version-info 5:0:0
HwVV_ME_Analysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
HwVV_ME_Analysis_la_LIBADD = $(FASTJETLIBS) 
HwVV_ME_Analysis_la_SOURCES = \
VV_ME_Analysis.cc VV_ME_Analysis.h

pkglib_LTLIBRARIES += HwInclHiggsAnalysis.la
HwInclHiggsAnalysis_la_LDFLAGS = -module -version-info 5:0:0
HwInclHiggsAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
HwInclHiggsAnalysis_la_LIBADD = $(FASTJETLIBS) 
HwInclHiggsAnalysis_la_SOURCES = \
InclHiggsAnalysis.cc InclHiggsAnalysis.h
endif

pkglib_LTLIBRARIES += HwTevatronAnalysis.la
HwTevatronAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 1:0:0
HwTevatronAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwTevatronAnalysis_la_SOURCES = \
ZpTRun2.cc ZpTRun2.h \
ZpTRun1.cc ZpTRun1.h \
Wpt.cc Wpt.h \
Zrapidity.cc Zrapidity.h 
