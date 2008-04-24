SUBDIRS = BSM

pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la HwDecayAnalysis.la \
                     HwMPIAnalysis.la
HwAnalysis_la_LDFLAGS = -module -version-info 5:0:0
HwAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwAnalysis_la_LIBADD =
HwAnalysis_la_SOURCES = \
EventShapes.cc EventShapes.fh EventShapes.h EventShapes.icc \
EventShapesMasterAnalysis.cc EventShapesMasterAnalysis.fh \
EventShapesMasterAnalysis.h EventShapesMasterAnalysis.icc \
BasicConsistency.cc BasicConsistency.fh \
BasicConsistency.h BasicConsistency.icc \
MultiplicityCount.cc MultiplicityCount.fh \
MultiplicityCount.h MultiplicityCount.icc \
LEPBMultiplicity.cc LEPBMultiplicity.fh \
LEPBMultiplicity.h LEPBMultiplicity.icc \
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
HiggsJetAnalysis.icc HiggsJetAnalysis.cc \
DrellYanPT.h DrellYanPT.fh DrellYanPT.icc DrellYanPT.cc

pkglib_LTLIBRARIES += HwBbarAnalysis.la
HwBbarAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwBbarAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwBbarAnalysis_la_SOURCES = \
BbarAnalysis.cc BbarAnalysis.fh \
BbarAnalysis.h BbarAnalysis.icc

if HAVE_HEPMC
HwAnalysis_la_SOURCES += GraphvizPlot.h GraphvizPlot.fh \
GraphvizPlot.icc GraphvizPlot.cc HepMCFile.h HepMCFile.cc
HwAnalysis_la_CPPFLAGS += $(HEPMCINCLUDE)
HwAnalysis_la_LIBADD += $(HEPMCLIBS)
endif


HwLEPAnalysis_la_LDFLAGS = -module -version-info 3:0:0
HwLEPAnalysis_la_SOURCES = \
BFragmentationAnalysisHandler.cc BFragmentationAnalysisHandler.fh \
BFragmentationAnalysisHandler.h BFragmentationAnalysisHandler.icc \
SingleParticleAnalysis.cc SingleParticleAnalysis.fh\
SingleParticleAnalysis.h SingleParticleAnalysis.icc \
LEPEventShapes.cc LEPEventShapes.fh LEPEventShapes.h LEPEventShapes.icc \
IdentifiedParticleAnalysis.cc IdentifiedParticleAnalysis.fh \
IdentifiedParticleAnalysis.h IdentifiedParticleAnalysis.icc \
DrellYanDalitzAnalysis.cc DrellYanDalitzAnalysis.fh \
DrellYanDalitzAnalysis.h DrellYanDalitzAnalysis.icc \
BELLECharmAnalysis.h BELLECharmAnalysis.icc \
BELLECharmAnalysis.fh BELLECharmAnalysis.cc \
CLEOCharmAnalysis.h CLEOCharmAnalysis.icc \
CLEOCharmAnalysis.fh CLEOCharmAnalysis.cc

HwDecayAnalysis_la_LDFLAGS = -module -version-info 1:0:0
HwDecayAnalysis_la_SOURCES = \
TauTo2MesonAnalysis.h TauTo2MesonAnalysis.icc \
TauTo2MesonAnalysis.fh TauTo2MesonAnalysis.cc \
TauToLeptonsAnalysis.h TauToLeptonsAnalysis.icc \
TauToLeptonsAnalysis.fh TauToLeptonsAnalysis.cc \
TauTo3MesonAnalysis.h TauTo3MesonAnalysis.icc \
TauTo3MesonAnalysis.fh TauTo3MesonAnalysis.cc \
TauTo4MesonAnalysis.h TauTo4MesonAnalysis.icc \
TauTo4MesonAnalysis.fh TauTo4MesonAnalysis.cc \
Tau5Pion.h Tau5Pion.icc \
Tau5Pion.fh Tau5Pion.cc \
TauCorrelationAnalysis.h TauCorrelationAnalysis.icc \
TauCorrelationAnalysis.fh TauCorrelationAnalysis.cc \
EtaDecayAnalysis.h EtaDecayAnalysis.icc \
EtaDecayAnalysis.fh EtaDecayAnalysis.cc \
SemiLeptonicDecayAnalysis.h SemiLeptonicDecayAnalysis.icc \
SemiLeptonicDecayAnalysis.fh SemiLeptonicDecayAnalysis.cc \
OmegaPhi3PionAnalysis.h OmegaPhi3PionAnalysis.icc \
OmegaPhi3PionAnalysis.fh OmegaPhi3PionAnalysis.cc \
VffGammaAnalysis.h VffGammaAnalysis.icc \
VffGammaAnalysis.fh VffGammaAnalysis.cc  \
VPPGammaAnalysis.h VPPGammaAnalysis.icc \
VPPGammaAnalysis.fh VPPGammaAnalysis.cc  \
PScalarVectorFermionsAnalysis.h PScalarVectorFermionsAnalysis.icc \
PScalarVectorFermionsAnalysis.fh PScalarVectorFermionsAnalysis.cc \
VectorPScalarFermionsAnalysis.h VectorPScalarFermionsAnalysis.icc \
VectorPScalarFermionsAnalysis.fh VectorPScalarFermionsAnalysis.cc \
Pi4ElectronAnalysis.h Pi4ElectronAnalysis.icc \
Pi4ElectronAnalysis.fh Pi4ElectronAnalysis.cc  \
BtoSGammaAnalysis.h BtoSGammaAnalysis.icc \
BtoSGammaAnalysis.fh BtoSGammaAnalysis.cc   \
ZPhotonsAnalysis.h ZPhotonsAnalysis.icc \
ZPhotonsAnalysis.fh ZPhotonsAnalysis.cc \
DDalitzAnalysis.h DDalitzAnalysis.icc \
DDalitzAnalysis.fh DDalitzAnalysis.cc  \
SemiLeptonicDPiAnalysis.h SemiLeptonicDPiAnalysis.icc \
SemiLeptonicDPiAnalysis.fh SemiLeptonicDPiAnalysis.cc  \
OniumToOniumPiPiAnalysis.h OniumToOniumPiPiAnalysis.icc \
OniumToOniumPiPiAnalysis.fh OniumToOniumPiPiAnalysis.cc\
BeamCorrelationAnalysis.h BeamCorrelationAnalysis.icc \
BeamCorrelationAnalysis.fh BeamCorrelationAnalysis.cc  \
a1DecayAnalysis.h a1DecayAnalysis.icc \
a1DecayAnalysis.fh a1DecayAnalysis.cc  \
BranchingRatioAnalysis.h BranchingRatioAnalysis.icc \
BranchingRatioAnalysis.fh BranchingRatioAnalysis.cc

# analysis code which depends on ktjet
if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
HwLEPJetAnalysis_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwLEPJetAnalysis_la_SOURCES = \
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.fh \
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc \
LEPJetAnalysis.cc LEPJetAnalysis.fh \
LEPJetAnalysis.h LEPJetAnalysis.icc \
LeptonDalitzAnalysis.cc LeptonDalitzAnalysis.fh \
LeptonDalitzAnalysis.h LeptonDalitzAnalysis.icc \
TopDalitzAnalysis.cc TopDalitzAnalysis.fh \
TopDalitzAnalysis.h TopDalitzAnalysis.icc

#code which also depends on ROOT
if WANT_LIBROOT
pkglib_LTLIBRARIES += HwUEAnalysis.la
HwUEAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwUEAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) $(ROOTINCLUDE)
HwUEAnalysis_la_LIBADD = $(KTJETLIBS) $(ROOTLIBS)
HwUEAnalysis_la_SOURCES = \
UEobservables.cc UEobservables.fh \
UEobservables.h UEobservables.icc 
endif
endif

HwMPIAnalysis_la_LDFLAGS = -module -version-info 4:0:0
HwMPIAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMPIAnalysis_la_LIBADD =
HwMPIAnalysis_la_SOURCES = \
RFieldAnalysis.h \
RFieldAnalysis.cc RFieldAnalysis.icc \
MPIMultiplicity.h MPIMultiplicity.fh \
MPIMultiplicity.icc MPIMultiplicity.cc

pkglib_LTLIBRARIES += HwZpt.la
HwZpt_la_LDFLAGS = -module -version-info 4:0:0
HwZpt_la_CPPFLAGS = $(AM_CPPFLAGS)
HwZpt_la_SOURCES = \
ZpTRun2.cc ZpTRun2.h \
ZpTRun1.cc ZpTRun1.h 


pkglib_LTLIBRARIES += HwZrap.la
HwZrap_la_LDFLAGS = -module -version-info 4:0:0
HwZrap_la_CPPFLAGS = $(AM_CPPFLAGS)
HwZrap_la_SOURCES = \
Zrapidity.cc Zrapidity.h 