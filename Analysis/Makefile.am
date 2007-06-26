pkglib_LTLIBRARIES = HwAnalysis.la HwLEPAnalysis.la HwDecayAnalysis.la HwBSMAnalysis.la
HwAnalysis_la_LDFLAGS = -module -version-info 3:0:0
HwAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS)
HwAnalysis_la_LIBADD =
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
HiggsJetAnalysis.icc HiggsJetAnalysis.cc

if HAVE_HEPMC
HwAnalysis_la_SOURCES += GraphvizPlot.h GraphvizPlot.fh \
GraphvizPlot.icc GraphvizPlot.cc HepMCFile.h HepMCFile.icc HepMCFile.cc
HwAnalysis_la_CPPFLAGS += $(HEPMCINCLUDE)
HwAnalysis_la_LIBADD += $(HEPMCLIBS)
endif


HwLEPAnalysis_la_LDFLAGS = -module -version-info 2:0:0
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

HwDecayAnalysis_la_LDFLAGS = -module
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
OniumToOniumPiPiAnalysis.fh OniumToOniumPiPiAnalysis.cc 

# analysis code which depends on ktjet
if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwLEPJetAnalysis.la
HwLEPJetAnalysis_la_LDFLAGS = -module -version-info 3:0:0
HwLEPJetAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE)
HwLEPJetAnalysis_la_LIBADD = $(KTJETLIBS)
HwLEPJetAnalysis_la_SOURCES = \
LEPFourJetsAnalysis.cc LEPFourJetsAnalysis.fh \
LEPFourJetsAnalysis.h LEPFourJetsAnalysis.icc \
LeptonDalitzAnalysis.cc LeptonDalitzAnalysis.fh \
LeptonDalitzAnalysis.h LeptonDalitzAnalysis.icc \
TopDalitzAnalysis.cc TopDalitzAnalysis.fh \
TopDalitzAnalysis.h TopDalitzAnalysis.icc
endif

HwBSMAnalysis_la_LDFLAGS = -module
HwBSMAnalysis_la_SOURCES = \
SpinCorrAnalysis.cc SpinCorrAnalysis.h \
SpinCorrAnalysis.icc SpinCorrAnalysis.fh
