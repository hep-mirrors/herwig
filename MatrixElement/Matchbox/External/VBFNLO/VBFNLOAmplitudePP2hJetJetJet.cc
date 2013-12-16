#include "VBFNLOCommonBlocks.h"
#include "VBFNLOAmplitudePP2hJetJetJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;


VBFNLOAmplitudePP2hJetJetJet::VBFNLOAmplitudePP2hJetJetJet()
  : VBFNLOAmplitudeVVJJNeutralBase(), theDecayChannel(0), theNarrowWidth(true) {
}

VBFNLOAmplitudePP2hJetJetJet::~VBFNLOAmplitudePP2hJetJetJet(){}

int VBFNLOAmplitudePP2hJetJetJet::Njets() const {return 3;}

void VBFNLOAmplitudePP2hJetJetJet::VbfnloMe2(double pbar[14][4],int mePartonSign[14],double qbar[5],int meGluonSign,int nlo,double& uucc,double& uuss,double& ddcc,double& ddss,double& udsc,double& ducs) const{

  initProcess(nlo);

  QQHQQJ(pbar,mePartonSign,qbar,meGluonSign,
  	uucc,uuss,ddcc,ddss,udsc,ducs);

  // now get rid of the feared factor fac from vbfnlo    
  double fac = crossingSign();
  if (NDecayProducts() == 1)
    fac *= 16*Constants::pi*BLIPSIVNJ.RMG/(sqr(meMomenta()[5].m2()/GeV2-BLIPSIVNJ.RM2) + sqr(BLIPSIVNJ.RMG));

  //average over spins and colors
  // if (meGluonSign == 1) {fac*=4*9;}  //outgoing gluon
  // else fac*=4*24;

  //give result in units of sHat
  fac/=pow(lastSHat()/GeV2,NDecayProducts()+1);

  if (NDecayProducts() == 2)
    fac/=BranchingRatio();
  else if (NDecayProducts() == 4){
    double v[4][4];
    double decayAmplitude;
    L5MomToDouble( (meMomenta()[5]), &v[0][0]);
    L5MomToDouble( (meMomenta()[6]), &v[1][0]);
    L5MomToDouble( (meMomenta()[7]), &v[2][0]);
    L5MomToDouble( (meMomenta()[8]), &v[3][0]);
    if (DecayChannel() == 5)
      M2S_VVSUM(v,3,4,2,1,decayAmplitude);
    else if (DecayChannel() == 6)
      M2S_VVSUM(v,2,2,2,2,decayAmplitude);
    else if (DecayChannel() == 7)
      M2S_VVSUM(v,2,2,2,1,decayAmplitude);
    fac /= decayAmplitude;
  }

  // cerr << "fac= " << fac << "\n" << flush;
  // cerr << uucc << " " << uuss << " " << ddcc << " " << ddss << " " << udsc << " " << ducs << "\n" << flush;

  uucc/=fac;
  uuss/=fac;
  ddcc/=fac;
  ddss/=fac;
  udsc/=fac;
  ducs/=fac;
  
  return;
}

void VBFNLOAmplitudePP2hJetJetJet::initProcess(const int & nlo) const{

  if (mePartonData()[0]->id() == 21) SPLITCB.GLUONID = 0;
  else if (mePartonData()[1]->id() == 21) SPLITCB.GLUONID = 1;
  else SPLITCB.GLUONID = -1;

  if (nlo !=0) {
    generator()->logWarning(Exception() 
			    << "The matrix element '" << name() << "' "
			    << "is not capable of calculating virtual corrections "
			    << Exception::warning);
  }

  CGLOBALI.SIGN1 = 1;
  CGLOBALI.SIGN2 = 1;
  CGLOBALI.N_P = 5;
  if (NDecayProducts() == 1)
    CGLOBALI.N_V = 2;
  else
    CGLOBALI.N_V = NDecayProducts();
  CGLOBALI.NLO_LOOPS = 0;
  QQHQQGI(0);
  
  CSCALES.ALS[0][0] = lastAlphaS();
  CSCALES.ALS[0][1] = lastAlphaS();

}

double VBFNLOAmplitudePP2hJetJetJet::colourCorrelatedME2(pair<int,int>) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  return 0.0;
}

double VBFNLOAmplitudePP2hJetJetJet::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  return 0.0;

}

void VBFNLOAmplitudePP2hJetJetJet::doinit(){
  VBFNLOAmplitudeBase::doinit();
}

int VBFNLOAmplitudePP2hJetJetJet::NDecayProducts() const {
  if (DecayChannel() == 0) return 1;
  if (DecayChannel() < 5) return 2;
  if (DecayChannel() < 8) return 4;
  return 0;
}

double VBFNLOAmplitudePP2hJetJetJet::BranchingRatio() const {
  if (DecayChannel() == 1) {
    return BRANCH.BHGAM;
  }
  if (DecayChannel() == 2) {
    return BRANCH.BHMU;
  }
  if (DecayChannel() == 3) {
    return BRANCH.BHTAU;
  }
  if (DecayChannel() == 4) {
    return BRANCH.BHBB;
  }
  return 1;
}

// bool VBFNLOAmplitudePP2hJetJetJet::noDipole(int emitter, int emission, int spectator) const{
//   if (emission == 4){
//     if (emitter == 0 && spectator == 2) return false;
//     if (emitter == 2 && spectator == 0) return false;
//     if (emitter == 1 && spectator == 3) return false;
//     if (emitter == 3 && spectator == 1) return false;
//   }
//   else if (emission == 2){
//     if (emitter == 0 && spectator == 4) return false;
//     if (emitter == 4 && spectator == 0) return false;
//   }
//   else if (emission == 3){
//     if (emitter == 1 && spectator == 4) return false;
//     if (emitter == 4 && spectator == 1) return false;
//   }
//   return true;
// }

// int VBFNLOAmplitudePP2hJetJetJet::nDim() const {
//   if (phasespace()) return phasespace()->nDim(4+NDecayProducts()-1);
//   if (NDecayProducts() == 1) return 10;
//   if (NDecayProducts() == 2) return 13;
//   if (NDecayProducts() == 4) return 19;
//   return 0;
// }

AbstractClassDescription<VBFNLOAmplitudePP2hJetJetJet> VBFNLOAmplitudePP2hJetJetJet::initVBFNLOAmplitudePP2hJetJetJet;
// Definition of the static class description member.

void VBFNLOAmplitudePP2hJetJetJet::persistentOutput(PersistentOStream & os) const {
  os << theDecayChannel << theNarrowWidth;
}

void VBFNLOAmplitudePP2hJetJetJet::persistentInput(PersistentIStream & is, int) {
  is >> theDecayChannel >> theNarrowWidth;
}

void VBFNLOAmplitudePP2hJetJetJet::Init() {

  static ClassDocumentation<VBFNLOAmplitudePP2hJetJetJet> documentation
    ("VBFNLOAmplitudePP2hJetJetJet");

  static Switch<VBFNLOAmplitudePP2hJetJetJet,int> interfaceDecayChannel
    ("DecayChannel",
     "Choose the decay channel that is simulated by VBFNLO.",
     &VBFNLOAmplitudePP2hJetJetJet::theDecayChannel, 0, true, false);
  static SwitchOption interfaceDecayChannelStable
    (interfaceDecayChannel,
     "Stable",
     "No Higgs decay is calculated within VBFNLO.",
     0);
  static SwitchOption interfaceDecayChannelAA
    (interfaceDecayChannel,
     "H -> A A",
     "Higgs decay into photons",
     1);
  static SwitchOption interfaceDecayChannelMu
    (interfaceDecayChannel,
     "H -> mu mu",
     "Higgs decay into muons",
     2);
  static SwitchOption interfaceDecayChannelTau
    (interfaceDecayChannel,
     "H -> tau tau",
     "Higgs decay into taus",
     3);
  static SwitchOption interfaceDecayChannelBBar
    (interfaceDecayChannel,
     "H -> b bbar",
     "Higgs decay into b anti-b",
     4);
  static SwitchOption interfaceDecayChannelWW
    (interfaceDecayChannel,
     "H -> W+ W-",
     "Higgs decay into W bosons",
     5);
  static SwitchOption interfaceDecayChannelZZ_ll
    (interfaceDecayChannel,
     "H -> Z Z -> l lbar",
     "Higgs decay into Z bosons into lepton antilepton",
     6);
  static SwitchOption interfaceDecayChannelZZ_lnu
    (interfaceDecayChannel,
     "H -> Z Z -> l nu",
     "Higgs decay into Z bosons into lepton neutrino",
     7);

  static Switch<VBFNLOAmplitudePP2hJetJetJet,bool> interfaceNarrowWidth
    ("NarrowWidth",
     "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
     &VBFNLOAmplitudePP2hJetJetJet::theNarrowWidth, true, true, false);
  static SwitchOption interfaceNarrowWidthTrue
    (interfaceNarrowWidth,
     "True",
     "Calculate with narrow width approximation",
     true);
  static SwitchOption interfaceNarrowWidthFalse
    (interfaceNarrowWidth,
     "False",
     "No narrow width approximation",
     false);

}
