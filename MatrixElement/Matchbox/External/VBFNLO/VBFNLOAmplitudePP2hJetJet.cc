#include "VBFNLOCommonBlocks.h"
#include "VBFNLOAmplitudePP2hJetJet.h"
#include "VBFNLOMEPP2hJetJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOAmplitudePP2hJetJet::VBFNLOAmplitudePP2hJetJet()
  : VBFNLOAmplitudeVVJJNeutralBase(), theDecayChannel(0), theNarrowWidth(true) {
}

VBFNLOAmplitudePP2hJetJet::~VBFNLOAmplitudePP2hJetJet(){}

int VBFNLOAmplitudePP2hJetJet::Njets() const {return 2;}

IBPtr VBFNLOAmplitudePP2hJetJet::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOAmplitudePP2hJetJet::fullclone() const {
  return new_ptr(*this);
}

bool VBFNLOAmplitudePP2hJetJet::allowedProcess(const PDVector& data) const {
  PDPtr Wplus = getParticleData(ParticleID::Wplus); 
  PDPtr Wminus = getParticleData(ParticleID::Wminus); 
  PDPtr Z0 = getParticleData(ParticleID::Z0); 

  assert( data[0]->id() != 21 && data[1]->id() != 21 );
  if ( theCurrent == neutral
       && ( data[0] != data[2]
	    || data[1] != data[3] ) ) return false;
  if ( theCurrent == charged
       &&  ( data[0] != SU2Helper::SU2CC(data[2])
	     || data[1] != SU2Helper::SU2CC(data[3]) ) ) return false;
  if ( ( data[0]->id() < 0 && theIncoming1)
       || (data[1]->id() < 0 && theIncoming2) ) return false;
  if ( ( data[0]->id() > 0 && !theIncoming1)
       || (data[1]->id() > 0 && !theIncoming2) ) return false;
  return true;
}

Ptr<MatchboxMEBase>::ptr VBFNLOAmplitudePP2hJetJet::makeME(const PDVector&) const {
  return new_ptr(VBFNLOMEPP2hJetJet(theCurrent, theIncoming1, theIncoming2, theDecayChannel, theNarrowWidth));
}

void VBFNLOAmplitudePP2hJetJet::prepareMomenta(double (& pbar)[14][4],
					       double (&) [5],int) const{

  L5MomToDouble( (mePartonData()[0]->id() > 0 ? 
  		  meMomenta()[0]  : 
  		  meMomenta()[2]), &pbar[0][0]);

  L5MomToDouble( (mePartonData()[2]->id() > 0 ? 
  		  meMomenta()[2] : 
  		  meMomenta()[0]), &pbar[1][0]);

  L5MomToDouble( (mePartonData()[1]->id() > 0 ? 
  		  meMomenta()[1] : 
  		  meMomenta()[3]), &pbar[2][0]);

  L5MomToDouble( (mePartonData()[3]->id() > 0 ? 
  		  meMomenta()[3] : 
  		  meMomenta()[1]), &pbar[3][0]);

  pbar[4][0]=0;
  pbar[4][1]=0;
  pbar[4][2]=0;
  pbar[4][3]=0;

  if (NDecayProducts() == 1) {
    L5MomToDouble( (meMomenta()[4]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 2) {
    L5MomToDouble( (meMomenta()[4]+meMomenta()[5]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 4) {
    L5MomToDouble( (meMomenta()[4]+meMomenta()[5]+meMomenta()[6]+meMomenta()[7]), &pbar[5][0]);
  }
    
  return;
}

void VBFNLOAmplitudePP2hJetJet::VbfnloMe2(double pbar[14][4],int mePartonSign[14],double [5],int,int nlo,double& uucc,double& uuss,double& ddcc,double& ddss,double& udsc,double& ducs) const{

  initProcess(nlo);
  
  double tree[6];

  QQHQQ(pbar,mePartonSign,nlo,
  	uucc,uuss,ddcc,ddss,udsc,ducs,tree);

  //average over spins
  double fac = 1.0;//4*9;

  // get rid of the feared factor fac from vbfnlo if we use Herwig to decay the Higgs Boson
  if (NDecayProducts() == 1)
    fac *= 16*Constants::pi*BLIPSIVNJ.RMG/(sqr(meMomenta()[4].m2()/GeV2-BLIPSIVNJ.RM2) + sqr(BLIPSIVNJ.RMG));

  //give result in units of sHat
  fac/=pow(lastSHat()/GeV2,NDecayProducts());

  //don't forget the branching ratio or the decay amplitude
  if (NDecayProducts() == 2)
    fac/=BranchingRatio();
  else if (NDecayProducts() == 4){
    double v[4][4];
    double decayAmplitude;
    L5MomToDouble( (meMomenta()[4]), &v[0][0]);
    L5MomToDouble( (meMomenta()[5]), &v[1][0]);
    L5MomToDouble( (meMomenta()[6]), &v[2][0]);
    L5MomToDouble( (meMomenta()[7]), &v[3][0]);
    if (theDecayChannel == 5){
      M2S_VVSUM(v,3,4,2,1,decayAmplitude);
    }
    else if (theDecayChannel == 6)
      M2S_VVSUM(v,2,2,2,2,decayAmplitude);
    else if (theDecayChannel == 7)
      M2S_VVSUM(v,2,2,2,1,decayAmplitude);
    fac /= decayAmplitude/(16 * Constants::pi * BKOPOU.XMG[5]);
  }

  uucc/=fac;
  uuss/=fac;
  ddcc/=fac;
  ddss/=fac;
  udsc/=fac;
  ducs/=fac;

  return;
}

void VBFNLOAmplitudePP2hJetJet::initProcess(const int&) const{

  CGLOBALI.SIGN1 = 1;
  CGLOBALI.SIGN2 = 1;
  if (NDecayProducts() == 1)
    CGLOBALI.N_V = 2;
  else
    CGLOBALI.N_V = NDecayProducts();
  CGLOBALI.N_P = 4;
  CGLOBALI.PS_DIMENSION = 7;
  CGLOBALI.NLO_LOOPS = 1;

  CSCALES.ALS[0][0] = lastAlphaS();
  CSCALES.ALS[0][1] = lastAlphaS();

  QQHQQGI(0);    
}

double VBFNLOAmplitudePP2hJetJet::colourCorrelatedME2(pair<int,int> ij) const {
  if ( ij.first == ij.second ||
       (ij != make_pair(0,2) &&
	ij != make_pair(1,3) &&
	ij != make_pair(2,0) &&
	ij != make_pair(3,1)) ) {
    generator()->logWarning(Exception() 
    			    << "A non-existing colour correlation was requested "
    			    << "from the matrix element '" << name() << "'."
    			    << "Asked for pair (" << ij.first << "," << ij.second
    			    << Exception::warning);
    return 0.0;
  }
  return -me2();
}

double VBFNLOAmplitudePP2hJetJet::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {
  generator()->logWarning(Exception() 
  			  << "A non-exisiting spin correlation was requested "
  			  << "from the matrix element '" << name() << "'."
  			  << Exception::warning);
  return 0.0;
}

bool VBFNLOAmplitudePP2hJetJet::noDipole(int emitter, int spectator) const{
  if (emitter == 0 && spectator == 2) return false;
  if (emitter == 2 && spectator == 0) return false;
  if (emitter == 1 && spectator == 3) return false;
  if (emitter == 3 && spectator == 1) return false;
  return true;
}


void VBFNLOAmplitudePP2hJetJet::doinit(){
  VBFNLOAmplitudeVVJJNeutralBase::doinit();
}

void VBFNLOAmplitudePP2hJetJet::doinitrun() {
  VBFNLOAmplitudeVVJJNeutralBase::doinitrun();
}

// void VBFNLOAmplitudePP2hJetJet::initPSGen() {

//   if (NDecayProducts() != 4) {
//     int bos=6;

//     BLIPSIVNJ.RM2=BKOPOU.XM2[bos-1];
//     BLIPSIVNJ.RMG=BKOPOU.XMG[bos-1];

//     BLIPSIVNJ.RM2MIN=sqrt(BLIPSIVNJ.RM2)*( 1 - 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2 );
//     if (BLIPSIVNJ.RM2MIN > 0) BLIPSIVNJ.RM2MIN = BLIPSIVNJ.RM2MIN*BLIPSIVNJ.RM2MIN;
//     else BLIPSIVNJ.RM2MIN = 0;

//     BLIPSIVNJ.RM2MAX=BLIPSIVNJ.RM2*pow( 1 + 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2, 2);
//     BLIPSIVNJ.RM2MAX = min(BLIPSIVNJ.RM2MAX,pow(sqrt(BLIPSIVNJ.RM2) + 200, 2));

//     BLIPSIVNJ.S=lastS()/GeV2;
//     BLIPSIVNJ.M2MIN=0.001;

//     //these should probably be taken from some cut class eventually
//     BLIPSIVNJ.YJMIN[0]=0;
//     BLIPSIVNJ.YJMIN[1]=0;

//     BLIPSIVNJ.YJMAX[0]=5;
//     BLIPSIVNJ.YJMAX[1]=5;

//     BLIPSIVNJ.PTJMIN[0]=20;
//     BLIPSIVNJ.PTJMIN[1]=20;

//     BLIPSIVNJ.EJMIN[0]=0;
//     BLIPSIVNJ.EJMIN[1]=0;

//     BLIPSIVNJ.INFOJ[0]=-1;
//     BLIPSIVNJ.INFOJ[1]=-1;
//   }
//   else {
//     int bos[3];
//     bos[0] = 6;
//     if (theDecayChannel == 5) {
//       bos[1] = 3;
//       bos[2] = 4;
//     }
//     else {
//       bos[1] = 2;
//       bos[2] = 2;
//     }

//     BLIPSIVVNJ.S = lastS()/GeV2;
//     BLIPSIVVNJ.M2MIN = 0.01;

//     BLIPSIVVNJ.IWIDTH[0] = 0;
//     BLIPSIVVNJ.IWIDTH[1] = 1;
//     BLIPSIVVNJ.IWIDTH[2] = 1;

//     for (int i = 0; i < 3; i++) {
//       BLIPSIVVNJ.RM2[i] = BKOPOU.XM2[bos[i]-1];
//       BLIPSIVVNJ.RMG[i] = BKOPOU.XMG[bos[i]-1];
//     }

//     BLIPSIVVNJ.RM2MIN[0] = sqrt(BLIPSIVVNJ.RM2[0])*( 1 - 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0] );
//     if (BLIPSIVVNJ.RM2MIN[0] > 0) BLIPSIVVNJ.RM2MIN[0] = BLIPSIVVNJ.RM2MIN[0]*BLIPSIVVNJ.RM2MIN[0];
//     else BLIPSIVVNJ.RM2MIN[0] = 0;
//     BLIPSIVVNJ.RM2MAX[0] = BLIPSIVVNJ.RM2[0]*pow( 1 + 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0], 2);
//     BLIPSIVVNJ.RM2MAX[0] = min(BLIPSIVVNJ.RM2MAX[0],BLIPSIVVNJ.S/4);

//     for (int i = 1; i < 3; i++) {
//       if (bos[i] == 2) {
// 	BLIPSIVVNJ.RM2MIN[i] = 225;
// 	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.RM2MAX[0];
//       }
//       else {
// 	BLIPSIVVNJ.RM2MIN[i] = 0.000001;
// 	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.S/4;
//       }
//     }

//     //these should probably be taken from some cut class eventually
//     BLIPSIVVNJ.YJMIN[0] = 0;
//     BLIPSIVVNJ.YJMIN[1] = 0;

//     BLIPSIVVNJ.YJMAX[0] = 5;
//     BLIPSIVVNJ.YJMAX[1] = 5;

//     BLIPSIVVNJ.PTJMIN[0] = 20;
//     BLIPSIVVNJ.PTJMIN[1] = 20;

//     BLIPSIVVNJ.EJMIN[0] = 0;
//     BLIPSIVVNJ.EJMIN[1] = 0;

//     BLIPSIVVNJ.INFOJ[0] = -1;
//     BLIPSIVVNJ.INFOJ[1] = -1;
//   }
//   return;
// }


int VBFNLOAmplitudePP2hJetJet::NDecayProducts() const {
  if (theDecayChannel == 0) return 1;
  if (theDecayChannel < 5) return 2;
  if (theDecayChannel < 8) return 4;
  return 0;
}

double VBFNLOAmplitudePP2hJetJet::BranchingRatio() const {
  if (theDecayChannel == 1) {
    return BRANCH.BHGAM;
  }
  if (theDecayChannel == 2) {
    return BRANCH.BHMU;
  }
  if (theDecayChannel == 3) {
    return BRANCH.BHTAU;
  }
  if (theDecayChannel == 4) {
    return BRANCH.BHBB;
  }
  return 1;
}

ClassDescription<VBFNLOAmplitudePP2hJetJet> VBFNLOAmplitudePP2hJetJet::initVBFNLOAmplitudePP2hJetJet;
// Definition of the static class description member.

void VBFNLOAmplitudePP2hJetJet::persistentOutput(PersistentOStream & os) const {
  os  << theDecayChannel << theNarrowWidth;
}

void VBFNLOAmplitudePP2hJetJet::persistentInput(PersistentIStream & is, int) {
  is  >> theDecayChannel >> theNarrowWidth;
}


void VBFNLOAmplitudePP2hJetJet::Init() {

  static ClassDocumentation<VBFNLOAmplitudePP2hJetJet> documentation
    ("VBFNLOAmplitudePP2hJetJet");

  static Switch<VBFNLOAmplitudePP2hJetJet,int> interfaceDecayChannel
    ("DecayChannel",
     "Choose the decay channel that is simulated by VBFNLO.",
     &VBFNLOAmplitudePP2hJetJet::theDecayChannel, 0, true, false);
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

  static Switch<VBFNLOAmplitudePP2hJetJet,bool> interfaceNarrowWidth
    ("NarrowWidth",
     "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
     &VBFNLOAmplitudePP2hJetJet::theNarrowWidth, true, true, false);
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
