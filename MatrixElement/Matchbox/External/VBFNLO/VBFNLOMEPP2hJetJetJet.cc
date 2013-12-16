#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEPP2hJetJetJet.h"
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


VBFNLOMEPP2hJetJetJet::VBFNLOMEPP2hJetJetJet()
  : VBFNLOMEVVJJNeutralBase(), theDecayChannel(0), theNarrowWidth(true) {
}

VBFNLOMEPP2hJetJetJet::VBFNLOMEPP2hJetJetJet(int current, bool incoming1, bool incoming2, int decay, bool nw)
  : VBFNLOMEVVJJNeutralBase(current, incoming1, incoming2), theDecayChannel(decay), theNarrowWidth(nw) {}

VBFNLOMEPP2hJetJetJet::~VBFNLOMEPP2hJetJetJet(){}

int VBFNLOMEPP2hJetJetJet::Njets() const {return 3;}



void VBFNLOMEPP2hJetJetJet::VbfnloMe2(double pbar[14][4],int mePartonSign[14],double qbar[5],int meGluonSign,int nlo,double& uucc,double& uuss,double& ddcc,double& ddss,double& udsc,double& ducs) const{

  initProcess(nlo);

  QQHQQJ(pbar,mePartonSign,qbar,meGluonSign,
  	uucc,uuss,ddcc,ddss,udsc,ducs);

  // now get rid of the feared factor fac from vbfnlo    
  double fac = 1;
  if (NDecayProducts() == 1)
    fac = 16*Constants::pi*BLIPSIVNJ.RMG/(sqr(lastMEMomenta()[5].m2()/GeV2-BLIPSIVNJ.RM2) + sqr(BLIPSIVNJ.RMG));

  //average over spins and colors
  if (meGluonSign == 1) {fac*=4*9;}  //outgoing gluon
  else fac*=4*24;

  //give result in units of sHat
  fac/=pow(lastSHat()/GeV2,NDecayProducts()+1);

  if (NDecayProducts() == 2)
    fac/=BranchingRatio();
  else if (NDecayProducts() == 4){
    double v[4][4];
    double decayAmplitude;
    L5MomToDouble( (lastMEMomenta()[5]), &v[0][0]);
    L5MomToDouble( (lastMEMomenta()[6]), &v[1][0]);
    L5MomToDouble( (lastMEMomenta()[7]), &v[2][0]);
    L5MomToDouble( (lastMEMomenta()[8]), &v[3][0]);
    if (DecayChannel() == 5)
      M2S_VVSUM(v,3,4,2,1,decayAmplitude);
    else if (DecayChannel() == 6)
      M2S_VVSUM(v,2,2,2,2,decayAmplitude);
    else if (DecayChannel() == 7)
      M2S_VVSUM(v,2,2,2,1,decayAmplitude);
    fac /= decayAmplitude;
  }

  uucc/=fac;
  uuss/=fac;
  ddcc/=fac;
  ddss/=fac;
  udsc/=fac;
  ducs/=fac;
  
  return;
}

void VBFNLOMEPP2hJetJetJet::initProcess(const int & nlo) const{

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

double VBFNLOMEPP2hJetJetJet::colourCorrelatedME2(pair<int,int>) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();
}

double VBFNLOMEPP2hJetJetJet::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

void VBFNLOMEPP2hJetJetJet::doinit(){
  VBFNLOMEBase::doinit();
  for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	q != theQuarkFlavours.end(); ++q ){
    if ( (**q).mass() != ZERO )
      Throw<InitException>() << "The matrix element '"
			     << name() << "' is only capable of "
			     << "producing massless quarks.";
    if ( (**q).id() < 1 || (**q).id() > 6) Throw<InitException>() << "Please insert  "
								  << "only quarks in the ME interface of" 
								  << name() <<".";				       
  }
  return;
}


void VBFNLOMEPP2hJetJetJet::initPSGen() {
 
  if (NDecayProducts() != 4) {
    int bos=6;
    
    BLIPSIVNJ.RM2=BKOPOU.XM2[bos-1];
    BLIPSIVNJ.RMG=BKOPOU.XMG[bos-1];
    
    BLIPSIVNJ.RM2MIN=sqrt(BLIPSIVNJ.RM2)*( 1 - 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2 );
    if (BLIPSIVNJ.RM2MIN > 0) BLIPSIVNJ.RM2MIN = BLIPSIVNJ.RM2MIN*BLIPSIVNJ.RM2MIN;
    else BLIPSIVNJ.RM2MIN = 0;
    
    BLIPSIVNJ.RM2MAX=BLIPSIVNJ.RM2*pow( 1 + 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2, 2);
    BLIPSIVNJ.RM2MAX = min(BLIPSIVNJ.RM2MAX,pow(sqrt(BLIPSIVNJ.RM2) + 200, 2));
    
    BLIPSIVNJ.S=lastS()/GeV2;
    BLIPSIVNJ.M2MIN=0.001;
    
    //these should probably be taken from some cut class eventually
    for (int i = 0; i < 4; i++){
      BLIPSIVNJ.YJMIN[i]=0;
      BLIPSIVNJ.YJMAX[i]=5;
      BLIPSIVNJ.PTJMIN[i]=10;
      BLIPSIVNJ.EJMIN[i]=0;
      BLIPSIVNJ.INFOJ[i]=-1;
    }
  }
  else {
    int bos[3];
    bos[0] = 6;
    if (DecayChannel() == 5) {
      bos[1] = 3;
      bos[2] = 4;
    }
    else {
      bos[1] = 2;
      bos[2] = 2;
    }

    BLIPSIVVNJ.S = lastS()/GeV2;
    BLIPSIVVNJ.M2MIN = 0.01;

    BLIPSIVVNJ.IWIDTH[0] = 0;
    BLIPSIVVNJ.IWIDTH[1] = 1;
    BLIPSIVVNJ.IWIDTH[2] = 1;

    for (int i = 0; i < 3; i++) {
      BLIPSIVVNJ.RM2[i] = BKOPOU.XM2[bos[i]-1];
      BLIPSIVVNJ.RMG[i] = BKOPOU.XMG[bos[i]-1];
    }

    BLIPSIVVNJ.RM2MIN[0] = sqrt(BLIPSIVVNJ.RM2[0])*( 1 - 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0] );
    if (BLIPSIVVNJ.RM2MIN[0] > 0) BLIPSIVVNJ.RM2MIN[0] = BLIPSIVVNJ.RM2MIN[0]*BLIPSIVVNJ.RM2MIN[0];
    else BLIPSIVVNJ.RM2MIN[0] = 0;
    BLIPSIVVNJ.RM2MAX[0] = BLIPSIVVNJ.RM2[0]*pow( 1 + 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0], 2);
    BLIPSIVVNJ.RM2MAX[0] = min(BLIPSIVVNJ.RM2MAX[0],BLIPSIVVNJ.S/4);

    for (int i = 1; i < 3; i++) {
      if (bos[i] == 2) {
	BLIPSIVVNJ.RM2MIN[i] = 225;
	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.RM2MAX[0];
      }
      else {
	BLIPSIVVNJ.RM2MIN[i] = 0.000001;
	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.S/4;
      }
    }
    //these should probably be taken from some cut class eventually
    for (int i = 0; i < 4; i++) {
      BLIPSIVVNJ.YJMIN[i] = 0;
      BLIPSIVVNJ.YJMAX[i] = 5;
      BLIPSIVVNJ.PTJMIN[i] = 20;
      BLIPSIVVNJ.EJMIN[i] = 0;
      BLIPSIVVNJ.INFOJ[i] = -1;
    }

  }

  return;
}



bool VBFNLOMEPP2hJetJetJet::generateKinematics(const double * r){

  if ( phasespace() ) {
    initPSGen();
    return MatchboxMEBase::generateKinematics(r);
  }

  int Fn = Njets();
  double Frd[100], Frn, Fk1[4],Fk2[4],Fq[5],Fd[10][4],Fd1[4],Fd2[4],Fp[4][4],Fx1,Fx2,Fw; //check dimensions of Fp
  int Fnw = 0;
  if (NarrowWidth()) Fnw =1;

  for (int i=0; i<nDim()-1; i++){
    Frd[i]=r[i];
  }
  Frn= r[nDim()-1];
  
  for (int i=0; i<nDim()-1; i++){
    //    cerr << "R(" << i+1 << ")=" << r[i] << "\n" << flush;
  }
  //  cerr << "RN=" << r[nDim()-1] << "\n" << flush;
  
  initPSGen();  
  if (NDecayProducts() == 1) {
    LIPSN0(Fn,Frd,Frn,Fk1,Fk2,Fq,Fp,Fx1,Fx2,Fw);
  }
  else if (NDecayProducts() == 2) {
    LIPSN(Fn,Frd,Frn,Fk1,Fk2,Fq,Fd1,Fd2,Fp,Fx1,Fx2,Fw,Fnw);
  }
  else if (NDecayProducts() == 4) {
    LIPSNVV(Fn,Frd,Frn,Fk1,Fk2,Fd,Fp,Fx1,Fx2,Fw);
  }
  if(Fw > 0 ){

    lastMEMomenta()[0]=(Lorentz5Momentum(Fk1[1]*GeV,Fk1[2]*GeV,Fk1[3]*GeV,Fk1[0]*GeV));
    lastMEMomenta()[1]=(Lorentz5Momentum(Fk2[1]*GeV,Fk2[2]*GeV,Fk2[3]*GeV,Fk2[0]*GeV));
    for (int i=0; i<Fn; i++){
      lastMEMomenta()[2+i]=(Lorentz5Momentum(Fp[i][1]*GeV,Fp[i][2]*GeV,Fp[i][3]*GeV,Fp[i][0]*GeV));
    }   
    if (NDecayProducts() == 1)
      lastMEMomenta()[5]=(Lorentz5Momentum(Fq[1]*GeV,Fq[2]*GeV,Fq[3]*GeV,Fq[0]*GeV));    
    else if (NDecayProducts() == 2) {
      lastMEMomenta()[5]=(Lorentz5Momentum(Fd1[1]*GeV,Fd1[2]*GeV,Fd1[3]*GeV,Fd1[0]*GeV));  
      lastMEMomenta()[6]=(Lorentz5Momentum(Fd2[1]*GeV,Fd2[2]*GeV,Fd2[3]*GeV,Fd2[0]*GeV));  
    }
    else if (NDecayProducts() == 4) {   
      lastMEMomenta()[5]=(Lorentz5Momentum(Fd[0][1]*GeV,Fd[0][2]*GeV,Fd[0][3]*GeV,Fd[0][0]*GeV));  
      lastMEMomenta()[6]=(Lorentz5Momentum(Fd[1][1]*GeV,Fd[1][2]*GeV,Fd[1][3]*GeV,Fd[1][0]*GeV));  
      lastMEMomenta()[7]=(Lorentz5Momentum(Fd[2][1]*GeV,Fd[2][2]*GeV,Fd[2][3]*GeV,Fd[2][0]*GeV));  
      lastMEMomenta()[8]=(Lorentz5Momentum(Fd[3][1]*GeV,Fd[3][2]*GeV,Fd[3][3]*GeV,Fd[3][0]*GeV));
    }

    //revert the multiplication with (hbar*c)^2/2   
    //division by sHat also takes place in LIPSN, need to divide by 
    //one more power of SHat
    //divide by Feynman x-es
    Energy2 thisSHat=(meMomenta()[0] + meMomenta()[1]).m2();
    Fw/=pow(thisSHat/GeV2,NDecayProducts());

    jacobian(Fw*2/3.8937E11/Fx1/Fx2);
    setScale();
    logGenerateKinematics(r);
    return true; 
  }
  else { 
    jacobian(0.0);
    return false;
  }
}

int VBFNLOMEPP2hJetJetJet::NDecayProducts() const {
  if (DecayChannel() == 0) return 1;
  if (DecayChannel() < 5) return 2;
  if (DecayChannel() < 8) return 4;
  return 0;
}

double VBFNLOMEPP2hJetJetJet::BranchingRatio() const {
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

bool VBFNLOMEPP2hJetJetJet::noDipole(int emitter, int emission, int spectator) const{
  if (emission == 4){
    if (emitter == 0 && spectator == 2) return false;
    if (emitter == 2 && spectator == 0) return false;
    if (emitter == 1 && spectator == 3) return false;
    if (emitter == 3 && spectator == 1) return false;
  }
  else if (emission == 2){
    if (emitter == 0 && spectator == 4) return false;
    if (emitter == 4 && spectator == 0) return false;
  }
  else if (emission == 3){
    if (emitter == 1 && spectator == 4) return false;
    if (emitter == 4 && spectator == 1) return false;
  }
  return true;
}

int VBFNLOMEPP2hJetJetJet::nDim() const {
  if (phasespace()) return phasespace()->nDim(4+NDecayProducts()-1);
  if (NDecayProducts() == 1) return 10;
  if (NDecayProducts() == 2) return 13;
  if (NDecayProducts() == 4) return 19;
  return 0;
}

AbstractClassDescription<VBFNLOMEPP2hJetJetJet> VBFNLOMEPP2hJetJetJet::initVBFNLOMEPP2hJetJetJet;
// Definition of the static class description member.

void VBFNLOMEPP2hJetJetJet::persistentOutput(PersistentOStream & os) const {
  os << theDecayChannel << theNarrowWidth;
}

void VBFNLOMEPP2hJetJetJet::persistentInput(PersistentIStream & is, int) {
  is >> theDecayChannel >> theNarrowWidth;
}

void VBFNLOMEPP2hJetJetJet::Init() {

  static ClassDocumentation<VBFNLOMEPP2hJetJetJet> documentation
    ("VBFNLOMEPP2hJetJetJet");

  // static Switch<VBFNLOMEPP2hJetJetJet,int> interfaceDecayChannel
  //   ("DecayChannel",
  //    "Choose the decay channel that is simulated by VBFNLO.",
  //    &VBFNLOMEPP2hJetJetJet::theDecayChannel, 0, true, false);
  // static SwitchOption interfaceDecayChannelStable
  //   (interfaceDecayChannel,
  //    "Stable",
  //    "No Higgs decay is calculated within VBFNLO.",
  //    0);
  // static SwitchOption interfaceDecayChannelAA
  //   (interfaceDecayChannel,
  //    "H -> A A",
  //    "Higgs decay into photons",
  //    1);
  // static SwitchOption interfaceDecayChannelMu
  //   (interfaceDecayChannel,
  //    "H -> mu mu",
  //    "Higgs decay into muons",
  //    2);
  // static SwitchOption interfaceDecayChannelTau
  //   (interfaceDecayChannel,
  //    "H -> tau tau",
  //    "Higgs decay into taus",
  //    3);
  // static SwitchOption interfaceDecayChannelBBar
  //   (interfaceDecayChannel,
  //    "H -> b bbar",
  //    "Higgs decay into b anti-b",
  //    4);
  // static SwitchOption interfaceDecayChannelWW
  //   (interfaceDecayChannel,
  //    "H -> W+ W-",
  //    "Higgs decay into W bosons",
  //    5);
  // static SwitchOption interfaceDecayChannelZZ_ll
  //   (interfaceDecayChannel,
  //    "H -> Z Z -> l lbar",
  //    "Higgs decay into Z bosons into lepton antilepton",
  //    6);
  // static SwitchOption interfaceDecayChannelZZ_lnu
  //   (interfaceDecayChannel,
  //    "H -> Z Z -> l nu",
  //    "Higgs decay into Z bosons into lepton neutrino",
  //    7);

  // static Switch<VBFNLOMEPP2hJetJetJet,bool> interfaceNarrowWidth
  //   ("NarrowWidth",
  //    "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
  //    &VBFNLOMEPP2hJetJetJet::theNarrowWidth, true, true, false);
  // static SwitchOption interfaceNarrowWidthTrue
  //   (interfaceNarrowWidth,
  //    "True",
  //    "Calculate with narrow width approximation",
  //    true);
  // static SwitchOption interfaceNarrowWidthFalse
  //   (interfaceNarrowWidth,
  //    "False",
  //    "No narrow width approximation",
  //    false);

}
