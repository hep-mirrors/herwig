
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEPP2ZJetJetJet.h"
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


VBFNLOMEPP2ZJetJetJet::VBFNLOMEPP2ZJetJetJet()
  : VBFNLOMEVVJJNeutralBase(), theNarrowWidth(true) {
}

VBFNLOMEPP2ZJetJetJet::~VBFNLOMEPP2ZJetJetJet(){}

int VBFNLOMEPP2ZJetJetJet::Njets() const {return 3;}



void VBFNLOMEPP2ZJetJetJet::VbfnloMe2(double pbar[14][4],int mePartonSign[6],double qbar[5],int meGluonSign,int nlo,int,double& uucc,double& uuss,double& ddcc,double& ddss,double& udsc,double& ducs) const{

  initProcess(nlo);

  QQZQQJ(pbar,mePartonSign,qbar,meGluonSign,
  	uucc,uuss,ddcc,ddss,udsc,ducs);

  double fac = 1;
  //average over spins and colors
  if (meGluonSign == 1) {fac*=4*9;}  //outgoing gluon
  else fac*=4*24;

  //give result in units of sHat
  fac/=pow(lastSHat()/GeV2,3);

  uucc/=fac;
  uuss/=fac;
  ddcc/=fac;
  ddss/=fac;
  udsc/=fac;
  ducs/=fac;
  
  return;
}

void VBFNLOMEPP2ZJetJetJet::initDiagramContainers() const {}

void VBFNLOMEPP2ZJetJetJet::initProcess(const int & nlo) const{

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
  CGLOBALI.N_V = 2;
  CGLOBALI.NLO_LOOPS = 0;

  int lflavr[2];
  if (abs(theDecayLeptons->id()) == 11
      || abs(theDecayLeptons->id()) == 13 
      || abs(theDecayLeptons->id()) == 15) {
    lflavr[0] = 2;
    lflavr[1] = 2;
  }
  else if (abs(theDecayLeptons->id()) == 12
	   || abs(theDecayLeptons->id()) == 14 
	   || abs(theDecayLeptons->id()) == 16) {
    lflavr[0] = 1;
    lflavr[1] = 1;
  }
  else throw ThePEG::Exception() << "The given decay products are not supported!"
				 << ThePEG::Exception::abortnow;


  QQBQQI(2,lflavr);
  CSCALES.ALS[0][0] = lastAlphaS();
  CSCALES.ALS[0][1] = lastAlphaS();
}

double VBFNLOMEPP2ZJetJetJet::colourCorrelatedME2(pair<int,int>) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

double VBFNLOMEPP2ZJetJetJet::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

void VBFNLOMEPP2ZJetJetJet::doinit(){
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


void VBFNLOMEPP2ZJetJetJet::initPSGen() {
 
  int bos=2;
  
  BLIPSIVNJ.RM2=BKOPOU.XM2[bos-1];
  BLIPSIVNJ.RMG=BKOPOU.XMG[bos-1];
  
  BLIPSIVNJ.RM2MIN=sqr(15.);
  BLIPSIVNJ.RM2MAX = sqr(5000.);
  
  BLIPSIVNJ.S=lastS()/GeV2;
  BLIPSIVNJ.M2MIN=0.001;
  
  //these should probably be taken from some cut class eventually
  for (int i = 0; i < 4; i++){
    BLIPSIVNJ.YJMIN[i]=0;
    BLIPSIVNJ.YJMAX[i]=5;
    BLIPSIVNJ.PTJMIN[i]=20;
    BLIPSIVNJ.EJMIN[i]=0;
    BLIPSIVNJ.INFOJ[i]=-1;
  }

  return;
}



bool VBFNLOMEPP2ZJetJetJet::generateKinematics(const double * r){

  if ( phasespace() ) {
    initPSGen();
    return MatchboxMEBase::generateKinematics(r);
  }

  int Fn = Njets();
  double Frd[100], Frn, Fk1[4],Fk2[4],Fq[5],Fd1[4],Fd2[4],Fp[4][4],Fx1,Fx2,Fw;
  int Fnw = 0;
  if (NarrowWidth()) Fnw =1;


  for (int i=0; i<nDim()-1; i++){
    Frd[i]=r[i];
  }
  Frn= r[nDim()-1];
  
  initPSGen();  
  LIPSN(Fn,Frd,Frn,Fk1,Fk2,Fq,Fd1,Fd2,Fp,Fx1,Fx2,Fw,Fnw);
  
  if(Fw > 0 ){

    lastMEMomenta()[0]=(Lorentz5Momentum(Fk1[1]*GeV,Fk1[2]*GeV,Fk1[3]*GeV,Fk1[0]*GeV));
    lastMEMomenta()[1]=(Lorentz5Momentum(Fk2[1]*GeV,Fk2[2]*GeV,Fk2[3]*GeV,Fk2[0]*GeV));
    for (int i=0; i<Fn; i++){
      lastMEMomenta()[2+i]=(Lorentz5Momentum(Fp[i][1]*GeV,Fp[i][2]*GeV,Fp[i][3]*GeV,Fp[i][0]*GeV));
    }   
    
    lastMEMomenta()[6]=(Lorentz5Momentum(Fd1[1]*GeV,Fd1[2]*GeV,Fd1[3]*GeV,Fd1[0]*GeV));  
    lastMEMomenta()[5]=(Lorentz5Momentum(Fd2[1]*GeV,Fd2[2]*GeV,Fd2[3]*GeV,Fd2[0]*GeV));  
    
    //revert the multiplication with (hbar*c)^2/2   
    //division by sHat also takes place in LIPSN, need to divide by 
    //one more power of SHat
    //divide by Feynman x-es
    Energy2 thisSHat=(meMomenta()[0] + meMomenta()[1]).m2();
    Fw/=sqr(thisSHat/GeV2);

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

bool VBFNLOMEPP2ZJetJetJet::noDipole(int emitter, int emission, int spectator) const{
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

int VBFNLOMEPP2ZJetJetJet::nDim() const {
  if (phasespace()) return phasespace()->nDim(5);
  return 13;
}

AbstractClassDescription<VBFNLOMEPP2ZJetJetJet> VBFNLOMEPP2ZJetJetJet::initVBFNLOMEPP2ZJetJetJet;
// Definition of the static class description member.

void VBFNLOMEPP2ZJetJetJet::persistentOutput(PersistentOStream & os) const {
  os << theDecayLeptons
     << theNarrowWidth;
}

void VBFNLOMEPP2ZJetJetJet::persistentInput(PersistentIStream & is, int) {
  is  >> theDecayLeptons 
      >> theNarrowWidth;
}

void VBFNLOMEPP2ZJetJetJet::Init() {

  static ClassDocumentation<VBFNLOMEPP2ZJetJetJet> documentation
    ("VBFNLOMEPP2ZJetJetJet");

  static Reference<VBFNLOMEPP2ZJetJetJet,ParticleData> interfaceDecayLeptons
    ("DecayLeptons",
     "The final state leptons for this process.",
     &VBFNLOMEPP2ZJetJetJet::theDecayLeptons, false, false, true, true, false);

  static Switch<VBFNLOMEPP2ZJetJetJet,bool> interfaceNarrowWidth
    ("NarrowWidth",
     "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
     &VBFNLOMEPP2ZJetJetJet::theNarrowWidth, true, true, false);
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
