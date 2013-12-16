
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEPP2WJetJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOMEPP2WJetJet::VBFNLOMEPP2WJetJet()
  : VBFNLOMEVVJJChargedBase(), theNarrowWidth(true) {
}

VBFNLOMEPP2WJetJet::~VBFNLOMEPP2WJetJet(){}

int VBFNLOMEPP2WJetJet::Njets() const {return 2;}

IBPtr VBFNLOMEPP2WJetJet::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEPP2WJetJet::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOMEPP2WJetJet::prepareMomenta(double (& pbar)[14][4],
					double (&qbar) [5],int) const{

  L5MomToDouble( (mePartonData()[0]->id() > 0 ? 
  		  lastMEMomenta()[0]  : 
  		  lastMEMomenta()[2]), &pbar[0][0]);

  L5MomToDouble( (mePartonData()[2]->id() > 0 ? 
  		  lastMEMomenta()[2] : 
  		  lastMEMomenta()[0]), &pbar[1][0]);

  L5MomToDouble( (mePartonData()[1]->id() > 0 ? 
  		  lastMEMomenta()[1] : 
  		  lastMEMomenta()[3]), &pbar[2][0]);

  L5MomToDouble( (mePartonData()[3]->id() > 0 ? 
  		  lastMEMomenta()[3] : 
  		  lastMEMomenta()[1]), &pbar[3][0]);

  L5MomToDouble( (lastMEMomenta()[4]), &pbar[4][0]);

  L5MomToDouble( (lastMEMomenta()[5]), &pbar[5][0]);

    
  return;
}

void VBFNLOMEPP2WJetJet::VbfnloMe2(double pbar[14][4],int mePartonSign[6],double [5],int,int nlo,int bosdec,double& uuX,double& ddX,double& Xcc,double& Xss) const{

  initProcess(nlo);

  if (theproducedBoson == getParticleData(ParticleID::Wplus))
    QQWPQQ(pbar,mePartonSign,nlo,
	   uuX,ddX,Xcc,Xss);
  else if (theproducedBoson == getParticleData(ParticleID::Wminus))
    QQWMQQ(pbar,mePartonSign,nlo,
	   uuX,ddX,Xcc,Xss);

  // cout << "nlo=" << nlo << endl;
  // cout << uucc<< "   " <<uuss<< "   " <<ddcc<< "   " <<ddss<< "   " <<udsc<< "   " <<ducs << endl;
  //  for (int i=0; i<nDim(); i++){
  //    //nDim()-1 because one random number is needed for Frn
  //    cout << "R("<<i+1<<")="<<theRandomNumbers[i]<<"d0"<<endl;;
  //  }
  //average over spins
  double fac = 4*9;



  //give result in units of sHat
  fac/=pow(lastSHat()/GeV2,2);




  uuX/=fac;
  ddX/=fac;
  Xcc/=fac;
  Xss/=fac;
  
  return;
}

void VBFNLOMEPP2WJetJet::initProcess(const int& nlo) const{

  CGLOBALI.SIGN1 = 1;
  CGLOBALI.SIGN2 = 1;
  CGLOBALI.N_V = 2;
  CGLOBALI.N_P = 4;
  CGLOBALI.PS_DIMENSION = 7;
  CGLOBALI.NLO_LOOPS = 1;


  //  double sca=renormalizationScale()/GeV2;
  // if (nlo == 0) { 
  //   CSCALES.ALS[0][0] = ALPHAS5_HARDWIRED(sca,0);
  //   CSCALES.ALS[0][1] = ALPHAS5_HARDWIRED(sca,0);
  // }
  // else {
  //   CSCALES.ALS[0][0] = ALPHAS5_HARDWIRED(sca,1);
  //   CSCALES.ALS[0][1] = ALPHAS5_HARDWIRED(sca,1);
  // }

  CSCALES.ALS[0][0] = lastAlphaS();
  CSCALES.ALS[0][1] = lastAlphaS();
  
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


  if (abs(nlo) == 1) {
    double p[5][4];
    double v[9][4];
    L5MomToDouble( (lastMEMomenta()[0]), &p[0][0]);
    L5MomToDouble( (lastMEMomenta()[1]), &p[1][0]);
    L5MomToDouble( (lastMEMomenta()[2]), &p[2][0]);
    L5MomToDouble( (lastMEMomenta()[3]), &p[3][0]);

    L5MomToDouble( (lastMEMomenta()[4]), &v[1][0]);
    L5MomToDouble( (lastMEMomenta()[5]), &v[0][0]);

    DEFINE_LEPTENS(2,0,1,p,v);
  }
    
}


// Energy2 VBFNLOMEPP2WJetJet::factorizationScale() const {
//   if ( theUserScale != ZERO )
//     return sqr(theUserScale);
//   return lastSHat();
// }

// Energy2 VBFNLOMEPP2WJetJet::renormalizationScale() const {
//   if ( theUserScale != ZERO )
//     return sqr(theUserScale);
//   return lastSHat();
// }


double VBFNLOMEPP2WJetJet::colourCorrelatedME2(pair<int,int> ij) const {
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
    lastME2(0.0);
    return lastME2();
  }
  return -me2();
}

double VBFNLOMEPP2WJetJet::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {
  generator()->logWarning(Exception() 
  			  << "A non-exisiting spin correlation was requested "
  			  << "from the matrix element '" << name() << "'."
  			  << Exception::warning);
  lastME2(0.0);
  return lastME2();
}

bool VBFNLOMEPP2WJetJet::noDipole(int emitter, int spectator) const{
  if (emitter == 0 && spectator == 2) return false;
  if (emitter == 2 && spectator == 0) return false;
  if (emitter == 1 && spectator == 3) return false;
  if (emitter == 3 && spectator == 1) return false;
  return true;
}


void VBFNLOMEPP2WJetJet::getDiagrams() const {

  if(theQuarkFlavours.size() == 0) throw ThePEG::Exception() << "Please insert pointers to quark flavours "
				<< "in the ME interface."
				<< ThePEG::Exception::abortnow;

  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);


  PDVector QuarksAndAntiQuarks1,QuarksAndAntiQuarks2;
  for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	q != theQuarkFlavours.end(); ++q ){
    if (requestedAsIncoming1(*q)) QuarksAndAntiQuarks1.push_back(*q);
    if (requestedAsIncoming1((**q).CC())) QuarksAndAntiQuarks1.push_back((**q).CC());
    if (requestedAsIncoming2(*q)) QuarksAndAntiQuarks2.push_back(*q);
    if (requestedAsIncoming2((**q).CC())) QuarksAndAntiQuarks2.push_back((**q).CC());
  }


  // get the decay products
  tcPDPtr dec1, dec2;

  if (theDecayLeptons->id() > 0){
    dec2 = theDecayLeptons->CC();
    dec1 = theDecayLeptons;
  }
  else {
    dec2 = theDecayLeptons;
    dec1 = theDecayLeptons->CC();
  }
  
  if (theExchangedBoson == Z0) {
    for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	  q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
      for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
	    q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,   1, *q1,  3, *q2, 2, Z0,  7,dec1,  7, dec2, -1)));
      }
    }
  }
  
  else if ( theExchangedBoson == Wplus || theExchangedBoson == Wminus) {
    for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	  q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
      for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
	    q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	Charge chargediff = (**q1).charge() - (**q2).charge() ;
	if ( chargediff != ZERO ) {
	  tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	  tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	  if ((**q1).charge()-(*q1prime).charge() == (*q2prime).charge()-(**q2).charge()){
	    if ( chargediff > ZERO) {
	      add(new_ptr((Tree2toNDiagram(4),*q1,Wplus,Wplus,*q2,   1,q1prime,  3,q2prime,  2,Z0,  7,dec1,  7, dec2, -2)));
	    }
	    else if ( chargediff < ZERO ) {
	      add(new_ptr((Tree2toNDiagram(4),*q1,Wminus,Wminus,*q2,   1,q1prime,  3,q2prime,  2,Z0,  7,dec1,  7, dec2, -2)));
	    }
	  }
	}
      }
    }
  }
  else {
    throw ThePEG::Exception() << "Please insert a pointer to W or Z "
			      << "boson as ExchangedBoson "
			      << "in the ME interface."
			      << ThePEG::Exception::abortnow;
  } 
}

Selector<MEBase::DiagramIndex> 
VBFNLOMEPP2WJetJet::diagrams(const DiagramVector & ) const {
  Selector<MEBase::DiagramIndex> sel;
  sel.insert(1,0);
  return sel;
}

Selector<const ColourLines *>
VBFNLOMEPP2WJetJet::colourGeometries(tcDiagPtr) const {

  static const ColourLines cUL      (" 1  5,  4  6");
  static const ColourLines cULbar   (" 1  5, -4 -6");
  static const ColourLines cUbarL   ("-1 -5,  4  6");
  static const ColourLines cUbarLbar("-1 -5, -4 -6");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0)
    sel.insert(1.0, &cUL);
  else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0)
    sel.insert(1.0, &cULbar);
  else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() > 0)
    sel.insert(1.0, &cUbarL);
  else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0)
    sel.insert(1.0, &cUbarLbar);
  return sel; 
}


void VBFNLOMEPP2WJetJet::doinit(){
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


void VBFNLOMEPP2WJetJet::initPSGen() {


  int bos=2;
  
  BLIPSIVNJ.RM2=BKOPOU.XM2[bos-1];
  BLIPSIVNJ.RMG=BKOPOU.XMG[bos-1];
  
  BLIPSIVNJ.RM2MIN=sqr(15.);
  BLIPSIVNJ.RM2MAX = sqr(5000.);
  
  BLIPSIVNJ.S=lastS()/GeV2;
  BLIPSIVNJ.M2MIN=0.01;
  
  //these should probably be taken from some cut class eventually
  BLIPSIVNJ.YJMIN[0]=0;
  BLIPSIVNJ.YJMIN[1]=0;
  
  BLIPSIVNJ.YJMAX[0]=5;
  BLIPSIVNJ.YJMAX[1]=5;
  
  BLIPSIVNJ.PTJMIN[0]=20;
  BLIPSIVNJ.PTJMIN[1]=20;
  
  BLIPSIVNJ.EJMIN[0]=0;
  BLIPSIVNJ.EJMIN[1]=0;
  
  BLIPSIVNJ.INFOJ[0]=-1;
  BLIPSIVNJ.INFOJ[1]=-1;
  

  return;
}



bool VBFNLOMEPP2WJetJet::generateKinematics(const double * r){

  if ( phasespace() ) {
    initPSGen();
    return MatchboxMEBase::generateKinematics(r);
  }

  int Fn = Njets();
  double Frd[100], Frn, Fk1[4],Fk2[4],Fq[5],Fd[10][4],Fd1[4],Fd2[4],Fp[2][4],Fx1,Fx2,Fw;
  int Fnw = 0;

  if (theNarrowWidth) Fnw =1;

  for (int i=0; i<nDim()-1; i++){
    //nDim()-1 because one random number is needed for Frn
    Frd[i]=r[i];
    theRandomNumbers[i]=r[i];
  }
  theRandomNumbers[nDim()-1] = r[nDim()-1];
  Frn = r[nDim()-1];
  
  initPSGen();

  

  LIPSN(Fn,Frd,Frn,Fk1,Fk2,Fq,Fd1,Fd2,Fp,Fx1,Fx2,Fw,Fnw);

  if(Fw > 0){
    lastMEMomenta()[0]=(Lorentz5Momentum(Fk1[1]*GeV,Fk1[2]*GeV,Fk1[3]*GeV,Fk1[0]*GeV));
    lastMEMomenta()[1]=(Lorentz5Momentum(Fk2[1]*GeV,Fk2[2]*GeV,Fk2[3]*GeV,Fk2[0]*GeV));
    for (int i=0; i<Fn; i++){
      lastMEMomenta()[2+i]=(Lorentz5Momentum(Fp[i][1]*GeV,Fp[i][2]*GeV,Fp[i][3]*GeV,Fp[i][0]*GeV));
    }

    lastMEMomenta()[5]=(Lorentz5Momentum(Fd1[1]*GeV,Fd1[2]*GeV,Fd1[3]*GeV,Fd1[0]*GeV));  
    lastMEMomenta()[4]=(Lorentz5Momentum(Fd2[1]*GeV,Fd2[2]*GeV,Fd2[3]*GeV,Fd2[0]*GeV));  
    
    // cout << "NDecayProducts()= " << NDecayProducts() << endl;
      // for (int i=0; i<Fn+4; i++){
      // 	cout<<"lastMEMomenta()["<<i<<"]="<<setprecision(32)<<lastMEMomenta()[i]/GeV<<endl;
      //  }
    //     cout<<"lastMEMomenta()[4+5]=" << (lastMEMomenta()[4]+lastMEMomenta()[5])/GeV << endl;
     //division by sHat takes place in LIPSN, but this is ok for 
    //hjj (2 to 3 process)
    //for processes with h decay, we need to divide by some powers of shat
    
    Energy2 thisSHat=(meMomenta()[0] + meMomenta()[1]).m2();
    Fw/=thisSHat/GeV2;
    
    //revert the multiplication with (hbar*c)^2/2   
    //and divide by Feynman x-es
    jacobian(Fw*2/3.8937E11/Fx1/Fx2);
    setScale();
    logGenerateKinematics(r);
    return true; 
  } else { 
    jacobian(0.0);
    return false;
  }
}

int VBFNLOMEPP2WJetJet::nDim() const {
  if (phasespace()) return phasespace()->nDim(4);
  return 10;
}

ClassDescription<VBFNLOMEPP2WJetJet> VBFNLOMEPP2WJetJet::initVBFNLOMEPP2WJetJet;
// Definition of the static class description member.



void VBFNLOMEPP2WJetJet::persistentOutput(PersistentOStream & os) const {
  os  << theExchangedBoson
      << theDecayLeptons << theNarrowWidth;
}

void VBFNLOMEPP2WJetJet::persistentInput(PersistentIStream & is, int) {
  is  >> theExchangedBoson
      >> theDecayLeptons >> theNarrowWidth;
}


void VBFNLOMEPP2WJetJet::Init() {

  static ClassDocumentation<VBFNLOMEPP2WJetJet> documentation
    ("VBFNLOMEPP2WJetJet");

  static Reference<VBFNLOMEPP2WJetJet,ParticleData> interfaceDecayLeptons
    ("DecayLeptons",
     "The final state leptons for this process.",
     &VBFNLOMEPP2WJetJet::theDecayLeptons, false, false, true, true, false);

  static Switch<VBFNLOMEPP2WJetJet,bool> interfaceNarrowWidth
    ("NarrowWidth",
     "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
     &VBFNLOMEPP2WJetJet::theNarrowWidth, true, true, false);
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
