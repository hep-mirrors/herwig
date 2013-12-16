#include "VBFNLOMEVVJJChargedBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "VBFNLOCommonBlocks.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/PDT/ParticleData.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;

VBFNLOMEVVJJChargedBase::VBFNLOMEVVJJChargedBase()
  : VBFNLOMEBase(), theIncoming1(true), theIncoming2(true) {}

VBFNLOMEVVJJChargedBase::~VBFNLOMEVVJJChargedBase(){}

void VBFNLOMEVVJJChargedBase::doinit() {}

double VBFNLOMEVVJJChargedBase::me2() const{

  // Energy2 M2 = sqr(getParticleData(ParticleID::Z0)->mass());
  // Energy2 G2 = sqr(getParticleData(ParticleID::Z0)->width());
  // Energy2 MH2 = sqr(getParticleData(ParticleID::h0)->mass());

  // vector<Lorentz5Momentum> meMomentaCM = meMomenta();

  // Boost tocm = (meMomentaCM[0]+meMomentaCM[1]).findBoostToCM();
  // if ( tocm.mag2() > Constants::epsilon ) {
  //   for ( vector<Lorentz5Momentum>::iterator p = meMomentaCM.begin();
  // 	  p != meMomentaCM.end(); ++p ) {
  //     p->boost(tocm);
  //   }
  // }

  // /*
  // if ( mePartonData().size() == 6 )
  //   cerr << "\n\n--------------------------------------------------------------------------------\n";
  // cerr << fixed << setprecision(8);
  // cerr << name() << " used in " << lastXComb().matrixElement()->name()
  //      << " check momenta ...\n";
  // Lorentz5Momentum sum = meMomenta()[0] + meMomenta()[1];
  // vector<Lorentz5Momentum>::const_iterator p = meMomenta().begin();
  // cPDVector::const_iterator pd = mePartonData().begin();
  // for ( ; p != meMomenta().end(); ++p, ++pd ) {
  //   cerr << (**pd).PDGName() << " : "
  // 	 << (*p/GeV) << " " << p->m2()/GeV2 << " " << sqrt(p->m2()/GeV2) << "\n";
  //   if ( p > meMomenta().begin() + 1 )
  //     sum -= *p;
  // }
  // cerr << "conserved? " << sum/GeV << "\n"
  //      << flush;
  // */

  // if ( mePartonData().size() == 5 ) {

  //   Lorentz5Momentum pq1,pq2,pqbar1,pqbar2;
  //   if ( mePartonData()[0]->id() < 0 ) {
  //     pq1 = -meMomentaCM[0];
  //     pqbar1 = meMomentaCM[2];
  //   } else {
  //     pq1 = meMomentaCM[2];
  //     pqbar1 = -meMomentaCM[0];
  //   }
  //   if ( mePartonData()[1]->id() < 0 ) {
  //     pq2 = -meMomentaCM[1];
  //     pqbar2 = meMomentaCM[3];
  //   } else {
  //     pq2 = meMomentaCM[3];
  //     pqbar2 = -meMomentaCM[1];
  //   }

  //   Energy2 Q2Upper = (pq1+pqbar1).m2();
  //   Energy2 Q2Lower = (pq2+pqbar2).m2();

  //   Energy4 qcd =
  //     (3.*4.*(pq1*pqbar1))*(3.*4.*(pq2*pqbar2));

  //   return
  //     (1./(4.*3.*3.))*(lastSHat()*MH2*qcd)/((sqr(Q2Upper-M2)+M2*G2)*(sqr(Q2Lower-M2)+M2*G2));

  // }
  // if ( mePartonData().size() == 6 &&
  //      mePartonData()[0]->id() != ParticleID::g &&
  //      mePartonData()[1]->id() != ParticleID::g ) {

  //   Lorentz5Momentum pq1,pq2,pqbar1,pqbar2,pg;
  //   if ( mePartonData()[0]->id() < 0 ) {
  //     pq1 = -meMomentaCM[0];
  //     pqbar1 = meMomentaCM[2];
  //   } else {
  //     pq1 = meMomentaCM[2];
  //     pqbar1 = -meMomentaCM[0];
  //   }
  //   if ( mePartonData()[1]->id() < 0 ) {
  //     pq2 = -meMomentaCM[1];
  //     pqbar2 = meMomentaCM[3];
  //   } else {
  //     pq2 = meMomentaCM[3];
  //     pqbar2 = -meMomentaCM[1];
  //   }
  //   pg = meMomentaCM[4];

  //   // upper line emitting
  //   Energy2 Q2Upper = (pq1+pqbar1+pg).m2();
  //   Energy2 Q2Lower = (pq2+pqbar2).m2();
  //   Energy4 qcd =
  //     (4.*16.*Constants::pi*lastSHat()*lastAlphaS())*
  //     ((sqr(pq1*pg)+sqr(pqbar1*pg)+2.*(pq1*pqbar1+pq1*pg)*(pq1*pqbar1+pqbar1*pg))/((pq1*pg)*(pqbar1*pg)))*
  //     (3.*4.*(pq2*pqbar2));
  //   double res =
  //     (lastSHat()*MH2*qcd)/((sqr(Q2Upper-M2)+M2*G2)*(sqr(Q2Lower-M2)+M2*G2));

  //   // lower line emitting
  //   Q2Upper = (pq1+pqbar1).m2();
  //   Q2Lower = (pq2+pqbar2+pg).m2();
  //   qcd =
  //     (4.*16.*Constants::pi*lastSHat()*lastAlphaS())*
  //     ((sqr(pq2*pg)+sqr(pqbar2*pg)+2.*(pq2*pqbar2+pq2*pg)*(pq2*pqbar2+pqbar2*pg))/((pq2*pg)*(pqbar2*pg)))*
  //     (3.*4.*(pq1*pqbar1));
  //   res +=
  //     (lastSHat()*MH2*qcd)/((sqr(Q2Upper-M2)+M2*G2)*(sqr(Q2Lower-M2)+M2*G2));

  //   return (1./(4.*3.*3.))*res;

  // }

  // assert(false);

  int mePartonSign[6]={1,1,1,1,1,1};

  int gluonIndex = 100;
  int meGluonSign = 1;

  for (int i = 0; i < mePartonData().size(); i++){
    if (mePartonData()[i]->id() == 21) gluonIndex = i;
  }
  
  //these are indices that make sure we get the proper subprocess
  //when we have gluons in the initial state for the real emission
  //matrix elements

  int U = 0;
  int L = 1;

  if (gluonIndex == 0) U = 4;
  if (gluonIndex == 1) L = 4;
  if (gluonIndex < 2) meGluonSign = -1;

  if (mePartonData()[U]->id() < 0) mePartonSign[0]=-1;
  if (mePartonData()[2]->id() < 0) mePartonSign[1]=-1;
  if (mePartonData()[L]->id() < 0) mePartonSign[2]=-1;
  if (mePartonData()[3]->id() < 0) mePartonSign[3]=-1;
  //  cout << "mePartonData()[4]->id()="<<mePartonData()[4]->id()<<"  mePartonData()[5]->id()="<<mePartonData()[5]->id()<<endl;
  if (mePartonData()[mePartonData().size()-1]->id() < 0) mePartonSign[4]=-1; //changed mePartonData indices here
  if (mePartonData()[mePartonData().size()-2]->id() < 0) mePartonSign[5]=-1; //because VBFNLO uses a strange convention

  SPLITCB.ALLSUBPROCS = false; 
  SPLITCB.SUBPROCID = getSubprocessID(gluonIndex);

  if (mePartonData()[0]->id() == 21) SPLITCB.GLUONID = 0;
  else if (mePartonData()[1]->id() == 21) SPLITCB.GLUONID = 1;
  else SPLITCB.GLUONID = -1;


  int nlo=0;
  int bosdec=0;

  //  double uucc, uuss, ddcc, ddss, udsc, ducs;
  double uuX, ddX, Xcc, Xss;

  double pbar[14][4];
  double qbar[5];
  prepareMomenta(pbar,qbar,gluonIndex);

  VbfnloMe2(pbar,mePartonSign,qbar,meGluonSign,nlo,bosdec,
	uuX,ddX,Xcc,Xss);
  
  double result[6];
  
  result[0]=uuX;
  result[1]=ddX;
  result[2]=Xcc;
  result[3]=Xss;

  lastME2(result[SPLITCB.SUBPROCID-1]);
  // cout << setprecision(32);
  // cout << "lastME2()=" << lastME2() << endl;
  return lastME2();
  
}

int VBFNLOMEVVJJChargedBase::getSubprocessID(int gluonIndex) const{

  //these are indices that make sure we get the proper subprocess
  //when we have gluons in the initial state for the real emission
  //matrix elements
  int U = 0;
  int L = 1;

  // cout<<mePartonData()[0]->id()<<"   "
  //     <<mePartonData()[1]->id()<<"   "
  //     <<mePartonData()[2]->id()<<"   "
  //     <<mePartonData()[3]->id()<<"   "
  //     <<mePartonData()[4]->id()<<"   ";

  if (gluonIndex == 0) U = 4;
  if (gluonIndex == 1) L = 4;

  if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
       SU2Helper::isSU2Up(mePartonData()[2]) 
       ){
    return 1;
  }
  else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
	    SU2Helper::isSU2Down(mePartonData()[2])  
	    ){
    return 2;
  }
  else if ( SU2Helper::isSU2Up(mePartonData()[L]) &&
	    SU2Helper::isSU2Up(mePartonData()[3]) 
	    ){
    return 3;
  }
  else if ( SU2Helper::isSU2Down(mePartonData()[L]) &&
	    SU2Helper::isSU2Down(mePartonData()[3]) 
	    ){
    return 4;
  }

  // if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
  //      SU2Helper::isSU2Up(mePartonData()[2]) &&
  //      SU2Helper::isSU2Up(mePartonData()[L]) &&
  //      SU2Helper::isSU2Up(mePartonData()[3]) 
  //      ){
  //   return 1;
  // }
  // else if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
  // 	    SU2Helper::isSU2Up(mePartonData()[2]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[3]) 
  // 	    ){
  //   return 2;
  // }
  // else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[2]) &&
  // 	    SU2Helper::isSU2Up(mePartonData()[L]) &&
  // 	    SU2Helper::isSU2Up(mePartonData()[3]) 
  // 	    ){
  //   return 3;
  // }
  // else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[2]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[3]) 
  // 	    ){
  //   return 4;
  // }
  // else {
  //   if (gluonIndex == 0) {
  //     if (mePartonData()[L]->id()>0){
  // 	if (SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 	    SU2Helper::isSU2Up(mePartonData()[3])) {
  // 	  return 5;
  // 	}
  // 	else if (SU2Helper::isSU2Up(mePartonData()[L]) &&
  // 		 SU2Helper::isSU2Down(mePartonData()[3])) {
  // 	  return 6;
  // 	}
  //     }
  //     else if (mePartonData()[L]->id()<0) {
  // 	if (SU2Helper::isSU2Up(mePartonData()[L]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[3])) {
  // 	  return 5;
  // 	}
  // 	else if (SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 		 SU2Helper::isSU2Up(mePartonData()[3])) {
  // 	  return 6;
  // 	}
  //     }
  //   }
  //   else if (gluonIndex == 1) {
  //     if (mePartonData()[U]->id()>0){
  // 	if (SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 	    SU2Helper::isSU2Up(mePartonData()[2])) {
  // 	  return 5;
  // 	}
  // 	else if (SU2Helper::isSU2Up(mePartonData()[U]) &&
  // 		 SU2Helper::isSU2Down(mePartonData()[2])) {
  // 	  return 6;
  // 	}
  //     }
  //     else if (mePartonData()[U]->id()<0) {
  // 	if (SU2Helper::isSU2Up(mePartonData()[U]) &&
  // 	    SU2Helper::isSU2Down(mePartonData()[2])) {
  // 	  return 5;
  // 	}
  // 	else if (SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 		 SU2Helper::isSU2Up(mePartonData()[2])) {
  // 	  return 6;
  // 	}
  //     }
  //   }
  //   else {	    
  //     if(mePartonData()[U]->id()*mePartonData()[L]->id() > 0) {
  // 	if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
  // 	     SU2Helper::isSU2Down(mePartonData()[2]) &&
  // 	     SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 	     SU2Helper::isSU2Up(mePartonData()[3]) 
  // 	     ){
  // 	  return 5;
  // 	}
  // 	else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 		  SU2Helper::isSU2Up(mePartonData()[2]) &&
  // 		  SU2Helper::isSU2Up(mePartonData()[L]) &&
  // 		  SU2Helper::isSU2Down(mePartonData()[3]) 
  // 		  ){
  // 	  return 6;
  // 	}
  //     }

  //     else if (mePartonData()[U]->id()<0){
  // 	if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 	     SU2Helper::isSU2Up(mePartonData()[2]) &&
  // 	     SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 	     SU2Helper::isSU2Up(mePartonData()[3]) 
  // 	     ){
  // 	  return 5;
  // 	}
  // 	else if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
  //   		  SU2Helper::isSU2Down(mePartonData()[2]) &&
  // 		  SU2Helper::isSU2Up(mePartonData()[L]) &&
  //   		  SU2Helper::isSU2Down(mePartonData()[3]) 
  // 		  ){
  // 	  return 6;
  // 	}
  //     }
  //     else if (mePartonData()[L]->id()<0){
  // 	if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
  // 	     SU2Helper::isSU2Down(mePartonData()[2]) &&
  // 	     SU2Helper::isSU2Up(mePartonData()[L]) &&
  // 	     SU2Helper::isSU2Down(mePartonData()[3]) 
  // 	     ){
  // 	  return 5;
  // 	}
  // 	else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
  // 		  SU2Helper::isSU2Up(mePartonData()[2]) &&
  // 		  SU2Helper::isSU2Down(mePartonData()[L]) &&
  // 		  SU2Helper::isSU2Up(mePartonData()[3]) 
  // 		  ){
  // 	  return 6;
  // 	}
      
  //     }
  //   }
  // }
  throw ThePEG::Exception() << "Could not determine subprocess ID."
			    << " Got the following particle IDs:"
			    <<  mePartonData()[0]->id() << " "
			    <<  mePartonData()[1]->id() << " "
			    <<  mePartonData()[2]->id() << " "
			    <<  mePartonData()[3]->id() << " "
			    << ThePEG::Exception::abortnow;  
}

Energy2 VBFNLOMEVVJJChargedBase::factorizationScale() const {
  return 10000*GeV2;
  if ( theUserScale != ZERO )
    return sqr(theUserScale);
  return lastSHat();
}

Energy2 VBFNLOMEVVJJChargedBase::renormalizationScale() const {
  return 10000*GeV2;
  if ( theUserScale != ZERO )
    return sqr(theUserScale);
  return lastSHat();
}



bool VBFNLOMEVVJJChargedBase::requestedAsIncoming1(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming1() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming1()) return true;
  if (p->id() < 0 && !Incoming1()) return true;
  return false;
}

bool VBFNLOMEVVJJChargedBase::requestedAsIncoming2(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming2() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming2()) return true;
  if (p->id() < 0 && !Incoming2()) return true;
  return false;
}

bool VBFNLOMEVVJJChargedBase::requestedAsIncoming1tc(tcPDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming1() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming1()) return true;
  if (p->id() < 0 && !Incoming1()) return true;
  return false;
}

bool VBFNLOMEVVJJChargedBase::requestedAsIncoming2tc(tcPDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming2() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming2()) return true;
  if (p->id() < 0 && !Incoming2()) return true;
  return false;
}

AbstractClassDescription<VBFNLOMEVVJJChargedBase> VBFNLOMEVVJJChargedBase::initVBFNLOMEVVJJChargedBase;
// Definition of the static class description member.

void VBFNLOMEVVJJChargedBase::persistentOutput(PersistentOStream & os) const {
  os << theIncoming1 << theIncoming2 << theExchangedBoson;
}

void VBFNLOMEVVJJChargedBase::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming1 >> theIncoming2 >> theExchangedBoson;
}


void VBFNLOMEVVJJChargedBase::Init() {

  static ClassDocumentation<VBFNLOMEVVJJChargedBase> documentation
    ("VBFNLOMEVVJJChargedBase");

  // static RefVector<VBFNLOMEVVJJChargedBase,ParticleData> interfaceQuarkFlavours
  //   ("QuarkFlavours",
  //    "The quark flavours for this matrix element.",
  //    &VBFNLOMEVVJJChargedBase::theQuarkFlavours, -1, false, false, true, true, false);


  // static Parameter<VBFNLOMEVVJJChargedBase,Energy> interfaceUserScale
  //   ("UserScale",
  //    "A user defined renormalization scale.",
  //    &VBFNLOMEVVJJChargedBase::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
  //    false, false, Interface::lowerlim);


  static Switch<VBFNLOMEVVJJChargedBase,bool> interfaceIncoming1
    ("Incoming1",
     "Set to true/false for parton 1 being particle/antiparticle",
     &VBFNLOMEVVJJChargedBase::theIncoming1, true, true, false);
  static SwitchOption interfaceIncoming1Particle
    (interfaceIncoming1,
     "Particle",
     "Parton 1 is considered to be a particle",
     true);
  static SwitchOption interfaceIncoming1Antiparticle
    (interfaceIncoming1,
     "Antiparticle",
     "Parton 1 is considered to be an antiparticle",
     false);

  static Switch<VBFNLOMEVVJJChargedBase,bool> interfaceIncoming2
    ("Incoming2",
     "Set to true/false for parton 2 being particle/antiparticle",
     &VBFNLOMEVVJJChargedBase::theIncoming2, true, true, false);
  static SwitchOption interfaceIncoming2Particle
    (interfaceIncoming2,
     "Particle",
     "Parton 2 is considered to be a particle",
     true);
  static SwitchOption interfaceIncoming2Antiparticle
    (interfaceIncoming2,
     "Antiparticle",
     "Parton 2 is considered to be an antiparticle",
     false);
  
  static Reference<VBFNLOMEVVJJChargedBase,ParticleData> interfaceExchangedBoson
    ("ExchangedBoson",
     "The exchanged vector boson for this process (only W or Z).",
     &VBFNLOMEVVJJChargedBase::theExchangedBoson, false, false, true, true, false);

}
