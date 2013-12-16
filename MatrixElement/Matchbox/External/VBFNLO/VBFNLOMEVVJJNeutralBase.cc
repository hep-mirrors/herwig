#include "VBFNLOMEVVJJNeutralBase.h"
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

VBFNLOMEVVJJNeutralBase::VBFNLOMEVVJJNeutralBase()
  : VBFNLOMEBase(), theCurrent(neutral), theIncoming1(false),
    theIncoming2(true) {}

VBFNLOMEVVJJNeutralBase::VBFNLOMEVVJJNeutralBase(int current, bool incoming1, bool incoming2)
  : VBFNLOMEBase(), theCurrent(current), theIncoming1(incoming1),
    theIncoming2(incoming2) {}

VBFNLOMEVVJJNeutralBase::~VBFNLOMEVVJJNeutralBase(){}

// double VBFNLOMEVVJJNeutralBase::me2() const{
//   cerr << "calling VBFNLOMEVVJJNeutralBase::me2()\n" << flush;
//   int mePartonSign[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};

//   int gluonIndex = 100;
//   int meGluonSign = 1;

//   for (int i = 0; i < mePartonData().size(); i++){
//     if (mePartonData()[i]->id() == 21) gluonIndex = i;
//   }
  
//   //these are indices that make sure we get the proper subprocess
//   //when we have gluons in the initial state for the real emission
//   //matrix elements

//   int U = 0;
//   int L = 1;

//   if (gluonIndex == 0) U = 4;
//   if (gluonIndex == 1) L = 4;
//   if (gluonIndex < 2) meGluonSign = -1;

//   if (mePartonData()[U]->id() < 0) mePartonSign[0]=-1;
//   if (mePartonData()[2]->id() < 0) mePartonSign[1]=-1;
//   if (mePartonData()[L]->id() < 0) mePartonSign[2]=-1;
//   if (mePartonData()[3]->id() < 0) mePartonSign[3]=-1;

//   //changed mePartonData indices here
//   //because of VBFNLO convention
//   if (mePartonData()[mePartonData().size()-1]->id() < 0) mePartonSign[4]=-1; 
//   if (mePartonData()[mePartonData().size()-2]->id() < 0) mePartonSign[5]=-1; 

//   SPLITCB.ALLSUBPROCS = false; 
//   SPLITCB.SUBPROCID = getSubprocessID(gluonIndex);

//   if (mePartonData()[0]->id() == 21) SPLITCB.GLUONID = 0;
//   else if (mePartonData()[1]->id() == 21) SPLITCB.GLUONID = 1;
//   else SPLITCB.GLUONID = -1;

//   int nlo=0;
//   CGLOBALL.DOVIRTUALS = false;

//   int bosdec=0;

//   double uucc, uuss, ddcc, ddss, udsc, ducs;

//   double pbar[14][4];
//   double qbar[5];
//   prepareMomenta(pbar,qbar,gluonIndex);

//   VbfnloMe2(pbar,mePartonSign,qbar,meGluonSign,nlo,bosdec,
// 	uucc,uuss,ddcc,ddss,udsc,ducs);
  
//   double result[6];
  
//   result[0]=uucc;
//   result[1]=uuss;
//   result[2]=ddcc;
//   result[3]=ddss;
//   result[4]=udsc;
//   result[5]=ducs;

//   lastME2(result[SPLITCB.SUBPROCID-1]);
//   for (int i = 0; i < lastMEMomenta().size(); i++) {
//     //  cerr << "p[" << i << "]=" << lastMEMomenta()[i]/GeV << "\n" << flush;
//   }
//   //  cerr << "lastME2()= " << lastME2() << "\n" << flush;
  
//   return lastME2();
  
// }

int VBFNLOMEVVJJNeutralBase::getSubprocessID(int gluonIndex) const{

  //these are indices that make sure we get the proper subprocess
  //when we have gluons in the initial state for the real emission
  //matrix elements
  int U = 0;
  int L = 1;

  if (gluonIndex == 0) U = 4;
  if (gluonIndex == 1) L = 4;


  if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
       SU2Helper::isSU2Up(mePartonData()[2]) &&
       SU2Helper::isSU2Up(mePartonData()[L]) &&
       SU2Helper::isSU2Up(mePartonData()[3]) 
       ){
    return 1;
  }
  else if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
	    SU2Helper::isSU2Up(mePartonData()[2]) &&
	    SU2Helper::isSU2Down(mePartonData()[L]) &&
	    SU2Helper::isSU2Down(mePartonData()[3]) 
	    ){
    return 2;
  }
  else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
	    SU2Helper::isSU2Down(mePartonData()[2]) &&
	    SU2Helper::isSU2Up(mePartonData()[L]) &&
	    SU2Helper::isSU2Up(mePartonData()[3]) 
	    ){
    return 3;
  }
  else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
	    SU2Helper::isSU2Down(mePartonData()[2]) &&
	    SU2Helper::isSU2Down(mePartonData()[L]) &&
	    SU2Helper::isSU2Down(mePartonData()[3]) 
	    ){
    return 4;
  }
  else {
    if (gluonIndex == 0) {
      if (mePartonData()[L]->id()>0){
	if (SU2Helper::isSU2Down(mePartonData()[L]) &&
	    SU2Helper::isSU2Up(mePartonData()[3])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Up(mePartonData()[L]) &&
		 SU2Helper::isSU2Down(mePartonData()[3])) {
	  return 6;
	}
      }
      else if (mePartonData()[L]->id()<0) {
	if (SU2Helper::isSU2Up(mePartonData()[L]) &&
	    SU2Helper::isSU2Down(mePartonData()[3])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Down(mePartonData()[L]) &&
		 SU2Helper::isSU2Up(mePartonData()[3])) {
	  return 6;
	}
      }
    }
    else if (gluonIndex == 1) {
      if (mePartonData()[U]->id()>0){
	if (SU2Helper::isSU2Down(mePartonData()[U]) &&
	    SU2Helper::isSU2Up(mePartonData()[2])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Up(mePartonData()[U]) &&
		 SU2Helper::isSU2Down(mePartonData()[2])) {
	  return 6;
	}
      }
      else if (mePartonData()[U]->id()<0) {
	if (SU2Helper::isSU2Up(mePartonData()[U]) &&
	    SU2Helper::isSU2Down(mePartonData()[2])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Down(mePartonData()[U]) &&
		 SU2Helper::isSU2Up(mePartonData()[2])) {
	  return 6;
	}
      }
    }
    else {	    
      if(mePartonData()[U]->id()*mePartonData()[L]->id() > 0) {
	if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
	     SU2Helper::isSU2Down(mePartonData()[2]) &&
	     SU2Helper::isSU2Down(mePartonData()[L]) &&
	     SU2Helper::isSU2Up(mePartonData()[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
		  SU2Helper::isSU2Up(mePartonData()[2]) &&
		  SU2Helper::isSU2Up(mePartonData()[L]) &&
		  SU2Helper::isSU2Down(mePartonData()[3]) 
		  ){
	  return 6;
	}
      }

      else if (mePartonData()[U]->id()<0){
	if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
	     SU2Helper::isSU2Up(mePartonData()[2]) &&
	     SU2Helper::isSU2Down(mePartonData()[L]) &&
	     SU2Helper::isSU2Up(mePartonData()[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
    		  SU2Helper::isSU2Down(mePartonData()[2]) &&
		  SU2Helper::isSU2Up(mePartonData()[L]) &&
    		  SU2Helper::isSU2Down(mePartonData()[3]) 
		  ){
	  return 6;
	}
      }
      else if (mePartonData()[L]->id()<0){
	if ( SU2Helper::isSU2Up(mePartonData()[U]) &&
	     SU2Helper::isSU2Down(mePartonData()[2]) &&
	     SU2Helper::isSU2Up(mePartonData()[L]) &&
	     SU2Helper::isSU2Down(mePartonData()[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Down(mePartonData()[U]) &&
		  SU2Helper::isSU2Up(mePartonData()[2]) &&
		  SU2Helper::isSU2Down(mePartonData()[L]) &&
		  SU2Helper::isSU2Up(mePartonData()[3]) 
		  ){
	  return 6;
	}
      
      }
    }
  }
  throw ThePEG::Exception() << "Could not determine subprocess ID."
			    << " Got the following particle IDs:"
			    <<  mePartonData()[0]->id() << " "
			    <<  mePartonData()[1]->id() << " "
			    <<  mePartonData()[2]->id() << " "
			    <<  mePartonData()[3]->id() << " "
			    <<  mePartonData()[4]->id() << " "
			    << ThePEG::Exception::abortnow;  
}

// Energy2 VBFNLOMEVVJJNeutralBase::factorizationScale() const {
//   return 10000*GeV2;

//   // if ( theUserScale != ZERO )
//   //   return sqr(theUserScale);
//   // return lastSHat();
  
//   // Energy2 scale = 1*GeV2;
//   // for (int i = 2; i < mePartonData().size(); i++){
//   //   if (mePartonData()[i]->coloured()) scale = max(scale, meMomenta()[i].perp2());
//   // }

//   Energy2 scale = 1E12*GeV2;
//   for (int i = 2; i < mePartonData().size(); i++){
//     if (mePartonData()[i]->coloured()) {
//       Energy2 replace = max(meMomenta()[i].perp2(), 4*GeV2);
//       scale = min(scale, replace);
//     }
//   }

//   return scale;
// }

// Energy2 VBFNLOMEVVJJNeutralBase::renormalizationScale() const {
//   return 10000*GeV2;
//   // if ( theUserScale != ZERO )
//   //   return sqr(theUserScale);
//   // return lastSHat();

//   // Energy2 scale = 1*GeV2;
//   // for (int i = 2; i < mePartonData().size(); i++){
//   //   if (mePartonData()[i]->coloured()) scale = max(scale, meMomenta()[i].perp2());
//   // }
//   // return scale;

//   Energy2 scale = 1E12*GeV2;
//   for (int i = 2; i < mePartonData().size(); i++){
//     if (mePartonData()[i]->coloured()) {
//       Energy2 replace = max(meMomenta()[i].perp2(), 4*GeV2);
//       scale = min(scale, replace);
//     }
//   }
//   return scale;
// }

bool VBFNLOMEVVJJNeutralBase::requestedAsIncoming1(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming1() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming1()) return true;
  if (p->id() < 0 && !Incoming1()) return true;
  return false;
}

bool VBFNLOMEVVJJNeutralBase::requestedAsIncoming2(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming2() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming2()) return true;
  if (p->id() < 0 && !Incoming2()) return true;
  return false;
}

bool VBFNLOMEVVJJNeutralBase::allowedDiagram(DiagPtr di, int gluonIndex) const {

  Ptr<Tree2toNDiagram>::tcptr dptr = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>( di );

  PDPtr Wplus = getParticleData(ParticleID::Wplus); 
  PDPtr Wminus = getParticleData(ParticleID::Wminus); 
  PDPtr Z0 = getParticleData(ParticleID::Z0); 

  bool processMatches = false;
  //for (vector<PDVector>::const_iterator i = subProcesses().begin(); i !=subProcesses().end(); i++){
    bool subProcessMatches = true;
    PDVector::const_iterator j=subProcess().legs.begin();
    if ( (**j).id() != (dptr->incoming()).first->id() ) subProcessMatches = false;
    j++;
    if ( (**j).id() != (dptr->incoming()).second->id() ) subProcessMatches = false;
    unsigned int c = 0;
    for ( j=subProcess().legs.begin()+2; j !=subProcess().legs.end(); j++, c++){
      if ( (**j).id() != dptr->outgoing()[c]->id() ) subProcessMatches = false;
    }
    processMatches = processMatches || subProcessMatches;
    //}
  if (!processMatches) return false;
  
  //cerr << "processMatches=true\n" << flush;
  if ( gluonIndex == -1 ){
    assert( (dptr->incoming()).first->id() != 21 && (dptr->incoming()).second->id() != 21 );
    //    if ( (dptr->incoming()).first->id() != 2 || (dptr->incoming()).second->id() != 2 ) return false;
    if ( theCurrent == neutral
	 && ( (dptr->incoming()).first != dptr->outgoing()[0]
	      || (dptr->incoming()).second != dptr->outgoing()[1] ) ) return false;
    if ( theCurrent == charged
	 &&  ( (dptr->incoming()).first != SU2Helper::SU2CC(dptr->outgoing()[0])
	       || (dptr->incoming()).second != SU2Helper::SU2CC(dptr->outgoing()[1]) ) ) return false;
    if ( ( (dptr->incoming()).first->id() < 0 && theIncoming1)
	 || ((dptr->incoming()).second->id() < 0 && theIncoming2) ) return false;
    if ( ( (dptr->incoming()).first->id() > 0 && !theIncoming1)
	 || ((dptr->incoming()).second->id() > 0 && !theIncoming2) ) return false;
    // bool firstIsOnList = false;
    // bool secondIsOnList = false;
    // for ( PDVector::const_iterator q = theQuarkFlavours.begin();
    // 	  q != theQuarkFlavours.end(); ++q ){
    //   if ( abs((**q).id()) == abs((dptr->incoming()).first->id() ) ) firstIsOnList = true;
    //   if ( abs((**q).id()) == abs( (dptr->incoming()).second->id() ) ) secondIsOnList = true;
    // }
    // if (firstIsOnList && secondIsOnList) return true;
    // else return false;
    return true;
  }
  else if ( gluonIndex == 0) {
    //cerr << "enter GI=0\n" << flush;
    if ( (dptr->incoming()).first->id() != 21 ) return false;
    if ( theCurrent == neutral
	 && ( dptr->outgoing()[0] != dptr->outgoing()[2]->CC()
	      || (dptr->incoming()).second != dptr->outgoing()[1] ) ) return false;
    if ( theCurrent == charged
	 && ( dptr->outgoing()[0] != SU2Helper::SU2CC(dptr->outgoing()[2]->CC())
	      || (dptr->incoming()).second != SU2Helper::SU2CC(dptr->outgoing()[1]) ) ) return false;
    if ( ((dptr->incoming()).second->id() < 0 && theIncoming2) ) return false;
    if ( ((dptr->incoming()).second->id() > 0 && !theIncoming2) ) return false;
    // bool firstIsOnList = false;
    // bool secondIsOnList = false;
    // for ( PDVector::const_iterator q = theQuarkFlavours.begin();
    // 	  q != theQuarkFlavours.end(); ++q ){
    //   if ( abs((**q).id()) == abs(dptr->outgoing()[0]->id()) ) firstIsOnList = true;
    //   if ( abs((**q).id()) == abs( (dptr->incoming()).second->id() ) ) secondIsOnList = true;
    // }
    // if (firstIsOnList && secondIsOnList) return true;
    // else return false;
    return true;
    
  }
  else if ( gluonIndex == 1) {
    if ( (dptr->incoming()).second->id() != 21 ) return false;
    if ( theCurrent == neutral
	 && ( (dptr->incoming()).first != dptr->outgoing()[0]
	      || dptr->outgoing()[1] != dptr->outgoing()[2]->CC() ) ) return false;
    if ( theCurrent == charged
	 && ( (dptr->incoming()).first != SU2Helper::SU2CC(dptr->outgoing()[0])
	      || dptr->outgoing()[1] != SU2Helper::SU2CC(dptr->outgoing()[2]->CC()) ) ) return false;
    if ( ( (dptr->incoming()).first->id() < 0 && theIncoming1)) return false;
    if ( ( (dptr->incoming()).first->id() > 0 && !theIncoming1)) return false;
    // bool firstIsOnList = false;
    // bool secondIsOnList = false;
    // for ( PDVector::const_iterator q = theQuarkFlavours.begin();
    // 	  q != theQuarkFlavours.end(); ++q ){
    //   if ( abs((**q).id()) == abs( (dptr->incoming()).first->id() ) ) firstIsOnList = true;
    //   if ( abs((**q).id()) == abs(dptr->outgoing()[1]->id()) ) secondIsOnList = true;
    // }
    // if (firstIsOnList && secondIsOnList) return true;
    // else return false;
    return true;
  }
 else
    return false;
}

void VBFNLOMEVVJJNeutralBase::doinit(){
  VBFNLOMEBase::doinit();
}

void VBFNLOMEVVJJNeutralBase::doinitrun(){
  VBFNLOMEBase::doinitrun();
}

AbstractClassDescription<VBFNLOMEVVJJNeutralBase> VBFNLOMEVVJJNeutralBase::initVBFNLOMEVVJJNeutralBase;
// Definition of the static class description member.

void VBFNLOMEVVJJNeutralBase::persistentOutput(PersistentOStream & os) const {
  os << theIncoming1 << theIncoming2 << theCurrent;
}

void VBFNLOMEVVJJNeutralBase::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming1 >> theIncoming2 >> theCurrent;
}


void VBFNLOMEVVJJNeutralBase::Init() {

  static ClassDocumentation<VBFNLOMEVVJJNeutralBase> documentation
    ("VBFNLOMEVVJJNeutralBase");

  // static Switch<VBFNLOMEVVJJNeutralBase,bool> interfaceIncoming1
  //   ("Incoming1",
  //    "Set to true/false for parton 1 being particle/antiparticle.",
  //    &VBFNLOMEVVJJNeutralBase::theIncoming1, true, true, false);
  // static SwitchOption interfaceIncoming1Particle
  //   (interfaceIncoming1,
  //    "Particle",
  //    "Parton 1 is considered to be a particle.",
  //    true);
  // static SwitchOption interfaceIncoming1Antiparticle
  //   (interfaceIncoming1,
  //    "Antiparticle",
  //    "Parton 1 is considered to be an antiparticle.",
  //    false);

  // static Switch<VBFNLOMEVVJJNeutralBase,bool> interfaceIncoming2
  //   ("Incoming2",
  //    "Set to true/false for parton 2 being particle/antiparticle",
  //    &VBFNLOMEVVJJNeutralBase::theIncoming2, true, true, false);
  // static SwitchOption interfaceIncoming2Particle
  //   (interfaceIncoming2,
  //    "Particle",
  //    "Parton 2 is considered to be a particle.",
  //    true);
  // static SwitchOption interfaceIncoming2Antiparticle
  //   (interfaceIncoming2,
  //    "Antiparticle",
  //    "Parton 2 is considered to be an antiparticle.",
  //    false);

  // static Switch<VBFNLOMEVVJJNeutralBase,int> interfaceCurrent
  //   ("Current",
  //    "Choose the exchanged current for this matrix element.",
  //    &VBFNLOMEVVJJNeutralBase::theCurrent, 0, false, false);
  // static SwitchOption interfaceCurrentNeutral
  //   (interfaceCurrent,
  //    "Neutral",
  //    "Z0 and photon exchange are allowed.",
  //    neutral);
  // static SwitchOption interfaceCurrentCharged
  //   (interfaceCurrent,
  //    "Charged",
  //    "W+ and W- exchange are allowed.",
  //    charged);
}
