#include "VBFNLOAmplitudeVVJJNeutralBase.h"
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

VBFNLOAmplitudeVVJJNeutralBase::VBFNLOAmplitudeVVJJNeutralBase()
  : VBFNLOAmplitudeBase(), theCurrent(neutral), theIncoming1(true),
    theIncoming2(true) {}

VBFNLOAmplitudeVVJJNeutralBase::~VBFNLOAmplitudeVVJJNeutralBase(){}

double VBFNLOAmplitudeVVJJNeutralBase::me2() const{
  int mePartonSign[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  int gluonIndex = 100;
  int meGluonSign = 1;

  for (unsigned int i = 0; i < mePartonData().size(); i++){
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

  //changed mePartonData indices here
  //because of VBFNLO convention
  if (mePartonData()[mePartonData().size()-1]->id() < 0) mePartonSign[4]=-1; 
  if (mePartonData()[mePartonData().size()-2]->id() < 0) mePartonSign[5]=-1; 

  SPLITCB.ALLSUBPROCS = false; 
  SPLITCB.SUBPROCID = getSubprocessID(gluonIndex);

  if (mePartonData()[0]->id() == 21) SPLITCB.GLUONID = 0;
  else if (mePartonData()[1]->id() == 21) SPLITCB.GLUONID = 1;
  else SPLITCB.GLUONID = -1;

  int nlo=0;
  CGLOBALL.DOVIRTUALS = false;

  double uucc, uuss, ddcc, ddss, udsc, ducs;

  double pbar[14][4];
  double qbar[5];
  prepareMomenta(pbar,qbar,gluonIndex);
  VbfnloMe2(pbar,mePartonSign,qbar,meGluonSign,nlo,
	uucc,uuss,ddcc,ddss,udsc,ducs);
  
  double result[6];
  
  result[0]=uucc;
  result[1]=uuss;
  result[2]=ddcc;
  result[3]=ddss;
  result[4]=udsc;
  result[5]=ducs;


  // for (int i = 0; i < meMomenta().size(); i++) {
  //    cerr << "p[" << i << "]=" << meMomenta()[i]/GeV << "\n" << flush;
  // }
  // cerr << "lastME2()= " << result[SPLITCB.SUBPROCID-1] << "\n" << flush;
  // cerr << "SPLITCB.SUBPROCID-1=" << SPLITCB.SUBPROCID-1 << "\n" << flush;
  // cerr << "gluonIndex = " << gluonIndex << "\n" << flush;
  // lastME2(result[SPLITCB.SUBPROCID-1]);
  return result[SPLITCB.SUBPROCID-1]; //lastME2();
  
}

double VBFNLOAmplitudeVVJJNeutralBase::oneLoopInterference() const {
  //for hjj production, the nlo correction can be simply evaluated as it factorizes with the born me:
  double born = me2();
  double cvirt = pow(Constants::pi,2)/3.0-7.0;
  return born*lastAlphaS()/Constants::pi*4.0/3.0*(cvirt-pow(Constants::pi,2)/3);

  //as soon as more complex processes are interfaced, the following should be used:

  // int mePartonSign[14]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  // if (mePartonData()[0]->id() < 0) mePartonSign[0]=-1;
  // if (mePartonData()[2]->id() < 0) mePartonSign[1]=-1;
  // if (mePartonData()[1]->id() < 0) mePartonSign[2]=-1;
  // if (mePartonData()[3]->id() < 0) mePartonSign[3]=-1;

  // if (mePartonData()[mePartonData().size()-1]->id() < 0) mePartonSign[4]=-1; //changed mePartonData indices here
  // if (mePartonData()[mePartonData().size()-2]->id() < 0) mePartonSign[5]=-1; //b

  // SPLITCB.ALLSUBPROCS = false; 
  // SPLITCB.SUBPROCID = getSubprocessID();
 
  // int nlo=-1;
  // CGLOBALL.DOVIRTUALS = true;

  // int bosdec=0;

  // double uucc, uuss, ddcc, ddss, udsc, ducs;

  // double pbar[14][4];
  // double qbar[5];
  // initProcess(nlo);
  // prepareMomenta(pbar,qbar);

  // VbfnloMe2(pbar,mePartonSign,qbar,1,nlo,bosdec,
  // 	uucc,uuss,ddcc,ddss,udsc,ducs);
  
  // double result[6];
  
  // result[0]=uucc;
  // result[1]=uuss;
  // result[2]=ddcc;
  // result[3]=ddss;
  // result[4]=udsc;
  // result[5]=ducs;

  // return result[SPLITCB.SUBPROCID-1] - born*(4.*Constants::pi/9.)*lastAlphaS() ;
}


bool VBFNLOAmplitudeVVJJNeutralBase::canHandle(const PDVector& data) const{

  if (data.size() == 5) {
    if (data[4]->id()!=25) return false;
    for (PDVector::const_iterator i = data.begin(); i != data.end(); i++)
      if ((**i).id()==21) return false;
    if (!allowedProcess(data)) return false;
    return true;
  }
  else if (data.size() == 6) {
    if (data[5]->id() != 25) return false;
    int gluonIndex = -1;
    for (PDVector::const_iterator i = data.begin(); i != data.end(); i++)
      if ((**i).id() == 21) {
	gluonIndex = i - data.begin();
	break;
      }
    if (gluonIndex > 1 && gluonIndex != 4) return false;
    //if (gluonIndex < 2) return false;
    if (!allowedProcess(data)) return false;
    return true;
  }
  return false;
}


int VBFNLOAmplitudeVVJJNeutralBase::getSubprocessID(const PDVector& data, int gluonIndex) const{

  //these are indices that make sure we get the proper subprocess
  //when we have gluons in the initial state for the real emission
  //matrix elements
  int U = 0;
  int L = 1;

  if (gluonIndex == 0) U = 4;
  if (gluonIndex == 1) L = 4;


  if ( SU2Helper::isSU2Up(data[U]) &&
       SU2Helper::isSU2Up(data[2]) &&
       SU2Helper::isSU2Up(data[L]) &&
       SU2Helper::isSU2Up(data[3]) 
       ){
    return 1;
  }
  else if ( SU2Helper::isSU2Up(data[U]) &&
	    SU2Helper::isSU2Up(data[2]) &&
	    SU2Helper::isSU2Down(data[L]) &&
	    SU2Helper::isSU2Down(data[3]) 
	    ){
    return 2;
  }
  else if ( SU2Helper::isSU2Down(data[U]) &&
	    SU2Helper::isSU2Down(data[2]) &&
	    SU2Helper::isSU2Up(data[L]) &&
	    SU2Helper::isSU2Up(data[3]) 
	    ){
    return 3;
  }
  else if ( SU2Helper::isSU2Down(data[U]) &&
	    SU2Helper::isSU2Down(data[2]) &&
	    SU2Helper::isSU2Down(data[L]) &&
	    SU2Helper::isSU2Down(data[3]) 
	    ){
    return 4;
  }
  else {
    if (gluonIndex == 0) {
      if (data[L]->id()>0){
	if (SU2Helper::isSU2Down(data[L]) &&
	    SU2Helper::isSU2Up(data[3])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Up(data[L]) &&
		 SU2Helper::isSU2Down(data[3])) {
	  return 6;
	}
      }
      else if (data[L]->id()<0) {
	if (SU2Helper::isSU2Up(data[L]) &&
	    SU2Helper::isSU2Down(data[3])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Down(data[L]) &&
		 SU2Helper::isSU2Up(data[3])) {
	  return 6;
	}
      }
    }
    else if (gluonIndex == 1) {
      if (data[U]->id()>0){
	if (SU2Helper::isSU2Down(data[U]) &&
	    SU2Helper::isSU2Up(data[2])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Up(data[U]) &&
		 SU2Helper::isSU2Down(data[2])) {
	  return 6;
	}
      }
      else if (data[U]->id()<0) {
	if (SU2Helper::isSU2Up(data[U]) &&
	    SU2Helper::isSU2Down(data[2])) {
	  return 5;
	}
	else if (SU2Helper::isSU2Down(data[U]) &&
		 SU2Helper::isSU2Up(data[2])) {
	  return 6;
	}
      }
    }
    else {	    
      if(data[U]->id()*data[L]->id() > 0) {
	if ( SU2Helper::isSU2Up(data[U]) &&
	     SU2Helper::isSU2Down(data[2]) &&
	     SU2Helper::isSU2Down(data[L]) &&
	     SU2Helper::isSU2Up(data[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Down(data[U]) &&
		  SU2Helper::isSU2Up(data[2]) &&
		  SU2Helper::isSU2Up(data[L]) &&
		  SU2Helper::isSU2Down(data[3]) 
		  ){
	  return 6;
	}
      }

      else if (data[U]->id()<0){
	if ( SU2Helper::isSU2Down(data[U]) &&
	     SU2Helper::isSU2Up(data[2]) &&
	     SU2Helper::isSU2Down(data[L]) &&
	     SU2Helper::isSU2Up(data[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Up(data[U]) &&
    		  SU2Helper::isSU2Down(data[2]) &&
		  SU2Helper::isSU2Up(data[L]) &&
    		  SU2Helper::isSU2Down(data[3]) 
		  ){
	  return 6;
	}
      }
      else if (data[L]->id()<0){
	if ( SU2Helper::isSU2Up(data[U]) &&
	     SU2Helper::isSU2Down(data[2]) &&
	     SU2Helper::isSU2Up(data[L]) &&
	     SU2Helper::isSU2Down(data[3]) 
	     ){
	  return 5;
	}
	else if ( SU2Helper::isSU2Down(data[U]) &&
		  SU2Helper::isSU2Up(data[2]) &&
		  SU2Helper::isSU2Down(data[L]) &&
		  SU2Helper::isSU2Up(data[3]) 
		  ){
	  return 6;
	}
      
      }
    }
  }
  return -1;
  // throw ThePEG::Exception() << "Could not determine subprocess ID."
  // 			    << " Got the following particle IDs:"
  // 			    <<  data[0]->id() << " "
  // 			    <<  data[1]->id() << " "
  // 			    <<  data[2]->id() << " "
  // 			    <<  data[3]->id() << " "
  // 			    <<  data[4]->id() << " "
  // 			    << ThePEG::Exception::abortnow;  
}


int VBFNLOAmplitudeVVJJNeutralBase::getSubprocessID(int gluonIndex) const{

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
  return -1;
  // throw ThePEG::Exception() << "Could not determine subprocess ID."
  // 			    << " Got the following particle IDs:"
  // 			    <<  mePartonData()[0]->id() << " "
  // 			    <<  mePartonData()[1]->id() << " "
  // 			    <<  mePartonData()[2]->id() << " "
  // 			    <<  mePartonData()[3]->id() << " "
  // 			    <<  mePartonData()[4]->id() << " "
  // 			    << ThePEG::Exception::abortnow;  
}

// Energy2 VBFNLOAmplitudeVVJJNeutralBase::factorizationScale() const {
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

// Energy2 VBFNLOAmplitudeVVJJNeutralBase::renormalizationScale() const {
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

bool VBFNLOAmplitudeVVJJNeutralBase::requestedAsIncoming1(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming1() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming1()) return true;
  if (p->id() < 0 && !Incoming1()) return true;
  return false;
}

bool VBFNLOAmplitudeVVJJNeutralBase::requestedAsIncoming2(PDPtr p) const {
  if (p->id() == 21) 
    throw ThePEG::Exception() << "isRequestedAsIncoming2() is not supposed to be used for gluons!"
			      << ThePEG::Exception::abortnow; 
  if (p->id() > 0 && Incoming2()) return true;
  if (p->id() < 0 && !Incoming2()) return true;
  return false;
}

bool VBFNLOAmplitudeVVJJNeutralBase::allowedDiagram(DiagPtr di, int gluonIndex) const {

  Ptr<Tree2toNDiagram>::tcptr dptr = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>( di );

  PDPtr Wplus = getParticleData(ParticleID::Wplus); 
  PDPtr Wminus = getParticleData(ParticleID::Wminus); 
  PDPtr Z0 = getParticleData(ParticleID::Z0); 

  if ( gluonIndex == -1 ){
    assert( (dptr->incoming()).first->id() != 21 && (dptr->incoming()).second->id() != 21 );
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
    return true;
  }
  else if ( gluonIndex == 0) {
    if ( (dptr->incoming()).first->id() != 21 ) return false;
    if ( theCurrent == neutral
	 && ( dptr->outgoing()[0] != dptr->outgoing()[2]->CC()
	      || (dptr->incoming()).second != dptr->outgoing()[1] ) ) return false;
    if ( theCurrent == charged
	 && ( dptr->outgoing()[0] != SU2Helper::SU2CC(dptr->outgoing()[2]->CC())
	      || (dptr->incoming()).second != SU2Helper::SU2CC(dptr->outgoing()[1]) ) ) return false;
    if ( ((dptr->incoming()).second->id() < 0 && theIncoming2) ) return false;
    if ( ((dptr->incoming()).second->id() > 0 && !theIncoming2) ) return false;
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
    return true;
  }
 else
    return false;
}

void VBFNLOAmplitudeVVJJNeutralBase::doinit(){
  VBFNLOAmplitudeBase::doinit();
}

void VBFNLOAmplitudeVVJJNeutralBase::doinitrun(){
  VBFNLOAmplitudeBase::doinitrun();
}

AbstractClassDescription<VBFNLOAmplitudeVVJJNeutralBase> VBFNLOAmplitudeVVJJNeutralBase::initVBFNLOAmplitudeVVJJNeutralBase;
// Definition of the static class description member.

void VBFNLOAmplitudeVVJJNeutralBase::persistentOutput(PersistentOStream & os) const {
  os << theIncoming1 << theIncoming2 << theCurrent;
}

void VBFNLOAmplitudeVVJJNeutralBase::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming1 >> theIncoming2 >> theCurrent;
}


void VBFNLOAmplitudeVVJJNeutralBase::Init() {

  static ClassDocumentation<VBFNLOAmplitudeVVJJNeutralBase> documentation
    ("VBFNLOAmplitudeVVJJNeutralBase");

  static Switch<VBFNLOAmplitudeVVJJNeutralBase,bool> interfaceIncoming1
    ("Incoming1",
     "Set to true/false for parton 1 being particle/antiparticle.",
     &VBFNLOAmplitudeVVJJNeutralBase::theIncoming1, true, true, false);
  static SwitchOption interfaceIncoming1Particle
    (interfaceIncoming1,
     "Particle",
     "Parton 1 is considered to be a particle.",
     true);
  static SwitchOption interfaceIncoming1Antiparticle
    (interfaceIncoming1,
     "Antiparticle",
     "Parton 1 is considered to be an antiparticle.",
     false);

  static Switch<VBFNLOAmplitudeVVJJNeutralBase,bool> interfaceIncoming2
    ("Incoming2",
     "Set to true/false for parton 2 being particle/antiparticle",
     &VBFNLOAmplitudeVVJJNeutralBase::theIncoming2, true, true, false);
  static SwitchOption interfaceIncoming2Particle
    (interfaceIncoming2,
     "Particle",
     "Parton 2 is considered to be a particle.",
     true);
  static SwitchOption interfaceIncoming2Antiparticle
    (interfaceIncoming2,
     "Antiparticle",
     "Parton 2 is considered to be an antiparticle.",
     false);

  static Switch<VBFNLOAmplitudeVVJJNeutralBase,int> interfaceCurrent
    ("Current",
     "Choose the exchanged current for this matrix element.",
     &VBFNLOAmplitudeVVJJNeutralBase::theCurrent, 0, false, false);
  static SwitchOption interfaceCurrentNeutral
    (interfaceCurrent,
     "Neutral",
     "Z0 and photon exchange are allowed.",
     neutral);
  static SwitchOption interfaceCurrentCharged
    (interfaceCurrent,
     "Charged",
     "W+ and W- exchange are allowed.",
     charged);
}
