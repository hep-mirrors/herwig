// -*- C++ -*-
//
// VBFNLOVirtualMEVVJJNeutral.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFNLOVirtualMEVVJJNeutral class.
//

#include "VBFNLOVirtualMEVVJJNeutral.h"
#include "VBFNLOMEVVJJNeutralBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/Throw.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "VBFNLOCommonBlocks.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;

VBFNLOVirtualMEVVJJNeutral::VBFNLOVirtualMEVVJJNeutral() 
  : MatchboxInsertionOperator(), theCurrent(neutral), theIncoming1(true),
    theIncoming2(true) {}

VBFNLOVirtualMEVVJJNeutral::~VBFNLOVirtualMEVVJJNeutral() {}

IBPtr VBFNLOVirtualMEVVJJNeutral::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOVirtualMEVVJJNeutral::fullclone() const {
  return new_ptr(*this);
}


double VBFNLOVirtualMEVVJJNeutral::me2() const {
  int mePartonSign[6]={1,1,1,1,1,1};

  if (mePartonData()[0]->id() < 0) mePartonSign[0]=-1;
  if (mePartonData()[2]->id() < 0) mePartonSign[1]=-1;
  if (mePartonData()[1]->id() < 0) mePartonSign[2]=-1;
  if (mePartonData()[3]->id() < 0) mePartonSign[3]=-1;

  if (mePartonData()[mePartonData().size()-1]->id() < 0) mePartonSign[4]=-1; //changed mePartonData indices here
  if (mePartonData()[mePartonData().size()-2]->id() < 0) mePartonSign[5]=-1; //b

  SPLITCB.ALLSUBPROCS = false; 
  SPLITCB.SUBPROCID = BornMEVVJJNeutral()->getSubprocessID();
 
  int nlo=-1;
  CGLOBALL.DOVIRTUALS = true;

  int bosdec=0;

  double uucc, uuss, ddcc, ddss, udsc, ducs;

  double pbar[14][4];
  double qbar[5];
  BornMEVVJJNeutral()->initProcess(nlo);
  BornMEVVJJNeutral()->prepareMomenta(pbar,qbar);

  nlo=0;
  CGLOBALL.DOVIRTUALS = false;
  // cerr << "calling from virtual\n" << flush;
  BornMEVVJJNeutral()->VbfnloMe2(pbar,mePartonSign,qbar,1,nlo,bosdec,
	uucc,uuss,ddcc,ddss,udsc,ducs);
  
  double result[6];
  
  result[0]=uucc;
  result[1]=uuss;
  result[2]=ddcc;
  result[3]=ddss;
  result[4]=udsc;
  result[5]=ducs;

  // double cvirt = pow(Constants::pi,2)/3.0-7.0;
  // //double cvirt = 5.0*pow(Constants::pi,2)/3.0-16.0;
  // // return born*lastAlphaS()/Constants::pi*4.0/3.0*(9.0-pow(Constants::pi,2)+cvirt);
  // cerr << "evaluated ME as " << result[SPLITCB.SUBPROCID-1]*lastAlphaS()/Constants::pi*4.0/3.0*(cvirt-pow(Constants::pi,2)/3) << "\n" << flush;
  // return result[SPLITCB.SUBPROCID-1]*lastAlphaS()/Constants::pi*4.0/3.0*(cvirt-pow(Constants::pi,2)/3);


  return result[SPLITCB.SUBPROCID-1] - lastBorn()->lastXComb().lastME2()*(4.*Constants::pi/9.)*lastAlphaS() ;
}

bool VBFNLOVirtualMEVVJJNeutral::apply(const cPDVector& partons) const {
  if (theIncoming1 && partons[0]->id() < 0) return false;
  if (!theIncoming1 && partons[0]->id() > 0) return false;
  if (theIncoming2 && partons[1]->id() < 0) return false;
  if (!theIncoming2 && partons[1]->id() > 0) return false;
					
  for (int i = 0; i < 4; i++){
    if (abs(partons[i]->id())==0) return false;
    if (abs(partons[i]->id())>4) return false;
  }
  if ( theCurrent == neutral){
    if (partons[0]->id() != partons[2]->id()) return false;
    if (partons[1]->id() != partons[3]->id()) return false;
  }
  if ( theCurrent == charged ){
    if (partons[0]->id() != SU2Helper::SU2CC(partons[2])->id()) return false;
    if (partons[1]->id() != SU2Helper::SU2CC(partons[3])->id()) return false;
  }
  return true;
    
}

void VBFNLOVirtualMEVVJJNeutral::persistentOutput(PersistentOStream & os) const {
  os << theCurrent << theIncoming1
     << theIncoming2;
}

void VBFNLOVirtualMEVVJJNeutral::persistentInput(PersistentIStream & is, int) {
  is >> theCurrent >> theIncoming1
     >> theIncoming2;
}

ClassDescription<VBFNLOVirtualMEVVJJNeutral> VBFNLOVirtualMEVVJJNeutral::initVBFNLOVirtualMEVVJJNeutral;
// Definition of the static class description member.

void VBFNLOVirtualMEVVJJNeutral::Init() {

  static ClassDocumentation<VBFNLOVirtualMEVVJJNeutral> documentation
    ("VBFNLOVirtualMEVVJJNeutral is a class for all types of virtual",
     "corrections to vvjj-type matrix elements.");

  // static Reference<VBFNLOVirtualMEVVJJNeutral,VBFNLOMEVVJJNeutralBase> interfaceBornME
  //   ("BornME",
  //    "The corresponding Born matrix element to be used",
  //    &VBFNLOVirtualMEVVJJNeutral::theBornMEVVJJNeutral, false, false, true, false, false);

  static Switch<VBFNLOVirtualMEVVJJNeutral,int> interfaceCurrent
    ("Current",
     "Choose the exchanged current for this matrix element.",
     &VBFNLOVirtualMEVVJJNeutral::theCurrent, 0, false, false);
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

  static Switch<VBFNLOVirtualMEVVJJNeutral,bool> interfaceIncoming1
    ("Incoming1",
     "Set to true/false for parton 1 being particle/antiparticle",
     &VBFNLOVirtualMEVVJJNeutral::theIncoming1, true, true, false);
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

  static Switch<VBFNLOVirtualMEVVJJNeutral,bool> interfaceIncoming2
    ("Incoming2",
     "Set to true/false for parton 2 being particle/antiparticle",
     &VBFNLOVirtualMEVVJJNeutral::theIncoming2, true, true, false);
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
}

