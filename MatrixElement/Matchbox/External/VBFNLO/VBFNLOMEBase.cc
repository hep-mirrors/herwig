// -*- C++ -*-
//
// VBFNLOMEBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFNLOMEBase class.
//

#include "VBFNLOMEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "VBFNLOCommonBlocks.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

VBFNLOMEBase::VBFNLOMEBase() 
  : MatchboxMEBase(), theUserScale(0*GeV), theInteger(0){
}

VBFNLOMEBase::~VBFNLOMEBase() {}

void VBFNLOMEBase::initCouplings(){

}

void VBFNLOMEBase::doinit() {
  MatchboxMEBase::doinit();
  BKOPIN.XMH = getParticleData(ParticleID::h0)->mass()/GeV;
  BKOPIN.XMT = getParticleData(ParticleID::t)->mass()/GeV;
  double bmass = 4.855;
  QUARKMASSES.XMB = bmass;//getParticleData(ParticleID::b)->mass()/GeV;
  QUARKMASSES.XMC = 1.65;//getParticleData(ParticleID::c)->mass()/GeV;
  BKOPIN.ALFAS = SM().alphaS();
  CGLOBALL.SLHA_SWITCH = false;
  SPLITCB.PRINTOUTPUT = false;

  //EWSCHEME=3 in VBFNLO
  BKOPIN.GF = (SM().fermiConstant())*GeV2;
  BKOPIN.XMW = getParticleData(ParticleID::Wplus)->mass()/GeV;
  BKOPIN.XMZ = getParticleData(ParticleID::Z0)->mass()/GeV;
  BKOPIN.ALFA = -1.; //SM().alphaEM();
  BKOPIN.SIN2W = -1.; //SM().sin2ThetaW();

  ANOM_SWITCH.WITH_ANOM = false;
  ANOMHIGGS.WITH_ANOMHIGGS = false;
  KK_SWITCH.WITH_KK = false;

  CGLOBALI.EWSCHEME = 3;
  
  CLEARWIDTHS();
  
  double e,g2,s,c,z,w,q,g;
  
  SETEWPARA(e,g2,s,c,z,w,q,g);

  KOPPLN(0,e,g2,s,c,z,w,q,g);
  CTRANS(bmass);
}

void VBFNLOMEBase::doinitrun() {
  MatchboxMEBase::doinitrun();
  BKOPIN.XMH = getParticleData(ParticleID::h0)->mass()/GeV;
  BKOPIN.XMT = getParticleData(ParticleID::t)->mass()/GeV;
  double bmass = 4.855;
  QUARKMASSES.XMB = bmass; //getParticleData(ParticleID::b)->mass()/GeV;
  QUARKMASSES.XMC = 1.65; //getParticleData(ParticleID::c)->mass()/GeV;
  BKOPIN.ALFAS = SM().alphaS();
  CGLOBALL.SLHA_SWITCH = false;
  SPLITCB.PRINTOUTPUT = false;

  //EWSCHEME=3 in VBFNLO
  BKOPIN.GF = (SM().fermiConstant())*GeV2;
  BKOPIN.XMW = getParticleData(ParticleID::Wplus)->mass()/GeV;
  BKOPIN.XMZ = getParticleData(ParticleID::Z0)->mass()/GeV;
  BKOPIN.ALFA = -1.; //SM().alphaEM();
  BKOPIN.SIN2W = -1.; //SM().sin2ThetaW();

  ANOM_SWITCH.WITH_ANOM = false;
  ANOMHIGGS.WITH_ANOMHIGGS = false;
  KK_SWITCH.WITH_KK = false;

  CGLOBALI.EWSCHEME = 3;

  CLEARWIDTHS();

  double e,g2,s,c,z,w,q,g;
  
  SETEWPARA(e,g2,s,c,z,w,q,g);

  KOPPLN(0,e,g2,s,c,z,w,q,g);

  CTRANS(bmass);
}


AbstractClassDescription<VBFNLOMEBase> VBFNLOMEBase::initVBFNLOMEBase;
// Definition of the static class description member.

void VBFNLOMEBase::persistentOutput(PersistentOStream & os) const {
  os << theQuarkFlavours << ounit(theUserScale,GeV) << theInteger;
}

void VBFNLOMEBase::persistentInput(PersistentIStream & is, int) {
  is >> theQuarkFlavours >> iunit(theUserScale,GeV) >> theInteger;
}


void VBFNLOMEBase::Init() {

  static ClassDocumentation<VBFNLOMEBase> documentation
    ("VBFNLOMEBase");

  static RefVector<VBFNLOMEBase,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &VBFNLOMEBase::theQuarkFlavours, -1, false, false, true, true, false);

  static Parameter<VBFNLOMEBase,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &VBFNLOMEBase::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<VBFNLOMEBase,int> interfaceInteger
    ("Integer",
     "Generic integer switch for debugging",
     &VBFNLOMEBase::theInteger, 1, 0, 0, 0,
     false, false, Interface::lowerlim);

}






