// -*- C++ -*-
//
// VBFNLOAmplitudeBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFNLOAmplitudeBase class.
//

#include "VBFNLOAmplitudeBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "VBFNLOCommonBlocks.h"
#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

VBFNLOAmplitudeBase::VBFNLOAmplitudeBase() 
  : MatchboxAmplitude(){
}

VBFNLOAmplitudeBase::~VBFNLOAmplitudeBase() {}

void VBFNLOAmplitudeBase::initCouplings(){}

void VBFNLOAmplitudeBase::doinit() {
  MatchboxAmplitude::doinit();
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

void VBFNLOAmplitudeBase::doinitrun() {
  MatchboxAmplitude::doinitrun();
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


//AbstractClassDescription<VBFNLOAmplitudeBase> VBFNLOAmplitudeBase::initVBFNLOAmplitudeBase;
// Definition of the static class description member.

void VBFNLOAmplitudeBase::persistentOutput(PersistentOStream &) const {}

void VBFNLOAmplitudeBase::persistentInput(PersistentIStream &, int) {}


void VBFNLOAmplitudeBase::Init() {

  static ClassDocumentation<VBFNLOAmplitudeBase> documentation
    ("VBFNLOAmplitudeBase");
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<VBFNLOAmplitudeBase,MatchboxAmplitude>
describeVBFNLOAmplitudeBase("Herwig::VBFNLOAmplitudeBase", "HwMatchbox.so HwMatchboxVBFNLO.so");
