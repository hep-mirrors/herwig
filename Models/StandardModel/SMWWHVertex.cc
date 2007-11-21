// -*- C++ -*-
//
// SMWWHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WWHVertex class.
//
#include "SMWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMWWHVertex::SMWWHVertex() 
  : _couplast(0.), _q2last(), _mw(), _zfact(), _sw() {
  // particles
  vector<int> first,second,third;
  first.push_back(24);first.push_back(23);
  second.push_back(-24);second.push_back(23);
  third.push_back(25);third.push_back(25);
  setList(first,second,third);
}

void SMWWHVertex::doinit() throw(InitException) {
  // parameters
  _theSM = generator()->standardModel();
  _mw=getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _sw = _theSM->sin2ThetaW();
  _zfact = 1./(1.-_sw);
  _sw=sqrt(_sw);
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  // base class
  VVSVertex::doinit();
}
    
void SMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _zfact << _sw;
}

void SMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _zfact >> _sw;
}
    
ClassDescription<SMWWHVertex>SMWWHVertex::initSMWWHVertex;
// Definition of the static class description member.

void SMWWHVertex::Init() {
  static ClassDocumentation<SMWWHVertex> documentation
    ("The SMWWHVertex class is the implementation"
     " of the helicity amplitude calculation for the coupling of the Standard"
     " Model electroweak gauge bosons to the Higgs.");
}

void SMWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr) {
  int ibos=abs(a->id());
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = UnitRemoval::InvE * 
	sqrt(4.0*Constants::pi*alpha)*_mw/_sw;
    _q2last=q2;
  }
  if(ibos==24)      setNorm(_couplast);
  else if(ibos==23) setNorm(_couplast*_zfact);
  else {
    throw HelicityConsistencyError() << "SMWWHVertex::setCoupling "
				     << "Invalid particles in WWH Vertex" 
				     << Exception::warning;
    setNorm(0.);
  }
}
