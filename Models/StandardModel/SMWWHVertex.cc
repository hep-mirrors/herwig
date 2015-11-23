// -*- C++ -*-
//
// SMWWHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
  : _couplast(0.), _q2last(ZERO), _mw(ZERO), _zfact(0.) {
  orderInGem(1);
  orderInGs(0);
}

void SMWWHVertex::doinit() {
  addToList(24,-24, 25);
  addToList(23, 23, 25);
  // parameters
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _zfact = 1./(1.-sin2ThetaW());
  // base class
  VVSVertex::doinit();
}
    
void SMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw,GeV) << _zfact;
}

void SMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw,GeV) >> _zfact;
}
    
ClassDescription<SMWWHVertex>SMWWHVertex::initSMWWHVertex;
// Definition of the static class description member.

void SMWWHVertex::Init() {
  static ClassDocumentation<SMWWHVertex> documentation
    ("The SMWWHVertex class is the implementation"
     " of the helicity amplitude calculation for the coupling of the Standard"
     " Model electroweak gauge bosons to the Higgs.");
}

void SMWWHVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr, tcPDPtr) {
  int ibos=abs(aa->id());
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = weakCoupling(q2) * UnitRemoval::InvE * _mw;
    _q2last=q2;
  }
  if(ibos==24)      norm(_couplast);
  else if(ibos==23) norm(_couplast*_zfact);
  else
    throw HelicityConsistencyError() << "SMWWHVertex::setCoupling "
				     << "Invalid particles in WWH Vertex" 
				     << Exception::runerror;
}
