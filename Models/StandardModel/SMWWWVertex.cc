// -*- C++ -*-
//
// SMWWWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWVertex class.
//

#include "SMWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SMWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _zfact; 
}

void SMWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _zfact;
}

ClassDescription<SMWWWVertex>
SMWWWVertex::initSMWWWVertex;
// Definition of the static class description member.

void SMWWWVertex::Init() {
  static ClassDocumentation<SMWWWVertex> documentation
    ("The SMWWWVertex class is the implementation of the "
     "Standard Model triple electroweak boson coupling.");
  
}
    
// couplings for the WWW vertex
void SMWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || 
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )          norm(_couplast);
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) )     norm(-_couplast);
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) || 
          (ida== 23 && idb==-24 && idc== 24) || 
          (ida== 24 && idb== 23 && idc==-24) )     norm(_couplast*_zfact);
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) || 
          (ida== 23 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 23 && idc== 24) )     norm(-_couplast*_zfact);
  else
    throw Helicity::HelicityConsistencyError() 
      << "SMWWWVertex::setCoupling "
      << "Invalid particles in WWW Vertex"
      << a->PDGName() << " " << b->PDGName() << " " << c->PDGName() 
      << Exception::runerror;
}

SMWWWVertex::SMWWWVertex() : _zfact(0.),_couplast(0.), 
			     _q2last(sqr(Constants::MaxEnergy)) {
  orderInGem(1);
  orderInGs(0);
}

void SMWWWVertex::doinit() {
  addToList(24, -24, 22);
  addToList(24, -24, 23);
  VVVVertex::doinit();
  // factor for the Z vertex
  double sw2=sin2ThetaW();
  _zfact = sqrt((1.-sw2)/sw2);
}
