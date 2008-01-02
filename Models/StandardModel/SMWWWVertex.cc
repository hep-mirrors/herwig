// -*- C++ -*-
//
// SMWWWVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
  os << _theSM << _zfact; 
}

void SMWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _zfact;
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
  if(q2!=_q2last) {
    _couplast = sqrt(4.0*Constants::pi*_theSM->alphaEM(q2));
    _q2last=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )          setNorm(_couplast);
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) )     setNorm(-_couplast);
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) || 
          (ida== 23 && idb==-24 && idc== 24) || 
          (ida== 24 && idb== 23 && idc==-24) )     setNorm(_couplast*_zfact);
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) || 
          (ida== 23 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 23 && idc== 24) )     setNorm(-_couplast*_zfact);
  else
    throw Helicity::HelicityConsistencyError() 
      << "SMWWWVertex::setCoupling "
      << "Invalid particles in WWW Vertex" 
      << Exception::runerror;
}

SMWWWVertex::SMWWWVertex() : _zfact(0.),_couplast(0.),_q2last(0.*GeV2) {
  // particles
  vector<int> first,second,third;
  first.push_back(24);
  second.push_back(-24);
  third.push_back(22);
  first.push_back(24);
  second.push_back(-24);
  third.push_back(23);
  setList(first,second,third);
}

void SMWWWVertex::doinit() throw(InitException) {
  orderInGem(1);
  orderInGs(0);
  VVVVertex::doinit();
  // factor for the Z vertex
  _theSM = generator()->standardModel();
  _zfact = sqrt((1.-_theSM->sin2ThetaW())/_theSM->sin2ThetaW());
}
