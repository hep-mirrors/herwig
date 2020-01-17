// -*- C++ -*-
//
// RSModelVVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelWWWGRVertex class.
//

#include "RSModelWWWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelWWWGRVertex::RSModelWWWGRVertex() 
  : kappa_(ZERO), _couplast(0.), 
    _q2last(ZERO), _zfact(0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void RSModelWWWGRVertex::doinit() {
  addToList(24,-24, 22, 39);
  addToList(24,-24, 23, 39);
  VVVTVertex::doinit();
  _zfact = sqrt((1.-sin2ThetaW())/sin2ThetaW());
  // set the graviton coupling 
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() 
      << "Must have RSModel in RSModelWWWGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelWWWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << _zfact;
}
void RSModelWWWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> _zfact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelWWWGRVertex,VVVTVertex>
describeHerwigRSModelWWWGRVertex("Herwig::RSModelWWWGRVertex", "HwRSModel.so");

void RSModelWWWGRVertex::Init() {
 static ClassDocumentation<RSModelWWWGRVertex> documentation
    ("The RSModelWWWGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the WWWGR vertex
void RSModelWWWGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
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
     (ida== 24 && idb== 22 && idc==-24) )
    norm(Complex(_couplast*kappa_*UnitRemoval::E));
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) ||
	  (ida== 22 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 22 && idc== 24) )
    norm(-Complex(_couplast*kappa_*UnitRemoval::E));
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) ||
	  (ida== 23 && idb==-24 && idc== 24) || 
	  (ida== 24 && idb== 23 && idc==-24) )
    norm(Complex(_couplast*_zfact*kappa_*UnitRemoval::E));
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) ||
	  (ida== 23 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 23 && idc== 24) )
    norm(-Complex(_couplast*_zfact*kappa_*UnitRemoval::E));
  else assert(false);
}
