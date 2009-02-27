// -*- C++ -*-
//
// RSModelVVVGRVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelVVVGRVertex class.
//

#include "RSModelVVVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelVVVGRVertex::RSModelVVVGRVertex() : _couplast(2), _q2last(2) {
  vector<long> first,second,third,fourth;
  first.push_back(24);
  second.push_back(-24);
  third.push_back(22);
  fourth.push_back(39);

  first.push_back(24);
  second.push_back(-24);
  third.push_back(23);
  fourth.push_back(39);

  first.push_back(21);
  second.push_back(21);
  third.push_back(21);
  fourth.push_back(39);

  setList(first,second,third,fourth);
  _theKappa=InvEnergy();
  _zfact=0.;
}

void RSModelVVVGRVertex::doinit() {
  VVVTVertex::doinit();
  _zfact = sqrt((1.-generator()->standardModel()->sin2ThetaW())/
		generator()->standardModel()->sin2ThetaW());
  // set the graviton coupling 
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(hwRS){_theKappa=2./hwRS->lambda_pi();}
  else{throw InitException();}
}

void RSModelVVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_theKappa,InvGeV) << _zfact;
}
void RSModelVVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_theKappa,InvGeV) >> _zfact;
}

ClassDescription<RSModelVVVGRVertex> RSModelVVVGRVertex::initRSModelVVVGRVertex;
// Definition of the static class description member.

void RSModelVVVGRVertex::Init() {
 static ClassDocumentation<RSModelVVVGRVertex> documentation
    ("The RSModelVVVGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the VVVGR vertex
void RSModelVVVGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				 tcPDPtr c, tcPDPtr) {
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(ida==21 && idb==21 && idc==21) {
    if(q2!=_q2last[1]) {
      _couplast[1] = strongCoupling(q2);
      _q2last[1]=q2;
    }
    setNorm(Complex(_couplast[1]*_theKappa*UnitRemoval::E));
  }
  else {
    if(q2!=_q2last[0]) {
      _couplast[0] = electroMagneticCoupling(q2);
      _q2last[0]=q2;
    }
      // W- W+ photon and cylic perms
    if((ida==-24 && idb== 24 && idc== 22) ||
       (ida== 22 && idb==-24 && idc== 24) || 
       (ida== 24 && idb== 22 && idc==-24) )
      setNorm(Complex(_couplast[0]*_theKappa*UnitRemoval::E));
    // W+ W- photon (anticylic perms of above)
    else if((ida== 24 && idb==-24 && idc== 22) ||
	    (ida== 22 && idb== 24 && idc==-24) || 
	    (ida==-24 && idb== 22 && idc== 24) )
      setNorm(-Complex(_couplast[0]*_theKappa*UnitRemoval::E));
    // W- W+ Z and cylic perms
    else if((ida==-24 && idb== 24 && idc== 23) ||
	    (ida== 23 && idb==-24 && idc== 24) || 
	    (ida== 24 && idb== 23 && idc==-24) )
      setNorm(Complex(_couplast[0]*_zfact*_theKappa*UnitRemoval::E));
    // W+ W- Z (anticylic perms of above)
    else if((ida== 24 && idb==-24 && idc== 23) ||
	    (ida== 23 && idb== 24 && idc==-24) || 
	    (ida==-24 && idb== 23 && idc== 24) )
      setNorm(-Complex(_couplast[0]*_zfact*_theKappa*UnitRemoval::E));
    else throw HelicityConsistencyError() << "RSModelVVVGRVertex::setCoupling " 
					  << "Invalid particles in VVVGR Vertex" 
					  << Exception::runerror;
  }
}
