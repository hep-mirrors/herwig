// -*- C++ -*-
//
// SMFFWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFWVertex class.
//

#include "SMFFWVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
    
SMFFWVertex::SMFFWVertex() : 
  _diagonal(false), _ckm(3,vector<Complex>(3,0.0)), 
  _couplast(0.), _q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SMFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _diagonal << _ckm;
}
  
void SMFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _diagonal >> _ckm;
}
  
void SMFFWVertex::doinit() {
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      bool isOff = iy/2 != (ix+1)/2;
      if ( isOff && _diagonal )
	continue;
      addToList(-ix, iy, -24);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix+1, -24);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      bool isOff = ix/2 != (iy+1)/2;
      if ( isOff && _diagonal )
	continue;
      addToList(-ix, iy, 24);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1, ix, 24);
  }
  ThePEG::Helicity::FFVVertex::doinit();
  if ( !_diagonal ) {
    Ptr<CKMBase>::transient_pointer CKM = generator()->standardModel()->CKM();
    // cast the CKM object to the HERWIG one
    ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
      hwCKM = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
					transient_const_pointer>(CKM);
    if(hwCKM) {
      vector< vector<Complex > > CKM;
      CKM = hwCKM->getUnsquaredMatrix(generator()->standardModel()->families());
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  _ckm[ix][iy]=CKM[ix][iy];
	}
      }
    }
    else {
      throw Exception() << "Must have access to the Herwig::StandardCKM object"
			<< "for the CKM matrix in SMFFWVertex::doinit()"
			<< Exception::runerror;
    }
  } else {
    _ckm = vector< vector<Complex > >(3,vector<Complex >(3,1.0));
  }
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SMFFWVertex,FFVVertex>
describeHerwigSMFFWVertex("Herwig::SMFFWVertex", "Herwig.so");
  
void SMFFWVertex::Init() {
  static ClassDocumentation<SMFFWVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the W boson to the Standard Model fermions");


  static Switch<SMFFWVertex,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &SMFFWVertex::_diagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
  
}
  
// coupling for FFW vertex
void SMFFWVertex::setCoupling(Energy2 q2, tcPDPtr aa, tcPDPtr bb, tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -sqrt(0.5)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  int ianti=abs(bb->id());
  // quarks
  if(iferm>=1 && iferm <=6) {
    int iu,id;
    // up type first
    if(iferm%2==0) {
      iu = iferm/2;
      id = (ianti+1)/2;
    }
    // down type first
    else {
      iu = ianti/2;
      id = (iferm+1)/2;
    }
    assert( iu>=1 && iu<=3 && id>=1 && id<=3);
    left(_ckm[iu-1][id-1]);
    right(0.);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    left(1.);
    right(0.);
  }
  else 
    assert(false);
}







