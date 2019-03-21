// -*- C++ -*-
//
// RSModelFFWGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFWGRVertex class.
//

#include "RSModelFFWGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;
using namespace ThePEG;

RSModelFFWGRVertex::RSModelFFWGRVertex() 
  : charge_(17,0.), gl_(17,0.), gr_(17,0.),
    ckm_(3,vector<Complex>(3,0.0)),
    couplast_(0.), q2last_(ZERO), kappa_(ZERO) {
  orderInGem(2);
  orderInGs (0);
}

void RSModelFFWGRVertex::doinit() {
  for(int ix=11;ix<17;++ix) {
    addToList(-ix,ix,22,39);
    addToList(-ix,ix,23,39);
  }
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,22,39);
    addToList(-ix,ix,23,39);
  }
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      addToList(-ix, iy, -24,39);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix+1, -24,39);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      addToList(-ix, iy, 24,39);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1, ix, 24,39);
  }
  FFVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) throw Exception() 
	      << "Must have RSModel in RSModelFFWGRVertex::doinit()"
	      << Exception::runerror;
  double sw2 = sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    charge_[2*ix-1]  = hwRS->ed();
    charge_[2*ix ]   = hwRS->eu();
    charge_[2*ix+9 ] = hwRS->ee();
    charge_[2*ix+10] = hwRS->enu();
    gl_[2*ix-1]  = fact*(hwRS->vd()  + hwRS->ad() );
    gl_[2*ix ]   = fact*(hwRS->vu()  + hwRS->au() );
    gl_[2*ix+9 ] = fact*(hwRS->ve()  + hwRS->ae() );
    gl_[2*ix+10] = fact*(hwRS->vnu() + hwRS->anu());
    gr_[2*ix-1]  = fact*(hwRS->vd()  - hwRS->ad() );
    gr_[2*ix ]   = fact*(hwRS->vu()  - hwRS->au() );
    gr_[2*ix+9 ] = fact*(hwRS->ve()  - hwRS->ae() );
    gr_[2*ix+10] = fact*(hwRS->vnu() - hwRS->anu());
  }
  kappa_ = 2./hwRS->lambda_pi();
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
	ckm_[ix][iy]=CKM[ix][iy];
      }
    }
  }
  else {
    throw Exception() << "Must have access to the Herwig::StandardCKM object"
		      << "for the CKM matrix in SMFFWVertex::doinit()"
		      << Exception::runerror;
  }
}

void RSModelFFWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << charge_ << gl_ << gr_ << ounit(kappa_,InvGeV) << ckm_;
}

void RSModelFFWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> charge_ >> gl_ >> gr_ >> iunit(kappa_,InvGeV) >> ckm_;
}

ClassDescription<RSModelFFWGRVertex> RSModelFFWGRVertex::initRSModelFFWGRVertex;
// Definition of the static class description member.

void RSModelFFWGRVertex::Init() {
  static ClassDocumentation<RSModelFFWGRVertex> documentation
    ("The RSModelFFWGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}

// FFWGR coupling
void RSModelFFWGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr bb,
				     tcPDPtr cc, tcPDPtr) {
  // work out the particles
  int iferm= abs(aa->id());
  int ibos = abs(cc->id());
  Complex coup;
  // overall factor
  assert( ibos >= 22 && ibos <= 24 );
  if( q2last_ != q2 || couplast_ == 0. ) {
    couplast_ = electroMagneticCoupling(q2);
    q2last_ = q2;
  }
  // photon
  if(ibos==22) {
    // alpha
    coup = Complex(UnitRemoval::E * kappa_ * couplast_);
    // _charge of particle
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    coup *= charge_[iferm];
    left (1.);
    right(1.);
  }
  // Z boson
  else if(ibos==23) {
    coup = Complex(UnitRemoval::E * kappa_ * couplast_);
    // _charge of particle
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    left (gl_[iferm]);
    right(gr_[iferm]);
  }
  else if(ibos==24) {
    coup = Complex(UnitRemoval::E * kappa_ * couplast_) * 
      sqrt(0.5) / sqrt(sin2ThetaW());
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
      left(ckm_[iu-1][id-1]);
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
  // set the coupling
  norm(coup);
}
