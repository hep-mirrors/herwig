// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFPVertex class.
//

#include "LHTPFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "LHTPModel.h"

using namespace Herwig;

IBPtr LHTPFFPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFPVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << charge_ << coupd_ << coupu_ << coupe_ << coupnu_ 
     << TPreFactor_ << sL_ << cL_ << sR_ << cR_;
}

void LHTPFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> charge_ >> coupd_ >> coupu_ >> coupe_ >> coupnu_ 
     >> TPreFactor_ >> sL_ >> cL_ >> sR_ >> cR_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFPVertex,FFVVertex>
describeHerwigLHTPFFPVertex("Herwig::LHTPFFPVertex", "HwLHTPModel.so");

void LHTPFFPVertex::Init() {

  static ClassDocumentation<LHTPFFPVertex> documentation
    ("The LHTPFFPVertex class implements the coupling"
     " of the charged fermions to the photon in the Little Higgs"
     " model with T-parity.");

}

LHTPFFPVertex::LHTPFFPVertex() : 
  charge_(37,0.0), coupLast_(0.), q2Last_(-1.*GeV2),
  coupd_(0.), coupu_(0.), coupe_(0.), coupnu_(0.),
  TPreFactor_(0.), sL_(0.), cL_(1.), sR_(0.), cR_(1.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHTPFFPVertex::doinit() {
  // interactions with the photon
  // the quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 22);
  }
  // the leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-8,   8, 22);
  // the T-odd quarks
  for(long ix = 4000001;ix < 4000007; ++ix) {
    addToList(-ix,     ix, 22);
  }
  // the T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-4000008,  4000008, 22);
  // interactions with A_H
  // quarks and T-odd quark
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix - 4000000, ix          , 32);
    addToList(-ix          , ix + 4000000, 32);
  }
  // leptons and T-odd leptons (both charged leptons and neutrinos)
  for(int ix = 11; ix < 17; ++ix ) {
    addToList(-ix - 4000000, ix          , 32);
    addToList(-ix          , ix + 4000000, 32);
  }
  // T+T-A_H
  addToList(-4000008,       8, 32);
  addToList(      -8, 4000008, 32);
  // T-tA_H
  addToList(-4000008,       6, 32);
  addToList(      -6, 4000008, 32);
  // T-tA_H
  addToList(-4000006,       8, 32);
  addToList(      -8, 4000006, 32);
  // charges
  for(int ix = 1; ix < 16; ++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) charge_[ix] = double(ptemp->iCharge())/3.;
  }
  for(int ix = 4000001; ix < 4000016; ++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) charge_[ix-3999980] = double(ptemp->iCharge())/3.;
  }
  // couplings to A_H
  double sw = generator()->standardModel()->sin2ThetaW();
  double cw = sqrt(1.-sw);
  sw = sqrt(sw);
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFPVertex::doinit()"
				   << Exception::runerror;
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  sL_ = model->sinThetaL();
  cL_ = model->cosThetaL();
  sR_ = model->sinThetaR();
  cR_ = model->cosThetaR();
  // couplings of fermion T-odd fermion A_H
  coupd_  = -0.1*(cH/cw-5.*sH/sw);
  coupu_  = -0.1*(cH/cw+5.*sH/sw);
  coupe_  = -0.1*(cH/cw-5.*sH/sw);
  coupnu_ = -0.1*(cH/cw+5.*sH/sw);
  // couplings of T+T- A_H
  TPreFactor_ = 0.4*cH/cw;
  FFVVertex::doinit();
}

// coupling for FFP vertex
void LHTPFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -electroMagneticCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  long iferm=abs(a->id()),ibos(c->id());
  if(ibos == ParticleID::gamma) {
    if(iferm < 20) {
      left (charge_[iferm]);
      right(charge_[iferm]);
    }
    else {
      iferm-=3999980;
      left (charge_[iferm]);
      right(charge_[iferm]);
    }
  }
  else if(ibos == 32) {
    long ianti = abs(b->id());
    if(iferm>4000000) swap(iferm,ianti);
    assert(iferm<4000000&&ianti>4000000);
    if( iferm == 6 || iferm == 8 ) {
      if     (iferm==6&&ianti==4000006) {
	left (cL_*coupu_);
	right(0.);
      }
      else if(iferm==6&&ianti==4000008) {
	left (-TPreFactor_*sL_);
	right(-TPreFactor_*sR_);
      }
      else if(iferm==8&&ianti==4000006) {
	left (sL_*coupu_);
	right(0.);
      }
      else if(iferm==8&&ianti==4000008) {
	left ( TPreFactor_*cL_);
	right( TPreFactor_*cR_);
      }
      else
	assert(false);
    }
    // quarks (inclding top)
    else if(iferm <= 5) {
      if(iferm % 2 == 0) left(coupu_);
      else               left(coupd_);
      right(0.);
    }
    // leptons
    else {
      if(iferm %2 == 0) left(coupnu_);
      else              left(coupe_ );
      right(0.);
    }
  }
  else
    assert(false);
}
