// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWHHVertex class.
//

#include "LHTPWHHVertex.h"
#include "LHTPModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHTPWHHVertex::LHTPWHHVertex() : 
  coupLast_(0.), q2Last_(ZERO), coup_(11) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHTPWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHTPWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWHHVertex,VSSVertex>
describeHerwigLHTPWHHVertex("Herwig::LHTPWHHVertex", "HwLHTPModel.so");

void LHTPWHHVertex::Init() {

  static ClassDocumentation<LHTPWHHVertex> documentation
    ("The LHTPWHHVertex class implements the coupling of a pair of Higgs"
     " bosons to an electroweak gauge boson in the Little"
     " Higgs model with T-parity.");

}

void LHTPWHHVertex::doinit() {
  // photon
  addToList( 22, 37,-37);
  addToList( 22, 38,-38);
  // Z0
  addToList( 23, 37,-37);
  addToList( 23, 38,-38);
  addToList( 23, 35, 36);
  // W+
  addToList( 24, 35,-37);
  addToList( 24, 36,-37);
  addToList( 24, 37,-38);
  // W-
  addToList(-24, 35, 37);
  addToList(-24, 36, 37);
  addToList(-24,-37, 38);
  // A_H
  addToList( 32, 25, 36);
  // Z_H
  addToList( 33, 25, 36);
  // W_H
  addToList( 34, 25,-37);
  addToList(-34, 25, 37);
  VSSVertex::doinit();
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double vf(model->vev()/model->f());
  coup_[ 0] = 1.;
  coup_[ 1] = 2.;
  coup_[ 2] =-sw/cw;
  coup_[ 3] = (1.-2.*sw2)/cw/sw;
  coup_[ 4] =-Complex(0.,1.)/cw/sw;
  coup_[ 5] = sqrt(0.5)/sw;
  coup_[ 6] = Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 7] = 1./sw;
  coup_[ 8] = Complex(0.,1)*sqrt(0.5)*vf/3./cw;
  coup_[ 9] = Complex(0.,1)*sqrt(0.5)*vf/3./sw;
  coup_[10] =-vf/6./sw;
}

void LHTPWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
				tcPDPtr particle2, tcPDPtr particle3) {
  if( q2 != q2Last_ || coupLast_==0.) {
    q2Last_ = q2;
    coupLast_ = electroMagneticCoupling(q2);
  }
  int ibos = particle1->id();
  int isc1 = particle2->id();
  int isc2 = particle3->id();
  if(ibos==ParticleID::gamma) {
    if(isc1==37) 
      norm(coup_[0]*coupLast_);
    else if(isc1==38)
      norm(coup_[1]*coupLast_);
    else if(isc1==-37) 
      norm(-coup_[0]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[1]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Z0) {
    if(isc1==37) 
      norm(coup_[2]*coupLast_);
    else if(isc1==38)
      norm(coup_[3]*coupLast_);
    else if(isc1==-37) 
      norm(-coup_[2]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[3]*coupLast_);
    else if(isc1==35)
      norm(coup_[4]*coupLast_);
    else if(isc2==35)
      norm(-coup_[4]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wplus) {
    if(isc1==35)
      norm(coup_[5]*coupLast_);
    else if(isc1==36)
      norm(coup_[6]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[7]*coupLast_);
    else if(isc2==35)
      norm(-coup_[5]*coupLast_);
    else if(isc2==36)
      norm(-coup_[6]*coupLast_);
    else if(isc2==-38)
      norm( coup_[7]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wminus) {
    if(isc1==35)
      norm(conj(coup_[5])*coupLast_);
    else if(isc1==36)
      norm(conj(coup_[6])*coupLast_);
    else if(isc1==38)
      norm(-conj(coup_[7])*coupLast_);
    else if(isc2==35)
      norm(-conj(coup_[5])*coupLast_);
    else if(isc2==36)
      norm(-conj(coup_[6])*coupLast_);
    else if(isc2==38)
      norm( conj(coup_[7])*coupLast_);
    else
      assert(false);
  }
  else if(ibos==32) {
    if(isc1==25)
      norm( coup_[8]*coupLast_);
    else if(isc2==25)
      norm(-coup_[8]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==33) {
    if(isc1==25)
      norm( coup_[9]*coupLast_);
    else if(isc2==25)
      norm(-coup_[9]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==34) {
    if(isc1==25)
      norm( coup_[10]*coupLast_);
    else if(isc2==25)
      norm(-coup_[10]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==-34) {
    if(isc1==25)
      norm( conj(coup_[10])*coupLast_);
    else if(isc2==25)
      norm(-conj(coup_[10])*coupLast_);
    else
      assert(false);
  }
  else
    assert(false);
}
