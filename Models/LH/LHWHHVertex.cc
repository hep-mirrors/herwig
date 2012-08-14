// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWHHVertex class.
//

#include "LHWHHVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWHHVertex::LHWHHVertex() : 
  couplast_(0.), q2last_(ZERO), coup_(24) {
  orderInGs(0);
  orderInGem(1);
  // neutral
  addToList( 22, 37,-37);
  addToList( 22, 38,-38);
  addToList( 32, 37,-37);
  addToList( 32, 38,-38);
  addToList( 23, 37,-37);
  addToList( 23, 38,-38);
  addToList( 33, 37,-37);
  addToList( 33, 38,-38);
  addToList( 32, 25, 36);
  addToList( 32, 35, 36);
  addToList( 23, 25, 36);
  addToList( 23, 35, 36);
  addToList( 33, 25, 36);
  addToList( 33, 35, 36);
  // W+
  addToList( 24, 25,-37);
  addToList( 24, 35,-37);
  addToList( 24, 36,-37);
  addToList( 24, 37,-38);
  addToList( 34, 25,-37);
  addToList( 34, 35,-37);
  addToList( 34, 36,-37);
  addToList( 34, 37,-38);
  // W-
  addToList(-24, 25, 37);
  addToList(-24, 35, 37);
  addToList(-24, 36, 37);
  addToList(-24,-37, 38);
  addToList(-34, 25, 37);
  addToList(-34, 35, 37);
  addToList(-34, 36, 37);
  addToList(-34,-37, 38);
}

IBPtr LHWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWHHVertex,VSSVertex>
describeHerwigLHWHHVertex("Herwig::LHWHHVertex", "HwLHModel.so");

void LHWHHVertex::Init() {

  static ClassDocumentation<LHWHHVertex> documentation
    ("The LHWHHVertex class implements the coupling of a pair of Higgs"
     " bosons to an electroweak gauge boson in the Little Higgs model.");

}

void LHWHHVertex::doinit() {
  VSSVertex::doinit();
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double s0   (model->sinTheta0());
  double sP   (model->sinThetaP());
  double sPlus(model->sinThetaPlus());
  coup_[ 0] = 0.5/sw*(sqrt(2.)*s0-sPlus);
  coup_[ 1] = sqrt(0.5)/sw;
  coup_[ 2] = Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 3] = 1./sw;
  coup_[ 4] = 0.;
  coup_[ 5] = 0.;
  coup_[ 6] = 1.;
  coup_[ 7] = 2.;
  coup_[ 8] = Complex(0.,0.5)/cw/sw*(sP-2.*s0);	       
  coup_[ 9] =-Complex(0.,1.)/cw/sw;
  coup_[10] =-sw/cw;
  coup_[11] = (1.-2.*sw2)/cw/sw;
  coup_[12] =-0.25/sw*(sqr(c)-sqr(s))/s/c*(sqrt(2.)*s0-sPlus);
  coup_[13] =-sqrt(0.5)/sw*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[14] =-Complex(0.,1.)*sqrt(0.5)*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[15] =-0.5*(sqr(c)-sqr(s))/s/c/sw;
  coup_[16] =-Complex(0.,0.5)/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp*(sP-2.*s0);
  coup_[17] = Complex(0.,1.)/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[18] =-0.5*(sqr(cp)-sqr(sp))/sp/cp/cw;
  coup_[19] =-0.5*(sqr(cp)-sqr(sp))/sp/cp/cw;
  coup_[20] =-Complex(0.,0.5)/sw*0.5*(sqr(c)-sqr(s))/s/c*(sP-2.*s0);
  coup_[21] = Complex(0.,1.)/sw*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[22] = 0.;
  coup_[23] =-0.5/sw*(sqr(c)-sqr(s))/s/c;
}

void LHWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2, tcPDPtr particle3) {
  if( q2 != q2last_ || couplast_==0.) {
    q2last_ = q2;
    couplast_ = electroMagneticCoupling(q2);
  }
  int ibos = particle1->id();
  int isc1 = particle2->id();
  int isc2 = particle3->id();
  if(ibos==ParticleID::gamma) {
    if(isc1==37) 
      norm(coup_[6]*couplast_);
    else if(isc1==38)
      norm(coup_[7]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[6]*couplast_);
    else if(isc1==-38)
      norm(-coup_[7]*couplast_);
    else
      assert(false);
  }
  if(ibos==32) {
    if(isc1==37) 
      norm(coup_[18]*couplast_);
    else if(isc1==38)
      norm(coup_[19]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[18]*couplast_);
    else if(isc1==-38)
      norm(-coup_[19]*couplast_);
    else if(isc1==25)
      norm(coup_[16]*couplast_);
    else if(isc1==35)
      norm(coup_[17]*couplast_);
    else if(isc2==25)
      norm(-coup_[16]*couplast_);
    else if(isc2==35)
      norm(-coup_[17]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Z0) {
    if(isc1==37) 
      norm(coup_[10]*couplast_);
    else if(isc1==38)
      norm(coup_[11]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[10]*couplast_);
    else if(isc1==-38)
      norm(-coup_[11]*couplast_);
    else if(isc1==25)
      norm(coup_[8]*couplast_);
    else if(isc1==35)
      norm(coup_[9]*couplast_);
    else if(isc2==25)
      norm(-coup_[8]*couplast_);
    else if(isc2==35)
      norm(-coup_[9]*couplast_);
    else
      assert(false);
  }
  else if(ibos==33) {
    if(isc1==37) 
      norm(coup_[22]*couplast_);
    else if(isc1==38)
      norm(coup_[23]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[22]*couplast_);
    else if(isc1==-38)
      norm(-coup_[23]*couplast_);
    else if(isc1==25)
      norm(coup_[20]*couplast_);
    else if(isc1==35)
      norm(coup_[21]*couplast_);
    else if(isc2==25)
      norm(-coup_[20]*couplast_);
    else if(isc2==35)
      norm(-coup_[21]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wplus) {
    if(isc1==25)
      norm(coup_[0]*couplast_);
    else if(isc1==35)
      norm(coup_[1]*couplast_);
    else if(isc1==36)
      norm(coup_[2]*couplast_);
    else if(isc1==37)
      norm(coup_[3]*couplast_);
    else if(isc2==25)
      norm(-coup_[0]*couplast_);
    else if(isc2==35)
      norm(-coup_[1]*couplast_);
    else if(isc2==36)
      norm(-coup_[2]*couplast_);
    else if(isc2==37)
      norm(-coup_[3]*couplast_);
    else
      assert(false);
  }
  else if(ibos==34) {
    if(isc1==25)
      norm(coup_[12]*couplast_);
    else if(isc1==35)
      norm(coup_[13]*couplast_);
    else if(isc1==36)
      norm(coup_[14]*couplast_);
    else if(isc1==37)
      norm(coup_[15]*couplast_);
    else if(isc2==25)
      norm(-coup_[12]*couplast_);
    else if(isc2==35)
      norm(-coup_[13]*couplast_);
    else if(isc2==36)
      norm(-coup_[14]*couplast_);
    else if(isc2==37)
      norm(-coup_[15]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wminus) {
    if(isc1==25)
      norm(conj(coup_[0])*couplast_);
    else if(isc1==35)
      norm(conj(coup_[1])*couplast_);
    else if(isc1==36)
      norm(conj(coup_[2])*couplast_);
    else if(isc1==37)
      norm(conj(coup_[3])*couplast_);
    else if(isc2==25)
      norm(-conj(coup_[0])*couplast_);
    else if(isc2==35)
      norm(-conj(coup_[1])*couplast_);
    else if(isc2==36)
      norm(-conj(coup_[2])*couplast_);
    else if(isc2==37)
      norm(-conj(coup_[3])*couplast_);
    else
      assert(false);
  }
  else if(ibos==-34) {
    if(isc1==25)
      norm(conj(coup_[12])*couplast_);
    else if(isc1==35)
      norm(conj(coup_[13])*couplast_);
    else if(isc1==36)
      norm(conj(coup_[14])*couplast_);
    else if(isc1==37)
      norm(conj(coup_[15])*couplast_);
    else if(isc2==25)
      norm(-conj(coup_[12])*couplast_);
    else if(isc2==35)
      norm(-conj(coup_[13])*couplast_);
    else if(isc2==36)
      norm(-conj(coup_[14])*couplast_);
    else if(isc2==37)
      norm(-conj(coup_[15])*couplast_);
    else
      assert(false);
  }
}
