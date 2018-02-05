// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWHHVertex class.
//

#include "SSWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MSSM.h"

using namespace Herwig;

SSWWHHVertex::SSWWHHVertex()  : couplast_(0.), q2last_(ZERO) {
  orderInGs(0);
  orderInGem(2);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SSWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SSWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void SSWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSWWHHVertex,VVSSVertex>
describeHerwigSSWWHHVertex("Herwig::SSWWHHVertex", "HwSusy.so");

void SSWWHHVertex::Init() {

  static ClassDocumentation<SSWWHHVertex> documentation
    ("The SSWWHHVertex class implements the coupling of two Higgs bosons and"
     "two electroweak vector bosons in the MSSM.");

}

void SSWWHHVertex::doinit() {
  int id[3]={25,35,36};
  for(unsigned int ix=0;ix<3;++ix) {
    addToList( 24,-24,id[ix],id[ix]);
    addToList( 23, 23,id[ix],id[ix]);
    addToList( 22, 24,id[ix],-37);
    addToList( 22,-24,id[ix], 37);
    addToList( 23, 24,id[ix],-37);
    addToList( 23,-24,id[ix], 37);
  }
  addToList( 24,-24, 37,-37);
  addToList( 23, 23, 37,-37);
  addToList( 22, 23, 37,-37);
  addToList( 22, 22, 37,-37);
  VVSSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw Exception() 
      << "SSWWHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  coup_.resize(11);
  double sw2 = sin2ThetaW();
  double sw  = sqrt(sw2);
  double cw2 = 1.-sw2;
  double cw  = sqrt(cw2);
  double c2w = cw2-sw2;
  double sinalp = sin(model->higgsMixingAngle());
  double cosalp = sqrt(1. - sqr(sinalp));
  double tanbeta = model->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt( 1. - sqr(sinbeta) );
  double sinbma = sinbeta*cosalp - cosbeta*sinalp;
  double cosbma = cosbeta*cosalp + sinbeta*sinalp;
  // WWHH
  coup_[0] = 0.5/sw2;
  // ZZH0H0
  coup_[1] = 0.5/sw2/cw2;
  // ZZH+H-
  coup_[2] = 0.5*sqr(c2w)/cw2/sw2;
  // Z W h0 H+
  coup_[3] =-0.5/cw*cosbma;
  // Z W H0 H+
  coup_[4] = 0.5/cw*sinbma;
  // Z W A0 H+
  coup_[5] =-Complex(0.,0.5)/cw;
  // A A H+H-
  coup_[6] = 2.;
  // A Z H+H-
  coup_[7] = c2w/sw/cw;
  // A W h0 H+
  coup_[8] = 0.5*cosbma/sw;
  // A W H0 H+
  coup_[9] =-0.5*sinbma/sw;
  // A W A0 H+
  coup_[10] = Complex(0.,0.5)/sw;
}

void SSWWHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			       tcPDPtr part3, tcPDPtr part4) {
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(electroMagneticCoupling(q2));
    q2last_=q2;
  }
  int ibos1 = part1->id(), ibos2 = part2->id();
  int isca1 = part3->id(), isca2 = part4->id();
  if(abs(ibos1)==abs(ibos2)) {
    if(abs(ibos1)==ParticleID::Wplus) {
      norm(couplast_*coup_[0]);
    }
    else if(ibos1==ParticleID::Z0) {
      if(abs(isca1)==ParticleID::Hplus)
	norm(couplast_*coup_[2]);
      else
	norm(couplast_*coup_[1]);
    }
    else if(ibos1==ParticleID::gamma) {
	norm(couplast_*coup_[6]);
    }
    else
      assert(false);
  }
  else if(abs(ibos1)==ParticleID::Wplus ||
	  abs(ibos2)==ParticleID::Wplus) {
    if(abs(ibos1)==ParticleID::gamma ||
       abs(ibos2)==ParticleID::gamma) {
      if(abs(isca1)==ParticleID::h0 ||
	 abs(isca2)==ParticleID::h0) {
	norm(couplast_*coup_[8]);
      }
      else if(abs(isca1)==ParticleID::H0 ||
	      abs(isca2)==ParticleID::H0) {
	norm(couplast_*coup_[9]);
      }
      else if(abs(isca1)==ParticleID::A0 ||
	      abs(isca2)==ParticleID::A0) {
	if(isca1==ParticleID::Hplus ||
	   isca2==ParticleID::Hplus) {
	  norm(couplast_*     coup_[10] );
	}
	else {
	  norm(couplast_*conj(coup_[10]));
	}
      }
      else
	assert(false);
    }
    else {
      if(abs(isca1)==ParticleID::h0 ||
	 abs(isca2)==ParticleID::h0) {
	norm(couplast_*coup_[3]);
      }
      else if(abs(isca1)==ParticleID::H0 ||
	      abs(isca2)==ParticleID::H0) {
	norm(couplast_*coup_[4]);
      }
      else if(abs(isca1)==ParticleID::A0 ||
	      abs(isca2)==ParticleID::A0) {
	if(isca1==ParticleID::Hplus ||
	   isca2==ParticleID::Hplus) {
	  norm(couplast_*     coup_[5] );
	}
	else {
	  norm(couplast_*conj(coup_[5]));
	}
      }
      else
	assert(false);
    }
  }
  else {
    norm(couplast_*coup_[7]);
  }
}
