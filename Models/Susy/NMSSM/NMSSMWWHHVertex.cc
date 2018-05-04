// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWWHHVertex class.
//

#include "NMSSMWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "NMSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NMSSMWWHHVertex::NMSSMWWHHVertex() : couplast_(0.),q2last_(ZERO),
				     sw_(0.), cw_(0.), sb_(0.), cb_(0.) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr NMSSMWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr NMSSMWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void NMSSMWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_ << cw_ << sb_ << cb_ << mixS_ << mixP_;
}

void NMSSMWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> sb_ >> cb_ >> mixS_ >> mixP_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMWWHHVertex,Helicity::VVSSVertex>
describeHerwigNMSSMWWHHVertex("Herwig::NMSSMWWHHVertex", "NMSSMWWHHVertex.so");

void NMSSMWWHHVertex::Init() {

  static ClassDocumentation<NMSSMWWHHVertex> documentation
    ("The NMSSMWWHHVertex class implements the coupling of"
     " two electroweak and two Higgs bosons in the NMSSM.");

}

void NMSSMWWHHVertex::doinit() {
  int scalar[3]={25,35,45};
  int pseudo[2]={36,46};
  // scalar higgs bosons
  for(unsigned int i=0;i<3;++i) {
    // pair of scalars
    for(unsigned int j=0;j<3;++j) {
      addToList( 24,-24,scalar[i],scalar[j]);
      addToList( 23, 23,scalar[i],scalar[j]);
    }
    // scalar charged
    addToList( 22, 24,scalar[i],-37);
    addToList( 22,-24,scalar[i], 37);
    addToList( 23, 24,scalar[i],-37);
    addToList( 23,-24,scalar[i], 37);
  }
  // pair of pseudoscalars
  for(unsigned int i=0;i<2;++i) {
    for(unsigned int j=0;j<2;++j) {
      addToList( 24,-24,pseudo[i],pseudo[j]);
      addToList( 23, 23,pseudo[i],pseudo[j]);
    }
    // pseudo charged
    addToList( 22, 24,pseudo[i],-37);
    addToList( 22,-24,pseudo[i], 37);
    addToList( 23, 24,pseudo[i],-37);
    addToList( 23,-24,pseudo[i], 37);
  }
  addToList( 24,-24, 37,-37);
  addToList( 23, 23, 37,-37);
  addToList( 22, 23, 37,-37);
  addToList( 22, 22, 37,-37);
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMWWHHVertex::doinit()"
			  << Exception::runerror;
  // sin and cos theta_W
  sw_ = sqrt(   sin2ThetaW());
  cw_ = sqrt(1.-sin2ThetaW());
  // get the mixing matrices
  mixS_ = model->CPevenHiggsMix();
  if(!mixS_) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWWHHVertex::doinit()" 
				   << Exception::runerror;
  mixP_ = model->CPoddHiggsMix();
  if(!mixP_) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMWWHHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  sb_ = sin(beta);
  cb_ = cos(beta);
  // base class
  VVSSVertex::doinit();
}
 
void NMSSMWWHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
				  tcPDPtr part3, tcPDPtr part4) {
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(electroMagneticCoupling(q2));
    q2last_=q2;
  }
  int ibos1 = part1->id(), ibos2 = part2->id();
  int isca1 = part3->id(), isca2 = part4->id();
  Complex fact(1.);
  if(abs(ibos1)==abs(ibos2)) {
    fact *= 0.5/sqr(sw_);
    if(ibos1==ParticleID::Z0) fact /= sqr(cw_);
    // charged
    if(abs(isca1)==37) {
      if(ibos1==ParticleID::Z0) fact *= sqr(sqr(cw_)-sqr(sw_));
      else if(ibos1==ParticleID::gamma) fact = 2.;
    }
    // pair of scalars
    else if((isca1-5)%10==0) {
      unsigned int i = (isca1-25)/10;
      unsigned int j = (isca2-25)/10;
      fact *= (*mixS_)(i,0)*(*mixS_)(j,0)+(*mixS_)(i,1)*(*mixS_)(j,1);
    }
    // pair of pseudoscalars
    else if((isca1-6)%10==0) {
      unsigned int i = (isca1-36)/10;
      unsigned int j = (isca2-36)/10;
      fact *= (*mixP_)(i,0)*(*mixP_)(j,0)+(*mixP_)(i,1)*(*mixP_)(j,1);
    }
    else
      assert(false);
  }
  else if(abs(ibos1)==ParticleID::Wplus ||
	  abs(ibos2)==ParticleID::Wplus) {
    if(abs(ibos1)==ParticleID::gamma ||
       abs(ibos2)==ParticleID::gamma) {
      fact *= -0.5/sw_;
    }
    else {
      fact *=  0.5/cw_;
    }
    if((isca1-5)%10==0) {
      unsigned int i = (isca1-25)/10;
      fact *= sb_*(*mixS_)(i,0) - cb_*(*mixS_)(i,1);
    }
    else if((isca2-5)%10==0) {
      unsigned int i = (isca2-25)/10;
      fact *= sb_*(*mixS_)(i,0) - cb_*(*mixS_)(i,1);
    }
    else if((isca1-6)%10==0) {
      unsigned int i = (isca1-36)/10;
      fact *= sb_*(*mixP_)(i,0) + cb_*(*mixP_)(i,1);
      fact *= isca2==ParticleID::Hplus ? -Complex(0.,1.) : Complex(0.,1.);
    }
    else if((isca2-6)%10==0) {
      unsigned int i = (isca2-36)/10;
      fact *= sb_*(*mixP_)(i,0) + cb_*(*mixP_)(i,1);
      fact *= isca1==ParticleID::Hplus ? -Complex(0.,1.) : Complex(0.,1.);
    }
    else
      assert(false);
  }
  else {
    fact = (sqr(cw_)-sqr(sw_))/cw_/sw_;
  }
  // set the coupling
  norm(couplast_*fact);
}
