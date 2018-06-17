// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWSSVertex class.
//

#include "SSWWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SSWWSSVertex::SSWWSSVertex() : sw_(0.), cw_(0.), q2last_(), couplast_(0.), 
			       ulast_(0), dlast_(0), gblast1_(0), gblast2_(0),
			       factlast_(0.) {
  colourStructure(ColourStructure::DELTA);
}

IBPtr SSWWSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSWWSSVertex::fullclone() const {
  return new_ptr(*this);
}

void SSWWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_  << cw_ << stau_ << stop_ << sbottom_;
}

void SSWWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> stau_ >> stop_ >> sbottom_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWWSSVertex,VVSSVertex>
describeHerwigSSWWSSVertex("Herwig::SSWWSSVertex", "HwSusy.so");

void SSWWSSVertex::Init() {

  static ClassDocumentation<SSWWSSVertex> documentation
    ("The SSWWSSVertex class implements the coupling of two"
     " weak bosons and two scalar fermions");

}

void SSWWSSVertex::doinit() {
  // gamma gamma, gamma Z0 and Z0 Z0 PAIRS
  for(long ib1=22;ib1<24;++ib1) {
    for(long ib2=22;ib2<24;++ib2) {
      // sleptons
      long istep = ib1==23&&ib2==23 ? 1 : 2;
      for(long ix=1000011;ix<1000017;ix+=istep) {
	addToList(ib1,ib2,ix,-ix);
      }
      for(long ix=2000011;ix<2000017;ix+=2) {
	addToList(ib1,ib2,ix,-ix);
      }
      // squarks
      for(long ix=1000001;ix<1000007;++ix) {
	addToList(ib1,ib2,ix,-ix);
      }
      for(long ix=2000001;ix<2000007;++ix) {
	addToList(ib1,ib2,ix,-ix);
      }
      // L/R mixing if Z
      if(ib2==23) {
	//L-Rbar stau
	addToList(ib1,ib2,1000015,-2000015);
	//Lbar-R stau
	addToList(ib1,ib2,-1000015,2000015);
	//L-Rbar stop
	addToList(ib1,ib2,1000006,-2000006);
	//Lbar-R stop
	addToList(ib1,ib2,-1000006,2000006);
	//L-Rbar sbottom
	addToList(ib1,ib2,1000005,-2000005);
	//Lbar-R sbottom
	addToList(ib1,ib2,-1000005,2000005);
      }
    }
    // W- gamma/Z0
    // LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(-24,ib1,ix+1,-ix);
    }
    // 1-2 stop sbottom
    addToList(-24,ib1,1000006,-2000005);
    // 2-1 stop sbottom
    addToList(-24,ib1,2000006,-1000005);
    // 2-2 stop sbottom
    addToList(-24,ib1,2000006,-2000005);
    // LL-sleptons
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(-24,ib1,-ix,ix+1);
    }
    //2-L stau
    addToList(-24,ib1,-2000015,1000016);
    // W+ gamma/Z0
    // LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(24,ib1,-(ix+1),ix);
    }
    // 1-2 stop sbottom
    addToList(24,ib1,-1000006,2000005);
    // 2-1 stop sbottom
    addToList(24,ib1,-2000006,1000005);
    // 2-2 stop sbottom
    addToList(24,ib1,-2000006,2000005);
    // LL-sleptons
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(24,ib1,ix,-ix-1);
    }
    //2-L stau
    addToList(24,ib1,2000015,-1000016);
  }
  // finally WW pairs
  // sleptons
  for(long ix=1000011;ix<=1000017;++ix) {
    addToList(24,-24,ix,-ix);
  }
  addToList(24,-24,1000015,-2000015);
  addToList(24,-24,2000015,-1000015);
  // squarks
  for(long ix=1000001;ix<=1000007;++ix) {
    addToList(24,-24,ix,-ix);
  }
  addToList(24,-24,1000005,-2000005);
  addToList(24,-24,2000005,-1000005);
  addToList(24,-24,1000006,-2000006);
  addToList(24,-24,2000006,-1000006);
  // couplings etc
  orderInGem(2);
  orderInGs(0);
  VVSSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSWWSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  sw_ = sqrt(sin2ThetaW());
  cw_ = sqrt( 1. - sqr(sw_) );
  stop_ = model->stopMix();
  sbottom_ = model->sbottomMix();
  stau_ = model->stauMix();
  if(!stop_ || !stau_ || !sbottom_)
    throw InitException() << "SSWWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << stop_ << " sbottom: " << sbottom_
			  << " stau: " << stau_ << Exception::abortnow;
}

void SSWWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3,
			       tcPDPtr part4) {
  long boson[2] = {abs(part1->id()),abs(part2->id())};
  assert( boson[0] == ParticleID::Wplus || boson[0] == ParticleID::Z0 ||
	  boson[0] == ParticleID::gamma );
  assert( boson[1] == ParticleID::Wplus || boson[1] == ParticleID::Z0 ||
	  boson[1] == ParticleID::gamma );
  long sf1(abs(part3->id())),sf2(abs(part4->id()));

  assert( (sf1 >= 1000001 && sf1 <= 1000006) 
  	  || (sf1 >= 1000011 && sf1 <= 1000016)
  	  || (sf1 >= 2000001 && sf1 <= 2000006)
  	  || (sf1 >= 2000011 && sf1 <= 2000016) );
  
  assert( (sf2 >= 1000001 && sf2 <= 1000006) 
  	  || (sf2 >= 1000011 && sf2 <= 1000016)
  	  || (sf2 >= 2000001 && sf2 <= 2000006)
  	  || (sf2 >= 2000011 && sf2 <= 2000016) );

  if( sf1 % 2 != 0 ) swap(sf1, sf2);
  if(boson[0]>boson[1]) swap(boson[0],boson[1]);
  if( sf1 != ulast_ || sf2 != dlast_ ||
      boson[0] != gblast1_ || boson[1] != gblast2_) {
    gblast1_ = boson[0];
    gblast2_ = boson[1];
    ulast_ = sf1;
    dlast_ = sf2;
    double ef1 = getParticleData(sf1)->charge()/eplus;
    double ef2 = getParticleData(sf2)->charge()/eplus;
    // photon pairs are simplest
    if(gblast1_ == ParticleID::gamma &&
       gblast2_ == ParticleID::gamma ) {
      factlast_ = 2.*sqr( ef1 );
    }
    // otherwise we need the helicity states
    else {
      //determine which helicity state
      unsigned int alpha(sf1/1000000 - 1), beta(sf2/1000000 - 1);
      //mixing factors
      Complex m1a(0.), m1b(0.);
      if( sf1 == ParticleID::SUSY_t_1 || sf1 == ParticleID::SUSY_t_2 )
  	m1a = (*stop_)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_b_1 || sf1 == ParticleID::SUSY_b_2 )
  	m1a = (*sbottom_)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_tau_1minus || 
  	       sf1 == ParticleID::SUSY_tau_2minus )
  	m1a = (*stau_)(alpha, 0);
      else
  	m1a = (alpha == 0) ? Complex(1.) : Complex(0.);
      
      if( sf2 == ParticleID::SUSY_t_1 || sf2 == ParticleID::SUSY_t_2 )
  	m1b = (*stop_)(beta, 0);
      else if( sf2 == ParticleID::SUSY_b_1 || sf2 == ParticleID::SUSY_b_2 )
  	m1b = (*sbottom_)(beta, 0);
      else if( sf2 == ParticleID::SUSY_tau_1minus || 
  	       sf2 == ParticleID::SUSY_tau_2minus )
  	m1b = (*stau_)(beta, 0);
      else
  	m1b = (beta == 0) ? Complex(1.) : Complex(0.);
      // if either boson is a W
      if(gblast2_==ParticleID::Wplus) {
	// WW
	if(gblast1_==ParticleID::Wplus) {
	  factlast_ = 0.5*m1a*m1b/sqr(sw_);
	}
	// gamma W
	else if(gblast1_==ParticleID::gamma) {
	  factlast_ = sqrt(0.5)*m1a*m1b/sw_*(ef1+ef2);
	}
	// Z0 W
	else if(gblast1_==ParticleID::Z0) {
	  factlast_ = -sqrt(0.5)*m1a*m1b/cw_*(ef1+ef2);
	}
      }
      else {
	// compute the Z coupling
	factlast_=1.;
	if(gblast1_==ParticleID::Z0||gblast2_==ParticleID::Z0) {
	  if( sf1 == ParticleID::SUSY_nu_eL || sf1 == ParticleID::SUSY_nu_muL ||
	      sf1 == ParticleID::SUSY_nu_tauL ) {
	    factlast_ = 0.5/cw_/sw_;
	  }
	  else {
	    double lmda =  sf2 % 2 == 0 ? 1. : -1.;
	    factlast_ = lmda*m1a*m1b;
	    if( alpha == beta) factlast_ -= 2.*ef1*sqr(sw_);
	    factlast_ *= 0.5/cw_/sw_; 
	  }
	}
	// photon Z
	if(gblast1_ == ParticleID::gamma &&
	   gblast2_ == ParticleID::Z0 ) {
	  factlast_ *= 2.*ef1; 
	}
	// Z pairs
	else if(gblast1_ == ParticleID::Z0 &&
		gblast2_ == ParticleID::Z0 ) {
	  factlast_ *= 2.*factlast_;
	}
      }
    }
  }
  if( q2 != q2last_ || couplast_==0. ) {
    q2last_ = q2;
    couplast_ = sqr(electroMagneticCoupling(q2));
  }
  norm(couplast_*factlast_);
}

