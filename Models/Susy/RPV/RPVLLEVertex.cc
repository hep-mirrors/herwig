// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVLLEVertex class.
//

#include "RPVLLEVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVLLEVertex::RPVLLEVertex() {
  orderInGem(1);
  orderInGs(0);
}

IBPtr RPVLLEVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVLLEVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVLLEVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stau_;
}

void RPVLLEVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stau_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVLLEVertex,FFSVertex>
describeHerwigRPVLLEVertex("Herwig::RPVLLEVertex", "HwSusy.so HwRPV.so");

void RPVLLEVertex::Init() {

  static ClassDocumentation<RPVLLEVertex> documentation
    ("The RPVLLEVertex class implements the trilinear LLE"
     " coupling in R-parity violating models.");

}

void RPVLLEVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVLLEVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaLLE();
  stau_   = rpv->stauMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
	// continue if zero
	if(lambda_[i][j][k]==0.) continue;
	// particles in the vertex
	// sneutrino
	addToList( -2*i-11,  2*k+11, -1000012-2*j );
	addToList(  2*i+11, -2*k-11,  1000012+2*j );
	// left slepton
	addToList( -2*i-12,  2*k+11, -1000011-2*j );
	addToList(  2*i+12, -2*k-11, +1000011+2*j );
	if(j==2) {
	  addToList( -2*i-12,  2*k+11, -2000011-2*j );
	  addToList(  2*i+12, -2*k-11, +2000011+2*j );
	}
	// right slepton
	addToList(  2*i+12,  2*j+11, -2000011-2*k );
	addToList( -2*i-12, -2*j-11, +2000011+2*k );
	if(k==2) {
	  addToList(  2*i+12,  2*j+11, -1000011-2*k );
	  addToList( -2*i-12, -2*j-11, +1000011+2*k );
	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVLLEVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix=1.;
  // sneutrino case
  if( abs(islep) == ParticleID::SUSY_nu_eL  || 
      abs(islep) == ParticleID::SUSY_nu_muL || 
      abs(islep) == ParticleID::SUSY_nu_tauL) {
    j = (abs(islep) - 1000012)/2;
    i = (abs(part1->id())-11)/2;
    k = (abs(part2->id())-11)/2;
    if(part1->id()*islep<0) swap(i,k);
    if(islep<0) {
      left (1.);
      right(0.);
    }
    else {
      left (0.);
      right(1.);
    }
  }
  // charged slepton case
  else if( abs(islep) == ParticleID::SUSY_e_Lminus  || 
	   abs(islep) == ParticleID::SUSY_mu_Lminus || 
	   abs(islep) == ParticleID::SUSY_e_Rminus  || 
	   abs(islep) == ParticleID::SUSY_mu_Rminus || 
	   abs(islep) == ParticleID::SUSY_tau_1minus|| 
	   abs(islep) == ParticleID::SUSY_tau_2minus) {
    // right charged slepton
    if(part1->id()*part2->id()>0) {
      if(abs(part1->id())%2==0) {
	i = (abs(part1->id())-12)/2;
	j = (abs(part2->id())-11)/2;
      }
      else {
	i = (abs(part2->id())-12)/2;
	j = (abs(part1->id())-11)/2;
      }
      if(abs(islep)>2000000) {
	k = (abs(islep)-2000011)/2;
	if(k==2) mix = (*stau_)(1,1);
      }
      else {
	assert(abs(islep)==ParticleID::SUSY_tau_1minus);
	k = 2;
	mix = (*stau_)(0,1);
      }
      if(islep>0) {
	left (-1.);
	right(0.);
      }
      else {
	left (0.);
	right(-1.);
      }
    }
    // left charged
    else {
      if(abs(part1->id())%2==0) {
	i = (abs(part1->id())-12)/2;
	k = (abs(part2->id())-11)/2;
      }
      else {
	i = (abs(part2->id())-12)/2;
	k = (abs(part1->id())-11)/2;
      }
      if(abs(islep)<2000000) {
	j = (abs(islep)-1000011)/2;
	if(j==2) mix = (*stau_)(0,0);
      }
      else {
	assert(abs(islep)==ParticleID::SUSY_tau_2minus);
	j = 2;
	mix = (*stau_)(1,0);
      }
      if(islep<0) {
	left (-1.);
	right(0.);
      }
      else {
	left (0.);
	right(-1.);
      }
    }
  }
  else
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  norm(mix*lambda_[i][j][k]);
}
  
