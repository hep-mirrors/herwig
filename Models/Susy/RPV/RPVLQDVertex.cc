// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVLQDVertex class.
//

#include "RPVLQDVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVLQDVertex::RPVLQDVertex() {
  orderInGem(1);
  orderInGs(0);
}

IBPtr RPVLQDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVLQDVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVLQDVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stop_ << sbot_ << stau_;
}

void RPVLQDVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stop_ >> sbot_ >> stau_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVLQDVertex,FFSVertex>
describeHerwigRPVLQDVertex("Herwig::RPVLQDVertex", "HwSusy.so HwRPV.so");

void RPVLQDVertex::Init() {

  static ClassDocumentation<RPVLQDVertex> documentation
    ("The RPVLQDVertex class implements the trilinear LQD"
     " coupling in R-parity violating models.");

}

void RPVLQDVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVLQDVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaLQD();
  stop_   = rpv->stopMix();
  sbot_   = rpv->sbottomMix();
  stau_   = rpv->stauMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
	// continue if zero
	if(lambda_[i][j][k]==0.) continue;
	// particles in the vertex
	// sneutrino
	addToList( -2*j-1 ,  2*k+1 , -1000012-2*i );
	addToList(  2*j+1 , -2*k-1 ,  1000012+2*i );
	// left slepton
	addToList( -2*j-2 ,  2*k+1 , -1000011-2*i );
	addToList(  2*j+2 , -2*k-1 , +1000011+2*i );
	if(i==2) {
	  addToList( -2*j-2 ,  2*k+1 , -2000011-2*i );
	  addToList(  2*j+2 , -2*k-1 , +2000011+2*i );
	}
	// up squark
	addToList( -2*i-11,  2*k+1 , -1000002-2*j );
	addToList(  2*i+11, -2*k-1 ,  1000002+2*j );
	if(j==2) {
	  addToList( -2*i-11,  2*k+1 , -2000002-2*j );
	  addToList(  2*i+11, -2*k-1 ,  2000002+2*j );
	}
	// left down squark
	addToList( -2*i-12,  2*k+1 , -1000001-2*j );
	addToList(  2*i+12, -2*k-1 ,  1000001+2*j );
	if(j==2) {
	  addToList( -2*i-12,  2*k+1 , -2000001-2*j );
	  addToList(  2*i+12, -2*k-1 ,  2000001+2*j );
	}
	// right down squark
	addToList(  2*i+12,  2*j+1 , -2000001-2*k );
	addToList( -2*i-12, -2*j-1 , +2000001+2*k );
	if(k==2) {
	  addToList(  2*i+12,  2*j+1 , -1000001-2*k );
	  addToList( -2*i-12, -2*j-1 , +1000001+2*k );
	}
	addToList(  2*i+11,  2*j+2 , -2000001-2*k );
	addToList( -2*i-11, -2*j-2 , +2000001+2*k );
	if(k==2) {
	  addToList(  2*i+11,  2*j+2 , -1000001-2*k );
	  addToList( -2*i-11, -2*j-2 , +1000001+2*k );
	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVLQDVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix(1.);
  // sneutrino case
  if( abs(islep) == ParticleID::SUSY_nu_eL  || 
      abs(islep) == ParticleID::SUSY_nu_muL || 
      abs(islep) == ParticleID::SUSY_nu_tauL) {
    i = (abs(islep)-1000012)/2;
    j = (abs(part1->id())-1)/2;
    k = (abs(part2->id())-1)/2;
    if( islep*part1->id() <0 ) swap(j,k);
    if(islep>0) {
      left ( 0.);
      right(-1.);
    }
    else {
      left (-1.);
      right( 0.);
    }
  }
  else if(abs(islep) == ParticleID::SUSY_e_Lminus|| 
	  abs(islep) == ParticleID::SUSY_mu_Lminus|| 
	  abs(islep) == ParticleID::SUSY_tau_1minus|| 
	  abs(islep) == ParticleID::SUSY_tau_2minus) { 
    if(abs(islep)<2000000) {
      i = (abs(islep)-1000011)/2;
      if(i==2) mix = (*stau_)(0,0);
    }
    else {
      i = 2;
      mix = (*stau_)(1,0);
      assert(abs(islep)==ParticleID::SUSY_tau_2minus);
    }
    if(abs(part1->id())%2==1) {
      j = (abs(part2->id())-2)/2;
      k = (abs(part1->id())-1)/2;
    }
    else {
      j = (abs(part1->id())-2)/2;
      k = (abs(part2->id())-1)/2;
    }
    if(islep<0) {
      left ( 1.);
      right( 0.);
    }
    else {
      left ( 0.);
      right( 1.);
    }
  }
  // left up squark
  else if(abs(islep) == ParticleID::SUSY_u_L  || 
	  abs(islep) == ParticleID::SUSY_c_L || 
	  abs(islep) == ParticleID::SUSY_t_1 || 
	  abs(islep) == ParticleID::SUSY_t_2) {
    if(abs(islep)<2000000) {
      j = (abs(islep)-1000002)/2;
      if(j==2) mix = (*stop_)(0,0); 
    }
    else {
      j = 2;
      mix = (*stop_)(1,0);
      assert(abs(islep)==ParticleID::SUSY_t_2);
    }
    if(part1->coloured()) {
      i = (abs(part2->id())-11)/2;
      k = (abs(part1->id())- 1)/2;
    }
    else {
      i = (abs(part1->id())-11)/2;
      k = (abs(part2->id())- 1)/2;
    }
    if(islep<0) {
      left ( 1.);
      right( 0.);
    }
    else {
      left ( 0.);
      right( 1.);
    }
  }
  else if(abs(islep) == ParticleID::SUSY_d_L || 
	  abs(islep) == ParticleID::SUSY_s_L || 
	  abs(islep) == ParticleID::SUSY_b_1 || 
	  abs(islep) == ParticleID::SUSY_b_2 ||
	  abs(islep) == ParticleID::SUSY_d_R || 
	  abs(islep) == ParticleID::SUSY_s_R) {
    // right down squark
    if(part1->id()*part2->id()>0) {
      if(part1->coloured()) {
	if(abs(part1->id())%2==1) {
	  i = (abs(part2->id())-12)/2;
	  j = (abs(part1->id())-1 )/2;
	  mix *= -1.;
	  assert(i>=0);
	}
	else {
	  i = (abs(part2->id())-11)/2;
	  j = (abs(part1->id())-2 )/2;
	  assert(i>=0);
	}
      }
      else {
	if(abs(part2->id())%2==1) {
	  i = (abs(part1->id())-12)/2;
	  j = (abs(part2->id())-1 )/2;
	  mix *= -1.;
	  assert(i>=0);
	}
	else {
	  i = (abs(part1->id())-11)/2;
	  j = (abs(part2->id())-2 )/2;
	  assert(i>=0);
	}
      }
      if(abs(islep)>2000000) {
	k = (abs(islep)-2000001)/2;
	if(k==2) mix *= (*sbot_)(1,1);
      }
      else {
	k = 2;
	mix *= (*sbot_)(0,1);
	assert(abs(islep)==ParticleID::SUSY_b_1);
      }
      if(islep<0) {
	left ( 0.);
	right( 1.);
      }
      else {
	left ( 1.);
	right( 0.);
      }
    }
    // left   down  squark
    else {
      if(abs(islep)<2000000) {
	j = (abs(islep)-1000001)/2;
	if(j==2) mix = (*sbot_)(0,0); 
      }
      else {
	j = 2;
	mix = (*sbot_)(1,0);
	assert(abs(islep)==ParticleID::SUSY_b_2);
      }
      if(part1->coloured()) {
	i = (abs(part2->id())-12)/2;
	k = (abs(part1->id())- 1)/2;
      }
      else {
	i = (abs(part1->id())-12)/2;
	k = (abs(part2->id())- 1)/2;
      }
      if(islep<0) {
	left (-1.);
	right( 0.);
      }
      else {
	left ( 0.);
	right(-1.);
      }
      assert(i>=0);
    }
  }
  else
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  norm(mix*lambda_[i][j][k]);
}
