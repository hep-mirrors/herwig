// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UDDVertex class.
//

#include "RPVUDDVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVUDDVertex::RPVUDDVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::EPS);
}

IBPtr RPVUDDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVUDDVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVUDDVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stop_ << sbot_;
}

void RPVUDDVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stop_ >> sbot_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVUDDVertex,FFSVertex>
describeHerwigRPVUDDVertex("Herwig::RPVUDDVertex", "HwSusy.so HwRPV.so");

void RPVUDDVertex::Init() {

  static ClassDocumentation<RPVUDDVertex> documentation
    ("The RPVUDDVertex class implements the trilinear UDD"
     " coupling in R-parity violating models.");

}

void RPVUDDVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVUDDVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaUDD();
  stop_   = rpv->stopMix();
  sbot_   = rpv->sbottomMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
 	// continue if zero
 	if(lambda_[i][j][k]==0.) continue;
 	// particles in the vertex
 	// right up squark
	if(j<k) {
	  addToList( -2*j-1 , -2*k-1 , -2000002-2*i );
	  addToList(  2*j+1 ,  2*k+1 , +2000002+2*i );
	}
	if(i==2) {
	  addToList( -2*j-1 , -2*k-1 , -1000002-2*i );
	  addToList(  2*j+1 ,  2*k+1 , +1000002+2*i );
	}
 	// right down squark
	addToList( -2*i-2 , -2*j-1 , -2000001-2*k );
	addToList(  2*i+2 ,  2*j+1 , +2000001+2*k );
	if(k==2) {
	  addToList( -2*i-2 , -2*j-1 , -1000001-2*k );
	  addToList(  2*i+2 ,  2*j+1 , +1000001+2*k );
 	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVUDDVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix(1.);
  // left up squark
  if(abs(islep) == ParticleID::SUSY_u_R  || 
     abs(islep) == ParticleID::SUSY_c_R || 
     abs(islep) == ParticleID::SUSY_t_1 || 
     abs(islep) == ParticleID::SUSY_t_2) {
    if(abs(islep)>2000000) {
      i = (abs(islep)-2000002)/2;
      if(i==2) mix = (*stop_)(1,1); 
    }
    else {
      i = 2;
      mix = (*stop_)(0,1);
      assert(abs(islep)==ParticleID::SUSY_t_1);
    }
    j = (abs(part1->id())- 1)/2;
    k = (abs(part2->id())- 1)/2;
  }
  else if(abs(islep) == ParticleID::SUSY_d_R || 
	  abs(islep) == ParticleID::SUSY_s_R || 
	  abs(islep) == ParticleID::SUSY_b_1 || 
	  abs(islep) == ParticleID::SUSY_b_2) {
    if(abs(islep)>2000000) {
      k = (abs(islep)-2000001)/2;
      if(k==2) mix = (*sbot_)(1,1);
    }
    else {
      k = 2;
      mix = (*sbot_)(0,1);
      assert(abs(islep)==ParticleID::SUSY_b_1);
    }
    if(abs(part1->id())%2==0) {
      i = (abs(part1->id())- 2)/2;
      j = (abs(part2->id())- 1)/2;
    }
    else {
      i = (abs(part2->id())- 2)/2;
      j = (abs(part1->id())- 1)/2;
      // mix *= -1.;
    }
  }
  else 
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  if(islep>0) {
    left (1.);
    right(0.);
  }
  else {
    left (0.);
    right(1.);
  }
  norm(mix*lambda_[i][j][k]);
}
