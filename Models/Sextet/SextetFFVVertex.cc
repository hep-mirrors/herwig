// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetFFVVertex class.
//

#include "SextetFFVVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetFFVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetFFVVertex::fullclone() const {
  return new_ptr(*this);
}

void SextetFFVVertex::persistentOutput(PersistentOStream & os) const {
  os << g2_ << g2p_ ;
}

void SextetFFVVertex::persistentInput(PersistentIStream & is, int) {
  is >> g2_ >> g2p_ ;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetFFVVertex,Helicity::FFVVertex>
describeSextetFFVVertex("Herwig::SextetFFVVertex", "HwSextetModel.so");

void SextetFFVVertex::Init() {

  static ClassDocumentation<SextetFFVVertex> documentation
    ("The SextetFFVVertex class implements the coupling of two "
     "fermions to a scalar sextet particle.");

}

void SextetFFVVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // extract the couplings
  g2_  = model->g2();
  g2p_= model->g2p();
  // add the enabled particles
  if(model->VectorDoubletY16Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g2_[ix]!=0.) {
        addToList(  -id, -iu, ParticleID::VectorDQY16P);
        addToList(  id,  iu, ParticleID::VectorDQY16Pbar);
        addToList( -id, -id, ParticleID::VectorDQY16M);
        addToList(  id,  id, ParticleID::VectorDQY16Mbar);
      }
    }
  }
  if(model->VectorDoubletY56Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g2p_[ix]!=0.) {
        addToList( -iu, -iu, ParticleID::VectorDQY56P);
        addToList(  iu,  iu, ParticleID::VectorDQY56Pbar);
        addToList( -id, -iu, ParticleID::VectorDQY56M);
        addToList(  id,  iu, ParticleID::VectorDQY56Mbar);
      }
    }
  }
  Helicity::FFVVertex::doinit();
}


void SextetFFVVertex::setCoupling(Energy2, tcPDPtr part1,
				  tcPDPtr part2,tcPDPtr part3) {
  long q1ID=(abs(part1->id())), q2ID=(abs(part2->id())),
    vDQID=(abs(part3->id()));
  //check scalar diquark
  assert( vDQID == ParticleID::VectorDQY16P || 
          vDQID == ParticleID::VectorDQY16M ||
          vDQID == ParticleID::VectorDQY56P ||
          vDQID == ParticleID::VectorDQY56M);
  //check quarks
  assert(!(q1ID>6) && !(q2ID>6));
  bool part1Up = (q1ID==2 || q1ID==4 || q1ID==6) ? true : false;
#ifndef NDEBUG
  bool part2Up = (q2ID==2 || q2ID==4 || q2ID==6) ? true : false;
#endif
  Complex cRight(1.,1.), cLeft(1.,1.), prefactor(1.,0.);

  if(vDQID==ParticleID::VectorDQY16P){
    //should be one up and down type
    assert((!part1Up && part2Up) || (part1Up && !part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;
    if(upType==2)
      cRight=Complex(g2_[0]);
    else if(upType==4)
      cRight=Complex(g2_[1]);
    else
      cRight=Complex(g2_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY16M) {
    //should be both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cRight=Complex(g2_[0]);
    else if(q1ID==2)
      cRight=Complex(g2_[1]);
    else
      cRight=Complex(g2_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY56P) {
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cRight=Complex(g2p_[0]);
    else if(q1ID==4)
      cRight=Complex(g2p_[1]);
    else
      cRight=Complex(g2p_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY56M){
    //should be one up and down type
    assert((!part1Up && part2Up) || (part1Up && !part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;
    if(upType==2)
      cRight=Complex(g2p_[0]);
    else if(upType==4)
      cRight=Complex(g2p_[1]);
    else
      cRight=Complex(g2p_[2]);
    cLeft=Complex(0.);
  }
  left(cLeft);
  right(cRight);
  norm(prefactor);
}



  


