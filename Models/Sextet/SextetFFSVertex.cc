// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetFFSVertex class.
//

#include "SextetFFSVertex.h"
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

IBPtr SextetFFSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetFFSVertex::fullclone() const {
  return new_ptr(*this);
}

void SextetFFSVertex::persistentOutput(PersistentOStream & os) const {
  os << g1L_ << g1R_ << g1pR_ << g1ppR_ << g3L_;
}

void SextetFFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> g1L_ >> g1R_ >> g1pR_ >> g1ppR_ >> g3L_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetFFSVertex,Helicity::FFSVertex>
describeSextetFFSVertex("Herwig::SextetFFSVertex", "HwSextetModel.so");

void SextetFFSVertex::Init() {

  static ClassDocumentation<SextetFFSVertex> documentation
    ("The SextetFFSVertex class implements the coupling of two "
     "fermions to a scalar sextet particle.");

}

void SextetFFSVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // extract the couplings
  g1L_     = model->g1L();
  g1R_     = model->g1R();
  g1pR_   = model->g1pR();
  g1ppR_ = model->g1ppR();
  g3L_     = model->g3L();
  // add the enabled particles
  if(model->ScalarSingletY43Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      if(g1ppR_[ix]!=0.) {
        addToList(  iu,  iu, ParticleID::ScalarDQSingletY43bar);
        addToList( -iu, -iu, ParticleID::ScalarDQSingletY43);
      }
    }
  }
  if(model->ScalarSingletY13Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g1L_[ix]!=0. || g1R_[ix]!=0.) {
        addToList(  id,  iu, ParticleID::ScalarDQSingletY13bar);
        addToList( -id, -iu, ParticleID::ScalarDQSingletY13);
      }
    }
  }
  if(model->ScalarSingletY23Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long id = 2*ix + 1;
      if(g1pR_[ix]!=0. ) {
        addToList( id, id,ParticleID::ScalarDQSingletY23bar);
        addToList(-id,-id,ParticleID::ScalarDQSingletY23);
      }
    }
  }
  if(model->ScalarTripletY13Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g3L_[ix]!=0. ) {
        addToList(  iu,  iu, ParticleID::ScalarDQTripletPbar);
        addToList( -iu, -iu, ParticleID::ScalarDQTripletP);
        addToList(  iu,  id, ParticleID::ScalarDQTriplet0bar);
        addToList( -iu, -id, ParticleID::ScalarDQTriplet0);
        addToList(  id,  id, ParticleID::ScalarDQTripletMbar);
        addToList( -id, -id, ParticleID::ScalarDQTripletM);
      }
    }
  }
  Helicity::FFSVertex::doinit();
}


void SextetFFSVertex::setCoupling(Energy2,tcPDPtr part1,
				  tcPDPtr part2,tcPDPtr part3) {
  long q1ID=(abs(part1->id())), q2ID=(abs(part2->id())),
    sDQID=(abs(part3->id()));
  //check scalar diquark
  assert( sDQID == ParticleID::ScalarDQSingletY43 || 
          sDQID == ParticleID::ScalarDQSingletY13 ||
          sDQID == ParticleID::ScalarDQSingletY23 ||
          sDQID == ParticleID::ScalarDQTripletP ||
          sDQID == ParticleID::ScalarDQTriplet0 ||
          sDQID == ParticleID::ScalarDQTripletM); 
  //check quarks
  assert(!(q1ID>6) && !(q2ID>6));
  bool part1Up = (q1ID==2 || q1ID==4 || q1ID==6);
  bool part2Up = (q2ID==2 || q2ID==4 || q2ID==6);
  Complex cRight, cLeft, prefactor(1.);
  if(sDQID==ParticleID::ScalarDQSingletY43){
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cRight=Complex(g1ppR_[0]);
    else if(q1ID==4)
      cRight=Complex(g1ppR_[1]);
    else
      cRight=Complex(g1ppR_[2]);
    cLeft=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQSingletY13){
    //should be one up one down type
    assert((part1Up && !part2Up) || (!part1Up && part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;   
    if(upType==2){
      cRight=Complex(g1R_[0]);
      cLeft=Complex(2.*g1L_[0]);
    }
    else if(upType==4){
      cRight=Complex(g1R_[1]);
      cLeft=Complex(2.*g1L_[1]);
    }
    else
      cRight=Complex(g1R_[2]);{
      cLeft=Complex(2.*g1L_[2]);
    }
  }
  if(sDQID==ParticleID::ScalarDQSingletY23){
    //should both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cRight=Complex(g1pR_[0]);
    else if(q1ID==3)
      cRight=Complex(g1pR_[1]);
    else
      cRight=Complex(g1pR_[2]);
    cLeft=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTripletP){
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cLeft=Complex(g3L_[0]);
    else if(q1ID==4)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);  
    cRight=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTriplet0){
    //should both one up and down type
    assert((part1Up && !part2Up) || (!part1Up && part2Up));
    //possibly doesn't couple
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;   
    if(upType==2)
      cLeft=Complex(g3L_[0]);
    else if(upType==4)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);  
    cRight=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTripletM){
    //should one both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cLeft=Complex(g3L_[0]);
    else if(q1ID==3)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);
    cRight=Complex(0.);
    prefactor=Complex(-1.); 
  }
  left(cLeft);
  right(cRight);
  norm(prefactor);
}
