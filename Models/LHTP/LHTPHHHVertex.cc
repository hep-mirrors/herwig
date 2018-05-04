// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPHHHVertex class.
//

#include "LHTPHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHTPHHHVertex::LHTPHHHVertex() : ratio_(ZERO), coupLast_(0.), q2Last_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHTPHHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPHHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(ratio_,GeV);
}

void LHTPHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ratio_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<LHTPHHHVertex,SSSVertex>
describeHerwigLHTPHHHVertex("Herwig::LHTPHHHVertex", "HwLHTPModel.so");

void LHTPHHHVertex::Init() {

  static ClassDocumentation<LHTPHHHVertex> documentation
    ("The LHTPHHHVertex class implements the trilinear Higgs boson"
     " self couplings in the Little Higgs model with T-parity");

}

void LHTPHHHVertex::doinit() {
  addToList(25,25, 25);
  addToList(25,35, 35);
  addToList(25,36, 36);
  addToList(25,37,-37);
  SSSVertex::doinit();
  ratio_ = sqr(getParticleData(ParticleID::h0)->mass())/
    getParticleData(ParticleID::Wplus)->mass();
}

void LHTPHHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr) {
  if(q2!=q2Last_||coupLast_==0.) {
    coupLast_ = weakCoupling(q2)*ratio_*UnitRemoval::InvE;
    q2Last_=q2;
  }
  long id = part2->id()!=ParticleID::h0 ? abs(part2->id()) : abs(part1->id());
  if(id==ParticleID::h0) norm(    -coupLast_);
  else if(id==35)        norm(  3.*coupLast_);
  else if(id==36)        norm(  3.*coupLast_);
  else if(id==37)        norm(1.5*coupLast_);
  else assert(false);
}
