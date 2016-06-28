// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ElectroWeakReweighter class.
//

#include "ElectroWeakReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ElectroWeakReweighter::ElectroWeakReweighter() {}

ElectroWeakReweighter::~ElectroWeakReweighter() {}

IBPtr ElectroWeakReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr ElectroWeakReweighter::fullclone() const {
  return new_ptr(*this);
}

void ElectroWeakReweighter::persistentOutput(PersistentOStream & os) const {
  os << EWCouplings_;
}

void ElectroWeakReweighter::persistentInput(PersistentIStream & is, int) {
  is >> EWCouplings_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ElectroWeakReweighter,ReweightBase>
  describeHerwigElectroWeakReweighter("Herwig::ElectroWeakReweighter", "HwMEEW.so");

void ElectroWeakReweighter::Init() {

  static ClassDocumentation<ElectroWeakReweighter> documentation
    ("There is no documentation for the ElectroWeakReweighter class");

  static Reference<ElectroWeakReweighter,EWCouplings> interfaceEWCouplings
    ("EWCouplings",
     "The object to calculate the electroweak couplings",
     &ElectroWeakReweighter::EWCouplings_, false, false, true, false, false);

}

double ElectroWeakReweighter::weight() const {
  EWCouplings_->initialize();
  staticEWCouplings_ = EWCouplings_;
  cerr << "aEM\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.1) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->aEM(scale) << "\n";
  }
  cerr << "aS\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->aS(scale) << "\n";
  }
  cerr << "y_t\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->y_t(scale) << "\n";
  }
  cerr << "lambda\n";
  for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->lambda(scale) << "\n";
  }
  cerr << "vev\n";
  for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->vev(scale)/GeV << "\n";
  }


  cerr << "testing got here ???\n";
  cerr <<  subProcess() << "\n";
  cerr << *subProcess() << "\n";
  cerr << subProcess()->outgoing()[0] << *subProcess()->outgoing()[0] << "\n";
  cerr << subProcess()->outgoing()[0]->spinInfo() << "\n";
  cerr << subProcess()->outgoing()[0]->spinInfo()->productionVertex() << "\n";
  assert(false);
  staticEWCouplings_ = tEWCouplingsPtr();
}
