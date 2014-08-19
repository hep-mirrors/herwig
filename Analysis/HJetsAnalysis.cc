// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HJetsAnalysis class.
//

#include "HJetsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HJetsAnalysis::HJetsAnalysis() {}

HJetsAnalysis::~HJetsAnalysis() {}

IBPtr HJetsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr HJetsAnalysis::fullclone() const {
  return new_ptr(*this);
}

void HJetsAnalysis::reconstructHardObjects(ParticleVector& parts) {
  ParticleVector::iterator higgs = parts.begin();
  for ( ; higgs != parts.end(); ++higgs ) {
    if ( (**higgs).id() == ParticleID::h0 )
      break;
  }
  if ( higgs == parts.end() )
    throw Exception() << "No Higgs found in HJetsAnalysis"
		      << Exception::abortnow;
  hardObjectMomentum("h") = (**higgs).momentum();
  parts.erase(higgs);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HJetsAnalysis::persistentOutput(PersistentOStream &) const {}

void HJetsAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HJetsAnalysis,Herwig::JetsPlusAnalysis>
  describeHerwigHJetsAnalysis("Herwig::HJetsAnalysis", "JetCuts.so HwJetsAnalysis.so");

void HJetsAnalysis::Init() {

  static ClassDocumentation<HJetsAnalysis> documentation
    ("There is no documentation for the HJetsAnalysis class");

}

