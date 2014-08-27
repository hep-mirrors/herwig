// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTJetsAnalysis class.
//

#include "TTJetsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

TTJetsAnalysis::TTJetsAnalysis() {}

TTJetsAnalysis::~TTJetsAnalysis() {}

IBPtr TTJetsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr TTJetsAnalysis::fullclone() const {
  return new_ptr(*this);
}

void TTJetsAnalysis::reconstructHardObjects(ParticleVector& parts) {
  ParticleVector::iterator t = parts.begin();
  for ( ; t != parts.end(); ++t ) {
    if ( (**t).id() == ParticleID::t )
      break;
  }
  if ( t == parts.end() )
    throw Exception() << "No ttbar pair found in TTJetsAnalysis"
		      << Exception::abortnow;
  hardObjectMomentum("T") = (**t).momentum();
  parts.erase(t);
  ParticleVector::iterator tbar = parts.begin();
  for ( ; tbar != parts.end(); ++tbar ) {
    if ( (**tbar).id() == ParticleID::tbar )
      break;
  }
  if ( tbar == parts.end() )
    throw Exception() << "No ttbar pair found in TTJetsAnalysis"
		      << Exception::abortnow;
  hardObjectMomentum("Tbar") = (**tbar).momentum();
  parts.erase(tbar);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TTJetsAnalysis::persistentOutput(PersistentOStream &) const {}

void TTJetsAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TTJetsAnalysis,Herwig::JetsPlusAnalysis>
  describeHerwigTTJetsAnalysis("Herwig::TTJetsAnalysis", "JetCuts.so HwJetsAnalysis.so");

void TTJetsAnalysis::Init() {

  static ClassDocumentation<TTJetsAnalysis> documentation
    ("There is no documentation for the TTJetsAnalysis class");

}

