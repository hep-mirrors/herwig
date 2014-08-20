// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZJetsAnalysis class.
//

#include "ZJetsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ZJetsAnalysis::ZJetsAnalysis() {}

ZJetsAnalysis::~ZJetsAnalysis() {}

IBPtr ZJetsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr ZJetsAnalysis::fullclone() const {
  return new_ptr(*this);
}

void ZJetsAnalysis::reconstructHardObjects(ParticleVector& parts) {
  ParticleVector::iterator eplus = parts.begin();
  for ( ; eplus != parts.end(); ++eplus ) {
    if ( (**eplus).id() == ParticleID::eplus )
      break;
  }
  if ( eplus == parts.end() )
    throw Exception() << "No e+e- pair found in ZJetsAnalysis"
		      << Exception::abortnow;
  LorentzMomentum peplus = (**eplus).momentum();
  parts.erase(eplus);
  ParticleVector::iterator eminus = parts.begin();
  for ( ; eminus != parts.end(); ++eminus ) {
    if ( (**eminus).id() == ParticleID::eminus )
      break;
  }
  if ( eminus == parts.end() )
    throw Exception() << "No e+e- pair found in ZJetsAnalysis"
		      << Exception::abortnow;
  LorentzMomentum peminus = (**eminus).momentum();
  parts.erase(eminus);
  hardObjectMomentum("Z") = peplus + peminus;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ZJetsAnalysis::persistentOutput(PersistentOStream &) const {}

void ZJetsAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ZJetsAnalysis,Herwig::JetsPlusAnalysis>
  describeHerwigZJetsAnalysis("Herwig::ZJetsAnalysis", "JetCuts.so HwJetsAnalysis.so");

void ZJetsAnalysis::Init() {

  static ClassDocumentation<ZJetsAnalysis> documentation
    ("There is no documentation for the ZJetsAnalysis class");

}

