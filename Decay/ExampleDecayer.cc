// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExampleDecayer class.
//

#include "ExampleDecayer.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
// #include "Pythia7/Interface/Reference.h" 
#include "Pythia7/EventRecord/Particle.h"
#include "Pythia7/PDT/PDT.h"   
#include "Pythia7/PDT/EnumParticles.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "DecayConfig.h"


using namespace Herwig;


ExampleDecayer::~ExampleDecayer() {}


void ExampleDecayer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void ExampleDecayer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<ExampleDecayer> ExampleDecayer::initExampleDecayer;
// Definition of the static class description member.

void ExampleDecayer::Init() {

  static ClassDocumentation<ExampleDecayer> documentation
    ("This is a dummy class just to provide an example.");

}


bool ExampleDecayer::accept(const DecayMode & dm) const {
  //*** WRITE THE CODE *** 
  return true;                   // fake: just to compile!
}


ParticleVector ExampleDecayer::decay(const DecayMode & dm, const Particle & part) const {
  //*** WRITE THE CODE *** 
  return ParticleVector();       // fake: just to compile!
}
