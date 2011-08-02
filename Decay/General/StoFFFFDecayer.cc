// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoFFFFDecayer class.
//

#include "StoFFFFDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"

using namespace Herwig;

StoFFFFDecayer::StoFFFFDecayer() {}

IBPtr StoFFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoFFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoFFFFDecayer::persistentOutput(PersistentOStream & os) const {
}

void StoFFFFDecayer::persistentInput(PersistentIStream & is, int) {
}

DescribeClass<StoFFFFDecayer,GeneralFourBodyDecayer>
describeStoFFFFDecayer("Herwig::StoFFFFDecayer", "Herwig.so");

void StoFFFFDecayer::Init() {

  static ClassDocumentation<StoFFFFDecayer> documentation
    ("The StoFFFFDecayer class performs the 4-body decays of scalar particles in BSM models");

}

double StoFFFFDecayer::me2(const int ichan, const Particle & part,
			   const ParticleVector & decay, MEOption meopt) const {
  assert(false);
  return 0.;
}
