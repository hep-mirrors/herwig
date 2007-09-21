// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ToIncomingCMS class.
//

#include "ToIncomingCMS.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ClusteringParticle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ToIncomingCMS.tcc"
#endif


using namespace Herwig;

ToIncomingCMS::~ToIncomingCMS() {}

NoPIOClassDescription<ToIncomingCMS> ToIncomingCMS::initToIncomingCMS;
// Definition of the static class description member.

void ToIncomingCMS::Init() {

  static ClassDocumentation<ToIncomingCMS> documentation
    ("Boosts the particles to the CMS of the new incoming partons after a "
     "clustering involving initial state particles has been performed.");

}

void ToIncomingCMS::initialize (const vector<tClusteringParticlePtr>& event) {
  // get the initial state particles
  vector<tClusteringParticlePtr> initial;
  for (vector<tClusteringParticlePtr>::const_iterator p = event.begin();
       p != event.end(); ++p)
    if ((**p).pData().partonId.state == ClusteringParticleState::initial)
      initial.push_back(*p);
  if(initial.size() != 2) {
    throw Exception() << "HwCKKW : ToIncomingCMS::initialize : Didn't get two incoming particles."
		      << Exception::eventerror;
  }
  _boost = (initial[0]->momentum() + initial[1]->momentum()).findBoostToCM();
}

void ToIncomingCMS::doTransform (tClusteringParticlePtr p) {
  p->momentum().transform(_boost);
}
