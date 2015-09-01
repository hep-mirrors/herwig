// -*- C++ -*-
//
// ShowerParticle.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <ostream>

using namespace Herwig;

PPtr ShowerParticle::clone() const {
  return new_ptr(*this);
}

PPtr ShowerParticle::fullclone() const {
  return new_ptr(*this);
}

ClassDescription<ShowerParticle> ShowerParticle::initShowerParticle;
// Definition of the static class description member.

void ShowerParticle::vetoEmission(ShowerPartnerType::Type, Energy scale) {
  scales_.QED         = min(scale,scales_.QED        );
  scales_.QED_noAO    = min(scale,scales_.QED_noAO   );
  scales_.QCD_c       = min(scale,scales_.QCD_c      );
  scales_.QCD_c_noAO  = min(scale,scales_.QCD_c_noAO );
  scales_.QCD_ac      = min(scale,scales_.QCD_ac     );
  scales_.QCD_ac_noAO = min(scale,scales_.QCD_ac_noAO);
}

void ShowerParticle::addPartner(EvolutionPartner in ) {
  partners_.push_back(in); 
}
