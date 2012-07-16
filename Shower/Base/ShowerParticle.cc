// -*- C++ -*-
//
// ShowerParticle.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PPtr ShowerParticle::clone() const {
  return new_ptr(*this);
}

PPtr ShowerParticle::fullclone() const {
  return new_ptr(*this);
}

ClassDescription<ShowerParticle> ShowerParticle::initShowerParticle;
// Definition of the static class description member.

void ShowerParticle::vetoEmission(ShowerPartnerType::Type type, Energy scale) {
  // for(map<ShowerPartnerType::Type ,pair<Energy,Energy> >::iterator
  // 	it=evolutionScales_.begin();it!=evolutionScales_.end();++it) {
  //   if(it->first==type) {
  //     it->second = make_pair(scale,scale);
  //   }
  //   else {
  //     if(it->second.first >scale) it->second.first  = scale;
  //     if(it->second.second>scale) it->second.second = scale;
  //   }
  // }
  assert(false);
}

void ShowerParticle::addPartner(EvolutionPartner in ) {
  partners_.push_back(in); 
}
