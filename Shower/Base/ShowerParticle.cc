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
  for(map<ShowerPartnerType::Type ,pair<Energy,Energy> >::iterator
	it=evolutionScales_.begin();it!=evolutionScales_.end();++it) {
    if(it->first==type || (it->first!=ShowerPartnerType::QED &&
			   type     !=ShowerPartnerType::QED)) {
      it->second = make_pair(scale,scale);
    }
    else
      assert(false);
  }
}
