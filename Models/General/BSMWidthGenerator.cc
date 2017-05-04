// -*- C++ -*-
//
// BSMWidthGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMWidthGenerator class.
//

#include "BSMWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.h"

using namespace Herwig;

IBPtr BSMWidthGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr BSMWidthGenerator::fullclone() const {
  return new_ptr(*this);
}

void BSMWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << theModes;
}

void BSMWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theModes;
}

ClassDescription<BSMWidthGenerator> BSMWidthGenerator::initBSMWidthGenerator;
// Definition of the static class description member.

void BSMWidthGenerator::Init() {

  static ClassDocumentation<BSMWidthGenerator> documentation
    ("A width generator for BSM particles.");

}

void BSMWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, 
				  unsigned int) {
  tcGeneralTwoBodyDecayerPtr dec = 
    dynamic_ptr_cast<tcGeneralTwoBodyDecayerPtr>(decayer);
  theModes.push_back(make_pair(mode, dec));
}

Energy BSMWidthGenerator::partial2BodyWidth(int iloc, Energy m0,
					    Energy m1, Energy m2) const {
  if( m0 < (m1 + m2) ) return ZERO;
  //need pointers to particles involved
  tcDMPtr dm = theModes[iloc].first;
  ParticleMSet::const_iterator pit = dm->products().begin();
  tcPDPtr parta = *pit;
  ++pit;
  tcPDPtr partb = *pit;
  int dummya(0);
  double dummyb(0.);
  if( !theModes[iloc].second->twoBodyMEcode(*dm, dummya, dummyb) ) 
    swap(parta, partb);
  return theModes[iloc].second->partialWidth(make_pair(dm->parent(), m0), 
					     make_pair(parta, m1),
					     make_pair(partb, m2));
}
