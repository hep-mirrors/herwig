// -*- C++ -*-
//
// BtoSGammaFlatEnergy.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaFlatEnergy class.
//

#include "BtoSGammaFlatEnergy.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NoPIOClassDescription<BtoSGammaFlatEnergy> 
BtoSGammaFlatEnergy::initBtoSGammaFlatEnergy;
// Definition of the static class description member.

void BtoSGammaFlatEnergy::Init() {

  static ClassDocumentation<BtoSGammaFlatEnergy> documentation
    ("The BtoSGammaFlatEnergy class implements a model of the hadronic mass"
     " in B to s gamma decays designed to give a flat energy spectrum for"
     " testing purposes.");

}

Energy BtoSGammaFlatEnergy::hadronicMass(Energy mb,Energy mquark) {
  Energy upper(min(mb,maxMass())),lower(max(minMass(),mquark));
  double rand(UseRandom::rnd());
  return sqrt(upper*upper*rand+(1.-rand)*lower*lower);
}

void BtoSGammaFlatEnergy::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BtoSGammaFlatEnergy " 
		    << name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() 
		    << "\";" << endl;
}
