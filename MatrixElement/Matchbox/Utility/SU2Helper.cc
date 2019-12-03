// -*- C++ -*-
//
// SU2Helper.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "SU2Helper.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/Repository.h"

using namespace Herwig;

tcPDPtr SU2Helper::SU2CC(tcPDPtr p, int familyShift) { 
  if ( !isInSU2Doublet(p) )
    return tcPDPtr();
  PID idUp = p->id() < 0 ? p->id() + 1 - 2*familyShift : p->id() - 1 + 2*familyShift;
  PID idDown = p->id() < 0 ? p->id() - 1 - 2*familyShift : p->id() + 1 + 2*familyShift;
  if ( !CurrentGenerator::isVoid() )
    return isSU2Up(p) ? 
      CurrentGenerator::current().getParticleData(idUp) :
      CurrentGenerator::current().getParticleData(idDown); 
  return isSU2Up(p) ? 
    Repository::defaultParticle(idUp) :
    Repository::defaultParticle(idDown); 
}

int SU2Helper::family(tcPDPtr p) {

  long id = abs(p->id());

  if ( id > 20 )
    return 0;

  if ( id > 10 )
    id = id - 10;

  if ( id % 2 != 0 )
    id = id + 1;

  return id/2;

}
