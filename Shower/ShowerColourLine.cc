// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerColourLine class.
//

#include "ShowerColourLine.h"
// #include "Pythia7/EventRecord/Particle.h"
#include "Pythia7/Config/algorithm.h"

using namespace Herwig;

ShowerColourLine::~ShowerColourLine() {}

void ShowerColourLine::addAntiColoured(tShoParPtr p) {
  theAntiColoured.push_back(p);
  p->setAntiColourLine(this);
}

void ShowerColourLine::addColoured(tShoParPtr p, bool anti) {
  if ( anti ) addAntiColoured(p);
  else {
    theColoured.push_back(p);
    p->setColourLine(this);
  }
}

void ShowerColourLine::removeAntiColoured(tShoParPtr p) {
  theAntiColoured.erase(find(range(theAntiColoured), p));
  p->setAntiColourLine(ShoColinePtr());
}

void ShowerColourLine::removeColoured(tShoParPtr p, bool anti) {
  if ( anti ) removeAntiColoured(p);
  else {
    theColoured.erase(find(range(theColoured), p));
    p->setColourLine(ShoColinePtr());
  }
}

