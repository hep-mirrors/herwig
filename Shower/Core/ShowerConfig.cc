// -*- C++ -*-
//
// ShowerConfig.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerConfig class.
//
#include "ShowerConfig.h"
#include "Core/Base/SudakovFormFactor.h"

using namespace Herwig;

BranchingElement::~BranchingElement() {}

BranchingElement::BranchingElement() {}

BranchingElement::BranchingElement(SudakovPtr sud, IdList part) : sudakov(sud), particles(part) {
  for(unsigned int ix=0;ix<part.size();++ix)
    conjugateParticles.push_back(part[ix]->CC() ? tcPDPtr(part[ix]->CC()) : part[ix]);
}

namespace ThePEG {

/** 
 * Output operator to allow the structure
 */
PersistentOStream & operator << (PersistentOStream & os, 
				 const Herwig::BranchingElement  & x) {
  os << x.sudakov << x.particles << x.conjugateParticles;
  return os;
}

/** 
 * Input operator to allow the structure
 */
PersistentIStream & operator >> (PersistentIStream & is, 
				 Herwig::BranchingElement  & x) {
  is >> x.sudakov >> x.particles >> x.conjugateParticles;
  return is;
}
 
}
