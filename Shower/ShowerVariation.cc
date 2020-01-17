// -*- C++ -*-
//
// ShowerVariation.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVariation class.
//

#include "ShowerVariation.h"

using namespace Herwig;

string ShowerVariation::fromInFile(const string& in) {
  // pretty simple for the moment, just to try
  // TODO make this better
  istringstream read(in);
  string where;
  read >> renormalizationScaleFactor >> factorizationScaleFactor >> where;
  if ( !read )
    return "something went wrong with: " + in;
  if ( where != "Hard" && where != "All" && where!= "Secondary" )
    return "The specified process for reweighting does not exist.\nOptions are: Hard, Secondary, All.";
  if ( where == "Hard" || where == "All" )
    firstInteraction = true;
  else
    firstInteraction = false;
  if ( where == "Secondary" || where == "All" )
    secondaryInteractions = true;
  else
    secondaryInteractions = false;
  return "";
}

void ShowerVariation::put(PersistentOStream& os) const {
  os << renormalizationScaleFactor << factorizationScaleFactor
     << firstInteraction << secondaryInteractions;
}

void ShowerVariation::get(PersistentIStream& is) {
  is >> renormalizationScaleFactor >> factorizationScaleFactor
     >> firstInteraction >> secondaryInteractions;
}
