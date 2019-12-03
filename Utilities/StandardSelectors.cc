// -*- C++ -*-
//
// StandardSelectors.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions for the Herwig StandardSelectors
//
#include "StandardSelectors.h"

using namespace Herwig;

bool WeakBHadronSelector::check(const Particle & p) const {
  long id = abs(p.id());
  if (!( ( id > 510 && id < 532 ) || ( id > 5121 && id < 5555 ) ) )
    return false; 
  if ( p.children().size()==1 && abs(p.children()[0]->id()) == id )
    return false;
  switch(id)
    {
    case  511: // B0
    case 521: // B+
    case 531: // B_s0
    case 5122: // Lambda_b0
    case 5132: // Xi_b-
    case 5142: // Xi_bc0
    case 5232: // Xi_b0
    case 5242: // Xi_bc+
    case 5332: // Omega_b
    case 5342: // Omega_bc0
    case 5412: // Xi'_bc0
    case 5414: // Xi*_bc0
    case 5422: // Xi'_bc+
    case 5424: // Xi*_bc+
    case 5432: // Omega'_bc0
    case 5434: // Omega*_bc0
    case 5442: // Omega_bcc+
    case 5444: // Omega*_bcc+
    case 5512: // Xi_bb-
    case 5514: // Xi*_bb-
    case 5522: // Xi_bb0
    case 5524: // Xi*_bb0
    case 5532: // Omega_bb-
    case 5534: // Omega*_bb-
    case 5542: // Omega_bbc0
    case 5544: // Omega*_bbc0
    case 5554: // Omega*_bbb-
      return true;
      break;
    default:
      return false;
    }
}
