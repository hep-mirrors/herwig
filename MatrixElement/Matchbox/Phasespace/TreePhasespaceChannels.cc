// -*- C++ -*-
//
// TreePhasespaceChannels.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TreePhasespaceChannels class.
//

#include "TreePhasespaceChannels.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;

TreePhasespaceChannels::TreePhasespaceChannels() {}

TreePhasespaceChannels::~TreePhasespaceChannels() {}

IBPtr TreePhasespaceChannels::clone() const {
  return new_ptr(*this);
}

IBPtr TreePhasespaceChannels::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TreePhasespaceChannels::persistentOutput(PersistentOStream & os) const {
  os << theChannelMap.size();
  for ( map<tStdXCombPtr,map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> > >::const_iterator k =
          theChannelMap.begin(); k != theChannelMap.end(); ++k ) {
    os << k->first << k->second.size();
    for ( map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> >::const_iterator l = 
            k->second.begin(); l != k->second.end(); ++l ) {
      os << l->first;
      l->second.first.put(os);
      l->second.second.put(os);
    }
  }
}

void TreePhasespaceChannels::persistentInput(PersistentIStream & is, int) {
  size_t nk; is >> nk;
  for ( size_t k = 0; k < nk; ++k ) {
    tStdXCombPtr xc; is >> xc;
    size_t nl; is >> nl;
    map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> > cm;
    for ( size_t l = 0; l < nl; ++l ) {
      Ptr<Tree2toNDiagram>::ptr ci; is >> ci;
      pair<PhasespaceTree,PhasespaceTree> cp; cp.first.get(is); cp.second.get(is);
      cm[ci] = cp;
    }
    theChannelMap[xc] = cm;
  }
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TreePhasespaceChannels,HandlerBase>
  describeHerwigTreePhasespaceChannels("Herwig::TreePhasespaceChannels", "Herwig.so");

void TreePhasespaceChannels::Init() {

  static ClassDocumentation<TreePhasespaceChannels> documentation
    ("Store channels for the tree phasespace.");

}

