// -*- C++ -*-
//
// TreePhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TreePhasespace class.
//

#include "TreePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;

TreePhasespace::TreePhasespace() 
  : x0(0.01), xc(1e-4) {
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
}

TreePhasespace::~TreePhasespace() {}

IBPtr TreePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr TreePhasespace::fullclone() const {
  return new_ptr(*this);
}

void TreePhasespace::prepare(tStdXCombPtr xco, bool) {

  theLastXComb = xco;

  lastChannelsIterator = channelMap().find(lastXCombPtr());

  if ( lastChannelsIterator == channelMap().end() ) {
    map<Ptr<Tree2toNDiagram>::ptr,PhasespaceTree> channels;
    for ( StandardXComb::DiagramVector::const_iterator d =
	    lastXComb().diagrams().begin(); d != lastXComb().diagrams().end(); ++d ) {
      PhasespaceTree tree;
      Ptr<Tree2toNDiagram>::ptr diag =
	dynamic_ptr_cast<Ptr<Tree2toNDiagram>::ptr>(*d);
      tree.setup(*diag);
      channels[diag] = tree;
    }
    channelMap()[lastXCombPtr()] = channels;
    lastChannelsIterator = channelMap().find(lastXCombPtr());
  }

  lastPhasespaceInfo.sHat = lastXComb().lastSHat();
  lastPhasespaceInfo.sqrtSHat = sqrt(lastXComb().lastSHat());
  lastPhasespaceInfo.weight = 1.;

}

double TreePhasespace::generateKinematics(const double* random,
					  vector<Lorentz5Momentum>& momenta) {

  size_t nchannels = lastXComb().diagrams().size();
  map<Ptr<Tree2toNDiagram>::ptr,PhasespaceHelpers::PhasespaceTree>::iterator ds =
    lastChannels().begin();
  advance(ds,(size_t)(random[0]*nchannels));
  Ptr<Tree2toNDiagram>::ptr channel = ds->first;
  ++random;

  lastPhasespaceInfo.rnd.numbers = random;
  lastPhasespaceInfo.rnd.nRnd = 3*momenta.size() - 10;

  cPDVector::const_iterator pd = mePartonData().begin();
  vector<Lorentz5Momentum>::iterator p = momenta.begin();
  for ( ; pd != mePartonData().end(); ++pd, ++p )
    p->setMass((**pd).mass());

  try {
    lastChannels()[channel].generateKinematics(lastPhasespaceInfo,momenta);
  } catch (Veto) {
    return 0.;
  }

  fillDiagramWeights();

  double sum = 0.;
  for ( map<Ptr<Tree2toNDiagram>::ptr,PhasespaceHelpers::PhasespaceTree>::const_iterator d
	  = lastChannels().begin(); d != lastChannels().end(); ++d )
    sum += diagramWeight(*(d->first));

  double piWeight = pow(2.*Constants::pi,(double)(3*(momenta.size()-2)-4));

  return nchannels*lastPhasespaceInfo.weight*diagramWeight(*channel)/(sum*piWeight);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void TreePhasespace::doinit() {
  MatchboxPhasespace::doinit();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
}

void TreePhasespace::doinitrun() {
  MatchboxPhasespace::doinitrun();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
}

void TreePhasespace::persistentOutput(PersistentOStream & os) const {
  os << theChannelMap << x0 << xc;
}

void TreePhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theChannelMap >> x0 >> xc;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TreePhasespace,MatchboxPhasespace>
  describeHerwigTreePhasespace("Herwig::TreePhasespace", "HwMatchbox.so");

void TreePhasespace::Init() {

  static ClassDocumentation<TreePhasespace> documentation
    ("TreePhasespace is a multichannel phasespace generator "
     "adapting to singularity structures as determined from the matrix "
     "elements diagrams.");


  static Reference<TreePhasespace,TreePhasespaceChannels> interfaceChannelMap
    ("ChannelMap",
     "Set the object storing the channels.",
     &TreePhasespace::theChannelMap, false, false, true, false, false);


  static Parameter<TreePhasespace,double> interfaceX0
    ("X0",
     "Set the cut below which flat virtuality sampling is imposed.",
     &TreePhasespace::x0, 0.01, 0.0, 0,
     false, false, Interface::lowerlim);


  static Parameter<TreePhasespace,double> interfaceXC
    ("XC",
     "Set the cut below which no virtualities are generated.",
     &TreePhasespace::xc, 1e-4, 0.0, 0,
     false, false, Interface::lowerlim);


}

