// -*- C++ -*-
//
// TreePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TreePhasespace class.
//

#include <sstream> 
#include <string> 
#include "TreePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;

TreePhasespace::TreePhasespace() 
  : x0(0.01), xc(1e-4), M0(ZERO), Mc(ZERO) {
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
  theIncludeMirrored = true;
}

TreePhasespace::~TreePhasespace() {}

IBPtr TreePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr TreePhasespace::fullclone() const {
  return new_ptr(*this);
}

void TreePhasespace::setXComb(tStdXCombPtr xco) {

  MatchboxPhasespace::setXComb(xco);

  lastChannelsIterator = channelMap().find(lastXCombPtr());

  if ( lastChannelsIterator == channelMap().end() ) {
    map<Ptr<Tree2toNDiagram>::ptr,pair<PhasespaceTree, PhasespaceTree> > channels;
    for ( auto const & d : lastXComb().diagrams()) {
      PhasespaceTree tree;
      Ptr<Tree2toNDiagram>::ptr diag =
	dynamic_ptr_cast<Ptr<Tree2toNDiagram>::ptr>(d);
      tree.setup(*diag);
      PhasespaceTree treeMirror;
      treeMirror.setupMirrored(*diag, diag->nSpace() - 1);
      channels[diag] = make_pair(tree,treeMirror);
    }
    channelMap()[lastXCombPtr()] = channels;
    lastChannelsIterator = channelMap().find(lastXCombPtr());
  }

}

double TreePhasespace::generateTwoToNKinematics(const double* random,
						vector<Lorentz5Momentum>& momenta) {

  lastPhasespaceInfo.sHat = lastXComb().lastSHat();
  lastPhasespaceInfo.sqrtSHat = sqrt(lastXComb().lastSHat());
  lastPhasespaceInfo.weight = 1.;

  size_t nchannels = lastXComb().diagrams().size();
  bool doMirror = (UseRandom::rnd() < 0.5) && theIncludeMirrored;
  map<Ptr<Tree2toNDiagram>::ptr,
      pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> >::iterator ds =
    lastChannels().begin();

  size_t i = (size_t)(random[0]*nchannels);
  advance(ds,i);

  Ptr<Tree2toNDiagram>::ptr channel = ds->first;
  ++random;
  
  lastPhasespaceInfo.rnd.numbers = random;
  lastPhasespaceInfo.rnd.nRnd = 3*momenta.size() - 10;
    
  try {
    if ( !doMirror )
      lastChannels()[channel].first.generateKinematics(lastPhasespaceInfo,momenta);
    else 
      lastChannels()[channel].second.generateKinematics(lastPhasespaceInfo,momenta);
  } catch (Veto) {
    return 0.;
  }
    
  if ( !matchConstraints(momenta) )
    return 0.;

  double flatCut = x0;
  if ( M0 != ZERO )
    flatCut = M0/sqrt(lastSHat());

  fillDiagramWeights(flatCut);

  double sum = 0.;
  for ( auto const & d : lastChannels())
    sum += diagramWeight(*(d.first));

  double piWeight = pow(2.*Constants::pi,(double)(3*(momenta.size()-2)-4));

  for ( auto & k : momenta )
    k.rescaleRho();

  return nchannels*lastPhasespaceInfo.weight*diagramWeight(*channel)/(sum*piWeight);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void TreePhasespace::doinit() {
  MatchboxPhasespace::doinit();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
}

void TreePhasespace::doinitrun() {
  MatchboxPhasespace::doinitrun();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
}

void TreePhasespace::persistentOutput(PersistentOStream & os) const {
  os << theChannelMap << x0 << xc 
     << ounit(M0,GeV) << ounit(Mc,GeV)
     << theIncludeMirrored
     << theLastXComb;
}

void TreePhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theChannelMap >> x0 >> xc 
     >> iunit(M0,GeV) >> iunit(Mc,GeV)
     >> theIncludeMirrored
     >> theLastXComb;
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
  lastChannelsIterator = channelMap().find(lastXCombPtr());
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TreePhasespace,MatchboxPhasespace>
  describeHerwigTreePhasespace("Herwig::TreePhasespace", "Herwig.so");

void TreePhasespace::Init() {

  static ClassDocumentation<TreePhasespace> documentation
    ("TreePhasespace is a multi-channel phase space generator "
     "adapting to singularity structures as determined from the matrix "
     "elements diagrams.");


  static Reference<TreePhasespace,TreePhasespaceChannels> interfaceChannelMap
    ("ChannelMap",
     "Set the object storing the channels.",
     &TreePhasespace::theChannelMap, false, false, true, false, false);
  interfaceChannelMap.rank(-1);


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


  static Parameter<TreePhasespace,Energy> interfaceM0
    ("M0",
     "Set the cut below which flat virtuality sammpling is imposed.",
     &TreePhasespace::M0, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<TreePhasespace,Energy> interfaceMC
    ("MC",
     "Set the cut below which no virtualities are generated.",
     &TreePhasespace::Mc, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<TreePhasespace,bool> interfaceIncludeMirrored
    ("IncludeMirrored",
     "Choose whether to include mirrored diagrams for PS generation",
     &TreePhasespace::theIncludeMirrored, true, true, false);
  static SwitchOption interfaceIncludeMirroredYes
    (interfaceIncludeMirrored,
     "Yes",
     "Use unmirrored and mirrored diagrams",
     true);
  static SwitchOption interfaceIncludeMirroredNo
    (interfaceIncludeMirrored,
     "No",
     "Use only unmirrored diagrams",
     false);
  interfaceIncludeMirrored.rank(-1);

}

