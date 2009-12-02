// -*- C++ -*-
//
// ClusterHadronizationHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterHadronizationHandler class.
//

#include "ClusterHadronizationHandler.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Interface/Parameter.h> 
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Handlers/EventHandler.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Handlers/Hint.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/PDT/PDT.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Utilities/Throw.h>
#include <ThePEG/PDT/RemnantDecayer.h>
#include "Herwig++/Utilities/EnumParticles.h"
#include "CluHadConfig.h"
#include "Cluster.h"  
#include <iostream>
#include <cassert>

using namespace Herwig;

void ClusterHadronizationHandler::persistentOutput(PersistentOStream & os) 
  const {
  os << _partonSplitter 
     << _clusterFinder
     << _colourReconnector
     << _clusterFissioner
     << _lightClusterDecayer
     << _clusterDecayer
     << ounit(_minVirtuality2,GeV2)
     << ounit(_maxDisplacement,mm)
     << _underlyingEventHandler;
}


void ClusterHadronizationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitter 
     >> _clusterFinder
     >> _colourReconnector
     >> _clusterFissioner
     >> _lightClusterDecayer
     >> _clusterDecayer
     >> iunit(_minVirtuality2,GeV2)
     >> iunit(_maxDisplacement,mm)
     >> _underlyingEventHandler;
}

ClassDescription<ClusterHadronizationHandler> ClusterHadronizationHandler::initClusterHadronizationHandler;
// Definition of the static class description member.


void ClusterHadronizationHandler::Init() {

  static ClassDocumentation<ClusterHadronizationHandler> documentation
    ("This is the main handler class for the Cluster Hadronization",
     "The hadronization was performed using the cluster model of \\cite{Webber:1983if}.",
     "%\\cite{Webber:1983if}\n"
     "\\bibitem{Webber:1983if}\n"
     "  B.~R.~Webber,\n"
     "  ``A QCD Model For Jet Fragmentation Including Soft Gluon Interference,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 238}, 492 (1984).\n"
     "  %%CITATION = NUPHA,B238,492;%%\n"
     // main manual
     );

  static Reference<ClusterHadronizationHandler,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
		      "A reference to the PartonSplitter object", 
		      &Herwig::ClusterHadronizationHandler::_partonSplitter,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
		      "A reference to the ClusterFinder object", 
		      &Herwig::ClusterHadronizationHandler::_clusterFinder,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ColourReconnector> 
    interfaceColourReconnector("ColourReconnector", 
		      "A reference to the ColourReconnector object", 
		      &Herwig::ClusterHadronizationHandler::_colourReconnector,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
		      "A reference to the ClusterFissioner object", 
		      &Herwig::ClusterHadronizationHandler::_clusterFissioner,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
		    "A reference to the LightClusterDecayer object", 
		    &Herwig::ClusterHadronizationHandler::_lightClusterDecayer,
		    false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
		       "A reference to the ClusterDecayer object", 
		       &Herwig::ClusterHadronizationHandler::_clusterDecayer,
		       false, false, true, false);

  static Parameter<ClusterHadronizationHandler,Energy2> interfaceMinVirtuality2 
    ("MinVirtuality2",
     "Minimum virtuality^2 of partons to use in calculating distances  (unit [GeV2]).",
     &ClusterHadronizationHandler::_minVirtuality2, GeV2, 0.1*GeV2, ZERO, 10.0*GeV2,false,false,false);
  
  static Parameter<ClusterHadronizationHandler,Length> interfaceMaxDisplacement 
    ("MaxDisplacement",
     "Maximum displacement that is allowed for a particle  (unit [millimeter]).",
     &ClusterHadronizationHandler::_maxDisplacement, mm, 1.0e-10*mm, 
     0.0*mm, 1.0e-9*mm,false,false,false);

  static Reference<ClusterHadronizationHandler,StepHandler> interfaceUnderlyingEventHandler
    ("UnderlyingEventHandler",
     "Pointer to the handler for the Underlying Event. "
     "Set to NULL to disable.",
     &ClusterHadronizationHandler::_underlyingEventHandler, false, false, true, true, false);
}

void ClusterHadronizationHandler::doinitrun() {
  HadronizationHandler::doinitrun();
  // The run initialization is used here to all Cluster to have access to the
  // ClusterHadronizationHandler class instance, via a static pointer.
  Cluster::setPointerClusterHadHandler(this);
}

namespace {
  void extractChildren(tPPtr p, set<PPtr> & all) {
    if (p->children().empty()) return;

    for (PVector::const_iterator child = p->children().begin();
	 child != p->children().end(); ++child) {
      all.insert(*child);
      extractChildren(*child, all);
    }
  }
}

void ClusterHadronizationHandler::
handle(EventHandler & ch, const tPVector & tagged,
       const Hint &) {
  useMe();
  PVector currentlist(tagged.begin(),tagged.end());
  // set the scale for coloured particles to just above the gluon mass squared
  // if less than this so they are classed as perturbative
  Energy2 Q02 = 1.01*sqr(getParticleData(ParticleID::g)->constituentMass());
  for(unsigned int ix=0;ix<currentlist.size();++ix) {
    if(currentlist[ix]->scale()<Q02) currentlist[ix]->scale(Q02);
  }
  // split the gluons
  _partonSplitter->split(currentlist);
  // form the clusters
  ClusterVector clusters =
    _clusterFinder->formClusters(currentlist);

  _clusterFinder->reduceToTwoComponents(clusters); 

  // perform colour reconnection if needed and then
  // decay the clusters into one hadron
  bool lightOK = false;
  short tried = 0;
  const ClusterVector savedclusters = clusters;
  tPVector finalHadrons; // only needed for partonic decayer
  while (!lightOK && tried++ < 10) {

    _colourReconnector->rearrange(ch,clusters);

    finalHadrons = _clusterFissioner->fission(clusters,isSoftUnderlyingEventON());

    lightOK = _lightClusterDecayer->decay(clusters,finalHadrons);

    if (!lightOK) {
      clusters = savedclusters;
      for_each(clusters.begin(), 
	       clusters.end(), 
	       mem_fun(&Particle::undecay));
    }
  } 
  if (!lightOK)
    throw Exception("CluHad::handle(): tried LightClusterDecayer 10 times!", 
		    Exception::eventerror);



  // decay the remaining clusters
  _clusterDecayer->decay(clusters,finalHadrons);

  // *****************************************
  // *****************************************
  // *****************************************

  StepPtr pstep = newStep();
  set<PPtr> allDecendants;
  for (tPVector::const_iterator it = tagged.begin();
       it != tagged.end(); ++it) {
    extractChildren(*it, allDecendants);
  }

  for(set<PPtr>::const_iterator it = allDecendants.begin();
      it != allDecendants.end(); ++it) {
    // this is a workaround because the set sometimes 
    // re-orders parents after their children
    if ((*it)->children().empty())
      pstep->addDecayProduct(*it);
    else {
      pstep->addDecayProduct(*it);
      pstep->addIntermediate(*it);
    }
  }

  // *****************************************
  // *****************************************
  // *****************************************

  // soft underlying event if needed
  if (isSoftUnderlyingEventON()) {
    assert(_underlyingEventHandler);
    ch.performStep(_underlyingEventHandler,Hint::Default());
  }
  // zero all positions
  // extract all particles from the event
  tEventPtr event=ch.currentEvent();
  vector<tPPtr> particles;
  particles.reserve(256);
  event->select(back_inserter(particles), ThePEG::AllSelector());
  // and the final-state particles
  set<tPPtr> finalstate;
  event->selectFinalState(inserter(finalstate));
  for(vector<tPPtr>::const_iterator pit=particles.begin();
      pit!=particles.end();++pit) {
    // if a final-state particle just zero production
    if(finalstate.find(*pit)!=finalstate.end()) {
      (**pit).setVertex(LorentzPoint());
    }
    // if not zero the lot
    else {
      (**pit).setVertex(LorentzPoint());
      (**pit).setLifeLength(LorentzDistance());
    }
  }
}
