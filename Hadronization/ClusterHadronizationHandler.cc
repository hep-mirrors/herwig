// -*- C++ -*-
//
// ClusterHadronizationHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
#include <ThePEG/Handlers/Hint.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/PDT/PDT.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Utilities/Throw.h>
#include "Herwig++/Utilities/EnumParticles.h"
#include "CluHadConfig.h"
#include "Cluster.h"  
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

ClusterHadronizationHandler * ClusterHadronizationHandler::currentHandler_ = 0;

DescribeClass<ClusterHadronizationHandler,HadronizationHandler>
describeClusterHadronizationHandler("Herwig::ClusterHadronizationHandler","");

IBPtr ClusterHadronizationHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterHadronizationHandler::fullclone() const {
  return new_ptr(*this);
}

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
  currentHandler_ = this;
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

    // no colour reconnection with baryon-number-violating (BV) clusters
    ClusterVector CRclusters, BVclusters;
    CRclusters.reserve( clusters.size() );
    BVclusters.reserve( clusters.size() );
    for (size_t ic = 0; ic < clusters.size(); ++ic) {
      ClusterPtr cl = clusters.at(ic);
      bool hasClusterParent = false;
      for (unsigned int ix=0; ix < cl->parents().size(); ++ix) {
        if (cl->parents()[ix]->id() == ParticleID::Cluster) {
          hasClusterParent = true;
          break;
        }
      }
      if (cl->numComponents() > 2 || hasClusterParent) BVclusters.push_back(cl);
      else CRclusters.push_back(cl);
    }

    // colour reconnection
    _colourReconnector->rearrange(CRclusters);

    // tag new clusters as children of the partons to hadronize
    _setChildren(CRclusters);

    // recombine vectors of (possibly) reconnected and BV clusters
    clusters.clear();
    clusters.insert( clusters.end(), CRclusters.begin(), CRclusters.end() );
    clusters.insert( clusters.end(), BVclusters.begin(), BVclusters.end() );
    
    // fission of heavy clusters
    // NB: during cluster fission, light hadrons might be produced straight away
    finalHadrons = _clusterFissioner->fission(clusters,isSoftUnderlyingEventON());

    lightOK = _lightClusterDecayer->decay(clusters,finalHadrons);

    // if the decay of the light clusters was not successful, undo the cluster
    // fission and decay steps and revert to the original state of the event
    // record
    if (!lightOK) {
      clusters = savedclusters;
      for_each(clusters.begin(), 
	       clusters.end(), 
	       mem_fun(&Particle::undecay));
    }
  } 
  if (!lightOK) {
    //    currentHandler_ = 0;
    throw Exception("CluHad::handle(): tried LightClusterDecayer 10 times!", 
		    Exception::eventerror);
  }


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
  // calculate positions
  // extract all particles from the event
  tEventPtr event=ch.currentEvent();
  vector<tPPtr> particles;
  particles.reserve(256);
  event->select(back_inserter(particles), ThePEG::AllSelector());
  for(vector<tPPtr>::const_iterator pit=particles.begin();
      pit!=particles.end(); ++pit) {
    if ((**pit).parents().empty()) 
      continue;
    tPPtr parent = (**pit).parents()[0];
    // fudged fix for the shower's technical vertices:
    // make lifelength for 1-1 vertices zero
    // \todo sort out shower inserted vertices properly
    if ( (**pit).id() == parent->id() 
	 && (**pit).parents().size() == 1
	 && parent->children().size() == 1 ) {
      parent->setLifeLength(LorentzDistance());
      // ??? (**pit).setVertex( parent->vertex() );
    }
    // fix the vertex position for particles from clusters
    bool inClusters = false;
    bool newVertex  = false;
    // iterate up the ancestry to find the origin of the parent cluster
    while( !parent->parents().empty() ) {
      bool parentIsCluster = parent->id() == ParticleID::Cluster;
      if ( !inClusters && parentIsCluster ) {
	inClusters = true;
      }
      else if ( inClusters && !parentIsCluster ) {
	newVertex = true;
	break;
      }
      parent = parent->parents()[0]; 
    }
    // parent is now the ancestor of the cluster chain
    if ( newVertex ) (**pit).setVertex( parent->vertex() );
  }
}

void ClusterHadronizationHandler::_setChildren(ClusterVector clusters) const {

  // erase existing information about the partons' children
  tPVector partons;
  for (ClusterVector::const_iterator cl = clusters.begin();
       cl != clusters.end(); cl++) {
    partons.push_back( (*cl)->colParticle() );
    partons.push_back( (*cl)->antiColParticle() );
  }
  for_each(partons.begin(), partons.end(), mem_fun(&Particle::undecay));

  // give new parents to the clusters: their constituents
  for (ClusterVector::iterator cl = clusters.begin();
       cl != clusters.end(); cl++) {
    (*cl)->colParticle()->addChild(*cl);
    (*cl)->antiColParticle()->addChild(*cl);
  }

}

