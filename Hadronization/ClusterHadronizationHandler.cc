// -*- C++ -*-
//
// ClusterHadronizationHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include <ThePEG/Interface/Switch.h>
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
#include "Herwig/Utilities/EnumParticles.h"
#include "CluHadConfig.h"
#include "Cluster.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include <ThePEG/Interface/Command.h>

using namespace Herwig;

ClusterHadronizationHandler * ClusterHadronizationHandler::currentHandler_ = 0;

DescribeClass<ClusterHadronizationHandler,HadronizationHandler>
describeClusterHadronizationHandler("Herwig::ClusterHadronizationHandler","Herwig.so");

IBPtr ClusterHadronizationHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterHadronizationHandler::fullclone() const {
  return new_ptr(*this);
}

void ClusterHadronizationHandler::persistentOutput(PersistentOStream & os)
  const {
  os << _partonSplitters << _clusterFinders << _colourReconnectors
     << _clusterFissioners << _lightClusterDecayers << _clusterDecayers
     << _partonSplitter << _clusterFinder << _colourReconnector
     << _clusterFissioner << _lightClusterDecayer << _clusterDecayer
     << reshuffle_ << reshuffleMode_ << gluonMassGenerator_
     << ounit(_minVirtuality2,GeV2) << ounit(_maxDisplacement,mm)
     << _underlyingEventHandler << _reduceToTwoComponents;
}


void ClusterHadronizationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitters >> _clusterFinders >> _colourReconnectors
     >> _clusterFissioners >> _lightClusterDecayers >> _clusterDecayers
     >> _partonSplitter >> _clusterFinder >> _colourReconnector
     >> _clusterFissioner >> _lightClusterDecayer >> _clusterDecayer
     >> reshuffle_ >> reshuffleMode_ >> gluonMassGenerator_
     >> iunit(_minVirtuality2,GeV2) >> iunit(_maxDisplacement,mm)
     >> _underlyingEventHandler >> _reduceToTwoComponents;
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

  static Reference<ClusterHadronizationHandler,GluonMassGenerator> interfaceGluonMassGenerator
    ("GluonMassGenerator",
     "Set a reference to a gluon mass generator.",
     &ClusterHadronizationHandler::gluonMassGenerator_, false, false, true, true, false);

  static Switch<ClusterHadronizationHandler,bool> interfaceReshuffle
    ("Reshuffle",
     "Perform reshuffling if constituent masses have not yet been included by the shower",
     &ClusterHadronizationHandler::reshuffle_, false, false, false);
  static SwitchOption interfaceReshuffleYes
    (interfaceReshuffle,
     "Global",
     "Do reshuffle.",
     true);
  static SwitchOption interfaceReshuffleNo
    (interfaceReshuffle,
     "No",
     "Do not reshuffle.",
     false);
   
  static Switch<ClusterHadronizationHandler,int> interfaceReshuffleMode
    ("ReshuffleMode",
     "Which mode is used for the reshuffling to constituent masses",
     &ClusterHadronizationHandler::reshuffleMode_, 0, false, false);
  static SwitchOption interfaceReshuffleModeGlobal
    (interfaceReshuffleMode,
     "Global",
     "Global reshuffling on all final state partons",
     0);
  static SwitchOption interfaceReshuffleModeColourConnected
    (interfaceReshuffleMode,
     "ColourConnected",
     "Separate reshuffling for colour connected partons",
     1);

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

  // functions to move the handlers so they are set
  
  
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

   static Switch<ClusterHadronizationHandler,bool> interfaceReduceToTwoComponents
    ("ReduceToTwoComponents",
     "Whether or not to reduce three component baryon-number violating clusters to two components before cluster splitting or leave"
     " this till after the cluster splitting",
     &ClusterHadronizationHandler::_reduceToTwoComponents, true, false, false);
  static SwitchOption interfaceReduceToTwoComponentsYes
    (interfaceReduceToTwoComponents,
     "BeforeSplitting",
     "Reduce to two components",
     true);
  static SwitchOption interfaceReduceToTwoComponentsNo
    (interfaceReduceToTwoComponents,
     "AfterSplitting",
     "Treat as three components",
     false);

  static Command<ClusterHadronizationHandler> interfaceUseHandlersForInteraction
    ("UseHandlersForInteraction",
     "Use the current set of handlers for the given interaction.",
     &ClusterHadronizationHandler::_setHandlersForInteraction, false);

}

void ClusterHadronizationHandler::rebind(const TranslationMap & trans) {
  for (auto iter : _partonSplitters) {
     _partonSplitters[iter.first] = trans.translate(iter.second);
  }
  for (auto iter : _clusterFinders) {
     _clusterFinders[iter.first] = trans.translate(iter.second);
  }
  for (auto iter : _colourReconnectors) {
     _colourReconnectors[iter.first] = trans.translate(iter.second);
  }
  for (auto iter : _clusterFissioners) {
     _clusterFissioners[iter.first] = trans.translate(iter.second);
  }
  for (auto iter : _lightClusterDecayers) {
     _lightClusterDecayers[iter.first] = trans.translate(iter.second);
  }
  for (auto iter : _clusterDecayers) {
     _clusterDecayers[iter.first] = trans.translate(iter.second);
  }
  HandlerBase::rebind(trans);
}

IVector ClusterHadronizationHandler::getReferences() {
  IVector ret = HandlerBase::getReferences();
  for (auto iter : _partonSplitters) {
     ret.push_back(iter.second);
  }
  for (auto iter : _clusterFinders) {
     ret.push_back(iter.second);
  }
  for (auto iter : _colourReconnectors) {
     ret.push_back(iter.second);
  }
  for (auto iter : _clusterFissioners) {
     ret.push_back(iter.second);
  }
  for (auto iter : _lightClusterDecayers) {
     ret.push_back(iter.second);
  }
  for (auto iter : _clusterDecayers) {
     ret.push_back(iter.second);
  }
  return ret;
}

void ClusterHadronizationHandler::doinit() {
  HadronizationHandler::doinit();
  // put current handlers into action for QCD to guarantee default behaviour
  if ( _partonSplitters.empty() ) {
    _partonSplitters[PDT::ColouredQCD] = _partonSplitter;
    _clusterFinders[PDT::ColouredQCD] = _clusterFinder;
    _colourReconnectors[PDT::ColouredQCD] = _colourReconnector;
    _clusterFissioners[PDT::ColouredQCD] = _clusterFissioner;
    _lightClusterDecayers[PDT::ColouredQCD] = _lightClusterDecayer;
    _clusterDecayers[PDT::ColouredQCD] = _clusterDecayer; 
  }
}

string ClusterHadronizationHandler::_setHandlersForInteraction(string interaction) {
  int interactionId = -1;
  if ( interaction == " QCD" ) {
    interactionId = PDT::ColouredQCD;
  } else {
    return "unknown interaction " + interaction;
  }
  _partonSplitters[interactionId] = _partonSplitter;
  _clusterFinders[interactionId] = _clusterFinder;
  _colourReconnectors[interactionId] = _colourReconnector;
  _clusterFissioners[interactionId] = _clusterFissioner;
  _lightClusterDecayers[interactionId] = _lightClusterDecayer;
  _clusterDecayers[interactionId] = _clusterDecayer; 
  return "";
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

  PVector theList(tagged.begin(),tagged.end());
  
  // set the scale for coloured particles to just above the gluon mass squared
  // if less than this so they are classed as perturbative

  // TODO: Should this be hard-coded here, or defined in inputs?
  map<int,int> interactions = {{PDT::ColouredQCD, ParticleID::g}};
  
  // PVector currentlist(tagged.begin(),tagged.end());
  map<int,PVector> currentlists;
  for ( const auto & p : theList ) {
    for ( auto i : interactions ) {
      if ( p->data().colouredInteraction() == i.first ) {
	currentlists[i.first].push_back(p);
    Energy2 Q02 = 1.01*sqr(getParticleData(i.second)->constituentMass());
	if(currentlists[i.first].back()->scale()<Q02) currentlists[i.first].back()->scale(Q02);
      }
    }
  }

  map<int,ClusterVector> clusters;

  // ATTENTION this needs to have changes to include the new hadron spectrum functionality (gluon mass generator)
  for ( auto i : interactions ) {

    // reshuffle to the constituents
    if ( reshuffle_ ) {
    
      vector<PVector> reshufflelists;

      if (reshuffleMode_==0){ // global reshuffling
	reshufflelists.push_back(currentlists[i.first]);
      }
      else if (reshuffleMode_==1){// colour connected reshuffling
	splitIntoColourSinglets(currentlists[i.first], reshufflelists, i.first);
      }

      for (auto currentlist : reshufflelists){
	// get available energy and energy needed for constituent mass shells
	LorentzMomentum totalQ;
	Energy needQ = ZERO;
	size_t nGluons = 0; // number of gluons for which a mass need be generated
	for ( auto p : currentlist ) {
	  totalQ += p->momentum();
	  if ( p->id() == i.second && gluonMassGenerator() ) {
	    ++nGluons;
	    continue;
	  }
	  needQ += p->dataPtr()->constituentMass();
	}
	Energy Q = totalQ.m();
	if ( needQ > Q )
	  throw Exception() << "cannot reshuffle to constituent mass shells" << Exception::eventerror;

	// generate gluon masses if needed
	list<Energy> gluonMasses;
	if ( nGluons && gluonMassGenerator() )
	  gluonMasses = gluonMassGenerator()->generateMany(nGluons,Q-needQ);

	// set masses for inidividual particles
	vector<Energy> masses;
	for ( auto p : currentlist ) {
	  if ( p->id() == i.second && gluonMassGenerator() ) {
	    list<Energy>::const_iterator it = gluonMasses.begin();
	    advance(it,UseRandom::irnd(gluonMasses.size()));
	    masses.push_back(*it);
	    gluonMasses.erase(it);
	  } 
	  else {
	    masses.push_back(p->dataPtr()->constituentMass());
	  }
	}

	// reshuffle to new masses
	reshuffle(currentlist,masses);

      }

    }

    // split the gluons
    if ( !currentlists[i.first].empty() ) {
      assert(_partonSplitters.find(i.first) != _partonSplitters.end());
      _partonSplitters[i.first]->split(currentlists[i.first]);
    }

    // form the clusters
    if ( !currentlists[i.first].empty() ) {
      assert(_clusterFinders.find(i.first) != _clusterFinders.end());
      clusters[i.first] = _clusterFinders[i.first]->formClusters(currentlists[i.first]);
      // reduce BV clusters to two components now if needed
      if(_reduceToTwoComponents)
	_clusterFinders[i.first]->reduceToTwoComponents(clusters[i.first]);
    }
  }

  // perform colour reconnection if needed and then
  // decay the clusters into one hadron
  for ( auto i : interactions ) {
    if ( clusters[i.first].empty() )
      continue;
    bool lightOK = false;
    short tried = 0;
    const ClusterVector savedclusters = clusters[i.first];
    tPVector finalHadrons; // only needed for partonic decayer
    while (!lightOK && tried++ < 10) {
      // no colour reconnection with baryon-number-violating (BV) clusters
      ClusterVector CRclusters, BVclusters;
      CRclusters.reserve( clusters[i.first].size() );
      BVclusters.reserve( clusters[i.first].size() );
      for (size_t ic = 0; ic < clusters[i.first].size(); ++ic) {
	ClusterPtr cl = clusters[i.first].at(ic);
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
      assert(_colourReconnectors.find(i.first) != _colourReconnectors.end());
      _colourReconnectors[i.first]->rearrange(CRclusters);


      // tag new clusters as children of the partons to hadronize
      _setChildren(CRclusters);
    
   
      // forms diquarks
      // we should already have used  _clusterFinder[i.first]
      _clusterFinders[i.first]->reduceToTwoComponents(CRclusters);
    
      // recombine vectors of (possibly) reconnected and BV clusters
      clusters[i.first].clear();
      clusters[i.first].insert( clusters[i.first].end(), CRclusters.begin(), CRclusters.end() );
      clusters[i.first].insert( clusters[i.first].end(), BVclusters.begin(), BVclusters.end() );

      // fission of heavy clusters
      // NB: during cluster fission, light hadrons might be produced straight away
      assert(_clusterFissioners.find(i.first) != _clusterFissioners.end());
      finalHadrons = _clusterFissioners[i.first]->fission(clusters[i.first],i.first == PDT::ColouredQCD ? isSoftUnderlyingEventON() : false);


      // if clusters not previously reduced to two components do it now
      if(!_reduceToTwoComponents)
	_clusterFinders[i.first]->reduceToTwoComponents(clusters[i.first]);
      assert(_lightClusterDecayers.find(i.first) != _lightClusterDecayers.end());
      lightOK = _lightClusterDecayers[i.first]->decay(clusters[i.first],finalHadrons);
      // if the decay of the light clusters was not successful, undo the cluster
      // fission and decay steps and revert to the original state of the event
      // record
      if (!lightOK) {
	clusters[i.first] = savedclusters;
	for_each(clusters[i.first].begin(),
		 clusters[i.first].end(),
		 std::mem_fn(&Particle::undecay));
      }
    }
    if (!lightOK) {
      throw Exception( "CluHad::handle(): tried LightClusterDecayer 10 times!",
		       Exception::eventerror);
    }

    // decay the remaining clusters
    assert(_clusterDecayers[i.first]);
    _clusterDecayers[i.first]->decay(clusters[i.first],finalHadrons);
    
  }

  // *****************************************
  // *****************************************
  // *****************************************

  bool finalStateCluster=false;
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
    if ((*it)->children().empty()){
      // If there is a cluster in the final state throw an event error
      if((*it)->id()==81) {
         finalStateCluster=true;
      }
      pstep->addDecayProduct(*it);
    }
    else {
      pstep->addDecayProduct(*it);
      pstep->addIntermediate(*it);
    }
  }

  // For very small center of mass energies it might happen that baryonic clusters cannot decay into hadrons
  if (finalStateCluster){
     throw Exception( "CluHad::Handle(): Cluster in the final state", 
                     Exception::eventerror);
  }
  // *****************************************
  // *****************************************
  // *****************************************

  // soft underlying event if needed
  if (isSoftUnderlyingEventON()) {
    assert(_underlyingEventHandler);
    ch.performStep(_underlyingEventHandler,Hint::Default());
  }
}


// Sets parent child relationship of all clusters with two components
// Relationships for clusters with more than two components are set elsewhere in the Colour Reconnector
void ClusterHadronizationHandler::_setChildren(const ClusterVector & clusters) const {
  // erase existing information about the partons' children
  tPVector partons;
  for ( const auto & cl : clusters ) {
    if ( cl->numComponents() > 2 ) continue;
    partons.push_back( cl->colParticle() );
    partons.push_back( cl->antiColParticle() );
  }
  // erase all previous information about parent child relationship
  for_each(partons.begin(), partons.end(), std::mem_fn(&Particle::undecay));

  // give new parents to the clusters: their constituents
  for ( const auto & cl : clusters ) {
    if ( cl->numComponents() > 2 ) continue;
    cl->colParticle()->addChild(cl);
    cl->antiColParticle()->addChild(cl);
  }
}

void ClusterHadronizationHandler::splitIntoColourSinglets(PVector copylist,
							  vector<PVector>& reshufflelists,
							  int){
 
  PVector currentlist;
  bool gluonloop;
  PPtr firstparticle, temp;
  reshufflelists.clear();
 
  while (copylist.size()>0){
    gluonloop=false;
    currentlist.clear();
 
    firstparticle=copylist.back();  
    copylist.pop_back();     
 
    if (!firstparticle->coloured()){
      continue; //non-coloured particles are not included
    }
     
    currentlist.push_back(firstparticle);
 
    //go up the anitColourLine and check if we are in a gluon loop
    temp=firstparticle;
    while( temp->hasAntiColour()){
      temp = temp->antiColourLine()->endParticle();
      if(temp==firstparticle){
	gluonloop=true;
	break;
      }
      else{
	currentlist.push_back(temp);
	copylist.erase(remove(copylist.begin(),copylist.end(), temp), copylist.end());
      }
    }
 
    //if not a gluon loop, go up the ColourLine
    if(!gluonloop){
      temp=firstparticle;
      while( temp->hasColour()){
	temp=temp->colourLine()->startParticle();
	currentlist.push_back(temp);
	copylist.erase(remove(copylist.begin(),copylist.end(), temp), copylist.end());
      }
    }
 
    reshufflelists.push_back(currentlist); 
  }
 
}
