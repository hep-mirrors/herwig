// -*- C++ -*-
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
#include "Remnant.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "CluHadConfig.h"
#include "Cluster.h"  
#include "Remnant.h"
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
     << _forcedSplitter
     << ounit(_minVirtuality2,GeV2)
     << ounit(_maxDisplacement,millimeter)
     << _softUnderlyingEventMode
     << _underlyingEventHandler;
}


void ClusterHadronizationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitter 
     >> _clusterFinder
     >> _colourReconnector
     >> _clusterFissioner
     >> _lightClusterDecayer
     >> _clusterDecayer
     >> _forcedSplitter
     >> iunit(_minVirtuality2,GeV2)
     >> iunit(_maxDisplacement,millimeter)
     >> _softUnderlyingEventMode
     >> _underlyingEventHandler;
}

ClassDescription<ClusterHadronizationHandler> ClusterHadronizationHandler::initClusterHadronizationHandler;
// Definition of the static class description member.


void ClusterHadronizationHandler::Init() {

  static ClassDocumentation<ClusterHadronizationHandler> documentation
    ("This is the main handler class for the Cluster Hadronization");

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

  static Reference<ClusterHadronizationHandler,ForcedSplitting> interfaceForcedSplitting
    ("ForcedSplitting",
     "Object responsible for the forced splitting of the Remnant",
     &ClusterHadronizationHandler::_forcedSplitter, false, false, true, false, false);

  static Parameter<ClusterHadronizationHandler,Energy2> interfaceMinVirtuality2 
    ("MinVirtuality2",
     "Minimum virtuality^2 of partons to use in calculating distances  (unit [GeV2]).",
     &ClusterHadronizationHandler::_minVirtuality2, GeV2, 0.1*GeV2, 0.0*GeV2, 10.0*GeV2,false,false,false);
  
  static Parameter<ClusterHadronizationHandler,Length> interfaceMaxDisplacement 
    ("MaxDisplacement",
     "Maximum displacement that is allowed for a particle  (unit [millimeter]).",
     &ClusterHadronizationHandler::_maxDisplacement, millimeter, 1.0e-10*millimeter, 
     0.0*millimeter, 1.0e-9*millimeter,false,false,false);

  static Switch<ClusterHadronizationHandler, bool> interfaceSoftUnderlyingEventMode
    ("OnOffSoftUnderlyingEventMode",
     "Choice of the soft underlying event switch mode.",
     &ClusterHadronizationHandler::_softUnderlyingEventMode, 0, false, false);
  static SwitchOption interfaceSoftUnderlyingEventMode0
    (interfaceSoftUnderlyingEventMode,"SoftUnderlyingEvent-OFF", 
     "soft underlying event is OFF", false);
  static SwitchOption interfaceSoftUnderlyingEventMode1
    (interfaceSoftUnderlyingEventMode,"SoftUnderlyingEvent-ON",
     "soft underlying event is ON", true);

  static Reference<ClusterHadronizationHandler,StepHandler> interfaceUnderlyingEventHandler
    ("UnderlyingEventHandler",
     "Pointer to the handler for the Underlying Event. "
     "This must be set when OnOffSoftUnderlyingEventMode is 1.",
     &ClusterHadronizationHandler::_underlyingEventHandler, false, false, true, true, false);
}


void ClusterHadronizationHandler::doinit() throw(InitException) {
  HadronizationHandler::doinit();
  // check that UEHandler is there when _softUnderlyingEventMode is true
  if (!_underlyingEventHandler && isSoftUnderlyingEventON())
    Throw<InitException>() << "OnOffSoftUnderlyingEventMode is set to 1, "
			   << "this requires the pointer\n" 
			   << name() << ":UnderlyingEventHandler to be set.";
}

void ClusterHadronizationHandler::doinitrun() {
  HadronizationHandler::doinitrun();
  // The run initialization is used here to all Cluster to have access to the
  // ClusterHadronizationHandler class instance, via a static pointer.
  Cluster::setPointerClusterHadHandler(this);
}

void ClusterHadronizationHandler::
handle(EventHandler & ch, const tPVector & tagged,
       const Hint &) throw(Veto, Stop, Exception) {
  ClusterVector clusters;
  StepPtr pstep = newStep();
  // split the remnants if needed
  tPVector partonsA=_forcedSplitter->split(tagged,pstep);
  // split the gluons
  tPVector partons=_partonSplitter->split(partonsA,pstep);
  // force a new step to form the clusters
  pstep = ch.newStep(this);
  _clusterFinder->formClusters(ch.currentCollision(),pstep,partons,clusters); 
  _clusterFinder->reduceToTwoComponents(pstep,clusters); 
  // perform colour reconnection if needed and then
  // decay the clusters into one hadron
  bool lightOK = false;
  short tried = 0;
  while (!lightOK && tried++ < 10) {
    // force a new step that can be erased
    pstep = ch.newStep(this);
    _colourReconnector->rearrange(ch,pstep,clusters);
    _clusterFissioner->fission(pstep, isSoftUnderlyingEventON());
    // force a new step that can be erased
    pstep = ch.newStep(this);
    lightOK = _lightClusterDecayer->decay(pstep);
    if (!lightOK) {
      ch.popStep(); 
      ch.popStep();
    }
  } 
  if (!lightOK)
    throw Exception("CluHad::handle(): tried LightClusterDecayer 10 times!", Exception::eventerror);
  
  // decay the remaining clusters
  _clusterDecayer->decay(pstep);

  // soft underlying event if needed
  if (isSoftUnderlyingEventON()) {
    assert(_underlyingEventHandler);
    ch.performStep(_underlyingEventHandler,Hint::Default());
  }
}
