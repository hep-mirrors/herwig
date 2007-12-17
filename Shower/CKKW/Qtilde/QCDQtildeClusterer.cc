// -*- C++ -*-
//
// QCDQtildeClusterer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QCDQtildeClusterer class.
//

#include "QCDQtildeClusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "Herwig++/Shower/CKKW/Reweighting/DefaultJetMeasure.h"
#include "Herwig++/Shower/CKKW/Clustering/ToIncomingCMS.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QCDQtildeClusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

QCDQtildeClusterer::~QCDQtildeClusterer() {}

void QCDQtildeClusterer::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void QCDQtildeClusterer::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<QCDQtildeClusterer> QCDQtildeClusterer::initQCDQtildeClusterer;
// Definition of the static class description member.

void QCDQtildeClusterer::Init() {

  static ClassDocumentation<QCDQtildeClusterer> documentation
    ("QCD clusterings in the qtilde formalism");

}


ClusteringPtr QCDQtildeClusterer::doScale (const vector<tClusteringParticlePtr>& children,
					   const vector<ClusteringParticlePtr>& parents,
					   const tClusteringConfigurationPtr& configuration) {

  EmitterSpectatorConfigurationPtr config =
    dynamic_ptr_cast<EmitterSpectatorConfigurationPtr>(configuration);

  if (!config)
    throw Exception() << "CKKW : QCDQtildeClusterer::doScale : expecting EmitterSpectatorConfiguration"
		      << Exception::runerror;

  JetMeasurePtr jetm = jetMeasure();

  DefaultJetMeasurePtr jetmeasure = 
    dynamic_ptr_cast<DefaultJetMeasurePtr>(jetm);

  if (!jetmeasure)
    throw Exception() << "CKKW : QCDQtildeClusterer::doScale : expecting DefaultJetMeasure"
		      << Exception::runerror;

  EmitterSpectatorClusteringPtr clustering =
    new_ptr(EmitterSpectatorClustering(children,
				       config->spectatorBeforeClusteringIndex(),
				       parents,config->spectatorAfterClusteringIndex(),
				       this,configuration));
  clustering->eventGenerator(generator());

  clustering->generateSudakovBasis();
  jetmeasure->sudakovBasis(clustering->sudakovBasis().first,clustering->sudakovBasis().second);

  Energy2 scale = jetmeasure->scale(children,parents,config);
  double z = jetmeasure->z(children,parents,config);

  clustering->scale(scale);
  clustering->alphaScale(scale);
  clustering->momentumFraction(z);

  if(clustering->emitter()->pData().partonId.state == ClusteringParticleState::initial) {

    clustering->emitter()->x(z * clustering->emitter()->x());
    clustering->postClustering(new_ptr(ToIncomingCMS()));

  }

  clustering->weight(1.);

  return clustering;

}

void QCDQtildeClusterer::doKinematics (const tClusteringPtr& clustering) {

  EmitterSpectatorClusteringPtr cl =
    dynamic_ptr_cast<EmitterSpectatorClusteringPtr>(clustering);

  if (!cl)
    throw Exception() << "CKKW : QCDQtildeClusterer::doKinematics : expecting EmitterSpectatorClustering"
		      << Exception::runerror;

  Energy2 scale = cl->scale();

  cl->emission().first->productionScale(scale);
  cl->emission().second->productionScale(scale);
  cl->emitter()->splittingScale(scale);

  Lorentz5Momentum emerging;
  if (cl->emission().first->pData().partonId.state == ClusteringParticleState::final &&
      cl->emission().second->pData().partonId.state == ClusteringParticleState::final)
    emerging = cl->emission().first->momentum() + cl->emission().second->momentum();
  if (cl->emission().first->pData().partonId.state == ClusteringParticleState::initial &&
      cl->emission().second->pData().partonId.state == ClusteringParticleState::final) {
    emerging = cl->emission().first->momentum() - cl->emission().second->momentum();
  }
  if (cl->emission().first->pData().partonId.state == ClusteringParticleState::final &&
      cl->emission().second->pData().partonId.state == ClusteringParticleState::initial) {
    emerging = - cl->emission().first->momentum() + cl->emission().second->momentum();
  }
  cl->emitter()->momentum(emerging);

  // nothing happened to the spectator
  cl->spectatorBeforeClustering()->productionScale(cl->spectatorBeforeClustering()->splittingScale());
  cl->spectatorBeforeClustering()->setNoReweight();
  cl->spectatorAfterClustering()->splittingScale(cl->spectatorBeforeClustering()->splittingScale());
  cl->spectatorAfterClustering()->momentum(cl->spectatorBeforeClustering()->momentum());

}


