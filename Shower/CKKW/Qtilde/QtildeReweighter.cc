// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildeReweighter class.
//

#include "QtildeReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtildeReweighter.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/Shower/CKKW/Clustering/Partitioner.h"
#include "Herwig++/Shower/CKKW/Clustering/EmitterSpectatorClustering.h"

using namespace Herwig;

QtildeReweighter::~QtildeReweighter() {}

void QtildeReweighter::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void QtildeReweighter::persistentInput(PersistentIStream & , int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<QtildeReweighter> QtildeReweighter::initQtildeReweighter;
// Definition of the static class description member.

void QtildeReweighter::Init() {

  static ClassDocumentation<QtildeReweighter> documentation
    ("A reweighter for ME/PS merging with the qtilde shower.");

}

void QtildeReweighter::analyzeHistory (CascadeHistory history) {

#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== QtildeReweighter::analyzeHistory" << endl
		     << "listing history given:" << endl
		     << "reconstrcuted particles" << endl;
  for(vector<ClusteringParticlePtr>::iterator c = history.reconstructed.begin(); c != history.reconstructed.end();
      ++c)
    (**c).debugDump(generator()->log());
  generator()->log() << "clusterings performed" << endl;
  for(list<ClusteringPtr>::iterator c = history.clusterings.begin(); c != history.clusterings.end();
      ++c)
    (**c).debugDump(generator()->log());

#endif

  // check for hadronic collision

  bool hadronic = false;

  // consider only coloured external lines

  vector<ClusteringParticlePtr> external;

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "looking for external particles to consider ...";
#endif
  for (vector<ClusteringParticlePtr>::iterator p = history.reconstructed.begin();
       p != history.reconstructed.end(); ++p) {
    if ((**p).children().empty()) {
	if ((**p).pData().colour != 0 || (**p).pData().antiColour != 0) {
	  external.push_back(*p);
	  if ((**p).pData().partonId.state == ClusteringParticleState::initial) hadronic = true;
	}
    }
    else break; // no more external ones from now on, see CascadeReconstructor.cc
  }

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << " found " << external.size() << endl;
#endif

  // for each, trace back until creation vertex  

  for (vector<ClusteringParticlePtr>::iterator p = external.begin(); p != external.end(); ++p) {

    // things we can't do yet

    // tops
    if (abs((**p).pData().partonId.PDGId) == 7)
      throw Exception() << "CKKW : QtildeReweighter::analyzeHistory : "
			<< "Top quarks cannot be handled yet." << Exception::runerror;

    // coloured stuff other than quarks and gluons
    if ((**p).pData().colour != (**p).pData().antiColour &&
	abs((**p).pData().partonId.PDGId) > 6 && (**p).pData().partonId.PDGId != 21)
      throw Exception() << "CKKW : QtildeReweighter::analyzeHistory : "
			<< "Coloured particles other than quarks or gluons cannot be handled yet."
			<< Exception::runerror;

    // things we can do ...

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "starting from" << endl;
  (**p).debugDump(generator()->log());
#endif
 
    // if not clustered at all, set shower scale from production scale
    if (!(**p).production()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "not clustered, shower scle / GeV2 -> " << (**p).productionScale()/GeV2  << endl;
#endif
      (**p).showerScale((**p).productionScale());
      continue;
    }

    ClusteringParticlePtr current = *p;

    while(true) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "looking at " << endl;
      current->debugDump(generator()->log());
#endif

      // hard process
      if (!current->production()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "hard process: stop."  << endl;
#endif
	break;
      }

      // if spectator, go on with the next clustering
      if (current == 
	  dynamic_ptr_cast<tcEmitterSpectatorClusteringPtr>
	  (current->production())->spectatorBeforeClustering()) {
	current = dynamic_ptr_cast<tcEmitterSpectatorClusteringPtr>
	  (current->production())->spectatorAfterClustering();
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "spectator, continuing"  << endl;
#endif
	continue;
      }

      // get emitter and emissions
      tClusteringParticlePtr emitter = dynamic_ptr_cast<tcEmitterSpectatorClusteringPtr>
	(current->production())->emitter();
      pair<tClusteringParticlePtr,tClusteringParticlePtr> emission =
	dynamic_ptr_cast<tcEmitterSpectatorClusteringPtr>
	(current->production())->emission();
      tClusteringParticlePtr partner;
      if (current == emission.first) partner = emission.second;
      else partner = emission.first;

      // singlet emission and coloured emitter
      if (partner->pData().colour == partner->pData().antiColour &&
	  emitter->pData().colour != emitter->pData().antiColour) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "singlet emitted, continuing"  << endl;
#endif
	current = emitter;
	continue;
      }

      // singlet emitter
      if (emitter->pData().colour == emitter->pData().antiColour) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "originated from singlet. stop."  << endl;
#endif
	// reset the hard scale according to the KtTilde measure
	vector<tClusteringParticlePtr> decay;
	decay.push_back(current);
	decay.push_back(partner);
	resolution()->hardScales(decay);
	break;
      }
      
      // final emission and initial emitter
      if (current->pData().partonId.state == ClusteringParticleState::final &&
	  emitter->pData().partonId.state == ClusteringParticleState::initial) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "final line originated from initial line. stop."  << endl;
#endif
	break;
      }

      // we proceed along the incoming line
      // in any case
      if (current->pData().partonId.state == ClusteringParticleState::initial &&
	  emitter->pData().partonId.state == ClusteringParticleState::initial) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "initial line originated from initial line. continuing"  << endl;
#endif
	current = emitter;
	continue;
      }

      // for a glue, go back until we find a harder glue or a quark

      // fs glue

      if (current->pData().partonId.state == ClusteringParticleState::final &&
	  current->pData().partonId.PDGId == 21) {

	// glue coming from quark
	if (abs(emitter->pData().partonId.PDGId) < 7) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	  generator()->log() << "glue from quark. stop."  << endl;
#endif
	  break;
	}

	// glue coming from a glue
	if (emitter->pData().partonId.PDGId == 21 && partner->pData().partonId.PDGId == 21) {

	  bool harder = false;
	  
	  if (!hadronic)
	    harder = partner->momentum().e() > current->momentum().e();
	  else
	    harder = partner->momentum().vect().perp2() > current->momentum().vect().perp2();
	    
	  if (harder) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	    generator()->log() << "harder glue found. stop."  << endl;
#endif
	    break;
	  }
	  else {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	    generator()->log() << "softer glue emitted. continue"  << endl;
#endif
	    current = emitter;
	    continue;
	  }
	}

      }

      // fs quark      

      if (current->pData().partonId.state == ClusteringParticleState::final &&
	  abs(current->pData().partonId.PDGId) < 7) {

	// quark coming from quark
	if (abs(emitter->pData().partonId.PDGId) < 7) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	  generator()->log() << "emission from quark line. continue"  << endl;
#endif
	  current = emitter;
	  continue;
	}

	// quark coming from gluon -> proceed with the gluon
	if (emitter->pData().partonId.PDGId == 21) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "quark from gluon. continuing with the glue"  << endl;
#endif
	  current = emitter;
	  continue;
	}

      }



    }

    // set the shower scale form current's production scale

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "production point for " << *p << " found." << endl
		       << "shower scale / GeV2 = " << current->productionScale()/GeV2 << endl;
#endif

    (**p).showerScale(current->productionScale());

  }

}

