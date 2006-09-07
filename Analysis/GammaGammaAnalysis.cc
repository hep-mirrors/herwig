// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaAnalysis class.
//

#include "GammaGammaAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GammaGammaAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GammaGammaAnalysis::~GammaGammaAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
    //Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
    //tmp += particle->children()[1]->next()->momentum();
    //tmp.rescaleMass();
    //return tmp;
  } 
}


void GammaGammaAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum p1, p2, ppair;  
  bool foundphotons = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id()==ParticleID::gamma) {
      // find the two hardest photons in the event
      if( getMomentum(*it).perp() > p2.perp() ) {
	if (getMomentum(*it).perp() > p1.perp()) {
	  p2 = p1;
	  p1 = getMomentum(*it);
	} else {
	  p2 = getMomentum(*it);
	}
      }
    }
  }

  //  cerr << "E1 = " <<  p1.e()/GeV << ", E1 = " <<  p1.e()/GeV << "\n";
  if (p1.perp()/GeV > 0 && p2.perp()/GeV > 0) foundphotons = true;

  ppair = p1 + p2;
  if (foundphotons) {
    _ptharder += p1.perp()/GeV;
    _ptsofter += p2.perp()/GeV;
    _ptpair += ppair.perp()/GeV;
    _Eharder += p1.e()/GeV;
    _Esofter += p2.e()/GeV;
    _Epair += ppair.e()/GeV;
    _rapharder += p1.rapidity();
    _rapsofter += p2.rapidity();
    _rappair += ppair.rapidity();
    _phiharder += p1.phi();
    _phisofter += p2.phi();
    _deltaphi += (p2.vect()).deltaPhi(p1.vect());
    _mpair += ppair.m()/GeV;
  } else {
    cerr << "Analysis/GammaGammaAnalysis: Found no hard photon in event " 
	 << event->number()  << ".\n";
    generator()->log() << "Analysis/GammaGammaAnalysis: " 
		       << "Found no hard photon in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

LorentzRotation GammaGammaAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void GammaGammaAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void GammaGammaAnalysis::analyze(tPPtr) {}

void GammaGammaAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void GammaGammaAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<GammaGammaAnalysis> GammaGammaAnalysis::initGammaGammaAnalysis;
// Definition of the static class description member.

void GammaGammaAnalysis::Init() {

  static ClassDocumentation<GammaGammaAnalysis> documentation
    ("There is no documentation for the GammaGammaAnalysis class");

}

