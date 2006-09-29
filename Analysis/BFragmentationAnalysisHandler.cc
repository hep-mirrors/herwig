// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BFragmentationAnalysisHandler class.
//

#include "BFragmentationAnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig++/Utilities/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BFragmentationAnalysisHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

BFragmentationAnalysisHandler::~BFragmentationAnalysisHandler() {}

void BFragmentationAnalysisHandler::analyze(tEventPtr event, long,
					    int loop, int state) 
{
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // extract the weakly decaying B hadrons using set to avoid double counting
  set<PPtr> allParticles;
  event->select(inserter(allParticles),WeakBHadronSelector());
  // convert to vector
  tPVector particles(allParticles.begin(),allParticles.end());
  // numerator
  _emax = 0.5*generator()->maximumCMEnergy();   
  analyze(particles); 
}

LorentzRotation BFragmentationAnalysisHandler::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void BFragmentationAnalysisHandler::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void BFragmentationAnalysisHandler::analyze(tPPtr part) 
{
  *_fragBxE  += part->momentum().e()/_emax;
  *_fragBxEa += part->momentum().e()/_emax;
}

void BFragmentationAnalysisHandler::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void BFragmentationAnalysisHandler::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<BFragmentationAnalysisHandler> BFragmentationAnalysisHandler::initBFragmentationAnalysisHandler;
// Definition of the static class description member.

void BFragmentationAnalysisHandler::Init() {

  static ClassDocumentation<BFragmentationAnalysisHandler> documentation
    ("The BFragmentationAnalysisHandler class performs analysis"
     " of the B fragmentation function");
  
}
