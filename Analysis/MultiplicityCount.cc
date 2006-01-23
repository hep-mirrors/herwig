// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiplicityCount class.
//

#include "MultiplicityCount.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/EventRecord/Event.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MultiplicityCount.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

MultiplicityCount::~MultiplicityCount() {}

void MultiplicityCount::analyze(tEventPtr event, long ieve, int loop, int state) {

  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  map <long,long> eventcount;
  eventcount.insert(make_pair(0,0));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {

    if((*it)->dataPtr()->charged()) 
      ++eventcount[0];
    
    long ID = abs( (*it)->id() );
    _finalstatecount.insert(make_pair(ID,0));
    ++_finalstatecount[ID];
  }

  // ========

  particles.clear();

  StepVector steps = event->primaryCollision()->steps();
  for ( StepVector::const_iterator it = steps.begin()+2;
	it != steps.end(); ++it ) {
    (**it).select(inserter(particles), ThePEG::AllSelector());
  }

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = abs( (*it)->id() );
    
    if (_data.find(ID) != _data.end()) {
      eventcount.insert(make_pair(ID,0));
      ++eventcount[ID];
    }
  }
  
  for(map<long,long>::const_iterator it = eventcount.begin();
      it != eventcount.end(); ++it) {
    _data[it->first].actualCount += it->second;
    _data[it->first].sumofsquares += sqr(double(it->second));
  }
}

LorentzRotation MultiplicityCount::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void MultiplicityCount::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void MultiplicityCount::analyze(tPPtr) {}

void MultiplicityCount::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void MultiplicityCount::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<MultiplicityCount> MultiplicityCount::initMultiplicityCount;
// Definition of the static class description member.

void MultiplicityCount::Init() {

  static ParVector<MultiplicityCount,long> interfaceparticlecodes
    ("ParticleCodes",
     "The PDG code for the particles",
     &MultiplicityCount::_particlecodes,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<MultiplicityCount,double> interfaceMultiplicity
    ("Multiplicity",
     "The multiplicity for the particle",
     &MultiplicityCount::_multiplicity,
     0, 0, 0, 0., 1000., false, false, true);

  static ParVector<MultiplicityCount,double> interfaceError
    ("Error",
     "The error on the multiplicity for the particle",
     &MultiplicityCount::_error,
     0, 0, 0, 0., 1000., false, false, true);

  static ParVector<MultiplicityCount,unsigned int> interfaceSpecies
    ("Species",
     "The type of particle",
     &MultiplicityCount::_species,
     0, 0, other, 0, other, false, false, true);

}

