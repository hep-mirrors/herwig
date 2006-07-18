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

MultiplicityCount::MultiplicityCount() 
{
  _data[0]   = MultiplicityInfo(20.924,0.117,lightMeson);
  _data[22]  = MultiplicityInfo(21.27, 0.60,lightMeson);
  _data[111] = MultiplicityInfo(9.59,0.33,lightMeson);
  _data[113] = MultiplicityInfo(1.295,0.125,lightMeson);
  _data[211] = MultiplicityInfo(17.04,0.25,lightMeson);
  _data[213] = MultiplicityInfo(2.40,0.43,lightMeson);
  _data[221] = MultiplicityInfo(0.956,0.049,lightMeson);
  _data[223] = MultiplicityInfo(1.083,0.088,lightMeson);
  _data[225] = MultiplicityInfo(0.168,0.021,other);
  _data[311] = MultiplicityInfo(2.027,0.025,lightMeson);
  _data[313] = MultiplicityInfo(0.761,0.032,strangeMeson);
  _data[315] = MultiplicityInfo(0.106,0.060,strangeMeson);
  _data[321] = MultiplicityInfo(2.319,0.079,strangeMeson);
  _data[323] = MultiplicityInfo(0.731,0.058,strangeMeson);
  _data[331] = MultiplicityInfo(0.152,0.030,lightMeson);
  _data[333] = MultiplicityInfo(0.097,0.007,strangeMeson);
  _data[335] = MultiplicityInfo(0.020,0.008,other);
  _data[411] = MultiplicityInfo(0.184,0.018,other);
  _data[413] = MultiplicityInfo(0.182,0.009,other);
  _data[421] = MultiplicityInfo(0.473,0.026,other);
  _data[431] = MultiplicityInfo(0.129,0.013,other);
  _data[433] = MultiplicityInfo(0.096,0.046,other);
  _data[443] = MultiplicityInfo(0.00544,0.00029,other);
  _data[2212] = MultiplicityInfo(0.991,0.054,lightBaryon);
  _data[2112] = MultiplicityInfo(0.991,0.054,lightBaryon);
  _data[2224] = MultiplicityInfo(0.088,0.034,lightBaryon);
  _data[2214] = MultiplicityInfo(0.000,0.000,lightBaryon);
  _data[2114] = MultiplicityInfo(0.000,0.000,lightBaryon);
  _data[3112] = MultiplicityInfo(0.083,0.011,lightBaryon);
  _data[3122] = MultiplicityInfo(0.373,0.008,lightBaryon);
  _data[3212] = MultiplicityInfo(0.074,0.009,lightBaryon);
  _data[3222] = MultiplicityInfo(0.099,0.015,lightBaryon);
  _data[3224] = MultiplicityInfo(0.0471,0.0046,lightBaryon);
  _data[3312] = MultiplicityInfo(0.0262,0.0010,lightBaryon);
  _data[3324] = MultiplicityInfo(0.0058,0.0010,lightBaryon);
  _data[3334] = MultiplicityInfo(0.00125,0.00024,lightBaryon);
  _data[4122] = MultiplicityInfo(0.077,0.016,other);
  _data[100443 ] = MultiplicityInfo(0.00229,0.00041,other);
  _data[9000211] = MultiplicityInfo(0.27,0.11,other);
  _data[10221  ] = MultiplicityInfo(0.142,0.011,other);
}

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

  StepVector steps = event->primaryCollision()->steps();

  particles.clear();
  steps[0]->select(inserter(particles), ThePEG::AllSelector());
  
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = (*it)->id();
    _collisioncount.insert(make_pair(ID,0));
    ++_collisioncount[ID];
  }

  // =======


  particles.clear();

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


NoPIOClassDescription<MultiplicityCount> MultiplicityCount::initMultiplicityCount;
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

