// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BasicConsistency class.
//

#include "BasicConsistency.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BasicConsistency.tcc"
#endif
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG;

void BasicConsistency::analyze(tEventPtr event, long, int, int) {
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));
    
  int charge(-event->incoming().first->dataPtr()->iCharge()
	     -event->incoming().second->dataPtr()->iCharge());

  Lorentz5Momentum 
    ptotal(-event->incoming().first->momentum()
	   -event->incoming().second->momentum());

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if(abs((*it)->id()) < 9) {
      cerr << "Had quarks in final state in event " 
	   << event->number()  << '\n';
      generator()->log() << "Had quarks in final state in event " 
			 << event->number()  << '\n'
			 << *event;
    }
    else if((**it).id()==ExtraParticleID::Cluster)
      {
	cerr << "Had clusters in final state in event " 
	     << event->number()  << '\n';
	generator()->log() << "Had clusters in final state in event " 
			   << event->number()  << '\n'
			   << *event;
      }
    charge += (*it)->dataPtr()->iCharge();
    ptotal += (*it)->momentum();
  }
  
  if (charge != 0) {
    cerr << "\nCharge imbalance by " << charge 
	 << "in event " << event->number()  << '\n';// << *event;
    generator()->log() << "Charge imbalance by " << charge 
		       << "in event " << event->number()  << '\n' 
		       << *event;
  }
  if (ptotal.mag() > 5.*MeV || abs(ptotal.t()) > 5.*MeV) {
    cerr << "\nMomentum imbalance by " << ptotal/GeV 
	 << " GeV in event " << event->number() << '\n';// << *event;
    generator()->log() <<"\nMomentum imbalance by " << ptotal/GeV 
		       << " GeV in event " << event->number() << '\n' 
		       << *event;
  }
  
  if (ptotal.mag() > _epsmom)
    _epsmom = ptotal.mag();

  if (abs(ptotal.t()) > _epsmom)
    _epsmom=abs(ptotal.t());

}

void BasicConsistency::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void BasicConsistency::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<BasicConsistency> BasicConsistency::initBasicConsistency;
// Definition of the static class description member.

void BasicConsistency::Init() {

  static ClassDocumentation<BasicConsistency> documentation
    ("The BasicConsistency analysis handler checks for momentum and charge conservation.");

}

