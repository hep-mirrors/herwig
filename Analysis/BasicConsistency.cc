// -*- C++ -*-
//
// BasicConsistency.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BasicConsistency class.
//

#include "BasicConsistency.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;
using namespace ThePEG;

BasicConsistency::BasicConsistency() 
  : _epsmom(ZERO),_checkquark(true), _checkcharge(true),
    _checkcluster(true), _checkBR(true)
{}

IBPtr BasicConsistency::clone() const {
  return new_ptr(*this);
}

IBPtr BasicConsistency::fullclone() const {
  return new_ptr(*this);
}

void BasicConsistency::analyze(tEventPtr event, long, int, int) {
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));
    
  int charge(-event->incoming().first->dataPtr()->iCharge()
	     -event->incoming().second->dataPtr()->iCharge());

  Lorentz5Momentum 
    ptotal(-event->incoming().first->momentum()
	   -event->incoming().second->momentum());

  const Energy beamenergy = ptotal.mag();

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if (_checkquark && (*it)->coloured()) {
      cerr << "Had quarks in final state in event " 
	   << event->number()  
	   << '\n';
      generator()->log() << "Had quarks in final state in event " 
			 << event->number()  
			 << '\n'
			 << *event;
    }
    else if( _checkcluster && (**it).id()==ExtraParticleID::Cluster) {
      cerr << "Had clusters in final state in event " 
	   << event->number()  
	   << '\n';
      generator()->log() << "Had clusters in final state in event " 
			 << event->number()  
			 << '\n'
			 << *event;
    }
    charge += (*it)->dataPtr()->iCharge();
    ptotal += (*it)->momentum();
    bool problem=false;
    LorentzDistance test;
    for(unsigned int ix=0;ix<5;++ix) {
      switch (ix) {
      case 0:
	test = (*it)->vertex();
	break;
      case 1:
	test = (*it)->labVertex();
	break;
      case 2:
	test = (*it)->decayVertex();
	break;
      case 3:
	test = (*it)->labDecayVertex();
	break;
      case 4:
	test = (*it)->lifeLength();
	break;
      }
      problem |= 
	isnan(test.x()/mm) || isnan(test.y()/mm) ||
	isnan(test.z()/mm) || isnan(test.t()/mm) ||
	isinf(test.x()/mm) || isinf(test.y()/mm) ||
	isinf(test.z()/mm) || isinf(test.t()/mm);
    }
    if(problem) {
      generator()->log() << "Problem with position of " << **it << "\n"
			 << (*it)->vertex()/mm << "\n"
			 << (*it)->labVertex()/mm << "\n"
			 << (*it)->decayVertex()/mm << "\n"
			 << (*it)->labDecayVertex()/mm << "\n"
			 << (*it)->lifeLength()/mm << "\n"; 
    }
  }
  
  if ( _checkcharge && charge != 0 ) {
    cerr << "\nCharge imbalance by " 
	 << charge 
	 << "in event " 
	 << event->number()  
	 << '\n';
    generator()->log() << "Charge imbalance by " 
		       << charge 
		       << "in event " 
		       << event->number()  
		       << '\n' 
		       << *event;
  }

  Energy mag = ptotal.mag();
  Energy ee  = ptotal.e();

  if (isnan(mag/MeV)) {
    cerr << "\nMomentum is 'nan'; " << ptotal/MeV 
	 << " MeV in event " << event->number() << '\n';
    generator()->log() <<"\nMomentum is 'nan'; " << ptotal/MeV 
		       << " MeV in event " << event->number() << '\n' 
		       << *event;
  }

  const Energy epsilonmax = 1.0e-5 * beamenergy;

  if (mag > epsilonmax || abs(ee) > epsilonmax) {
    cerr << "\nMomentum imbalance by " << ptotal/MeV 
	 << " MeV in event " << event->number() << '\n';
    generator()->log() <<"\nMomentum imbalance by " << ptotal/MeV 
		       << " MeV in event " << event->number() << '\n' 
		       << *event;
  }

  if (mag > _epsmom)
    _epsmom = mag;

  if (abs(ee) > _epsmom)
    _epsmom = abs(ee);

  if (abs(ptotal.x()) > _epsmom)
    _epsmom = abs(ptotal.x());

  if (abs(ptotal.y()) > _epsmom)
    _epsmom = abs(ptotal.y());

  if (abs(ptotal.z()) > _epsmom)
    _epsmom = abs(ptotal.z());

  particles.clear();

  event->select(inserter(particles), ThePEG::AllSelector());
  bool output=false;
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    bool problem=false;
    LorentzDistance test;
    for(unsigned int ix=0;ix<5;++ix) {
      switch (ix) {
      case 0:
	test = (*it)->vertex();
	break;
      case 1:
	test = (*it)->labVertex();
	break;
      case 2:
	test = (*it)->decayVertex();
	break;
      case 3:
	test = (*it)->labDecayVertex();
	break;
      case 4:
	test = (*it)->lifeLength();
	break;
      }
      problem |= isnan(test.mag2()/mm/mm) || isinf(test.mag2()/mm/mm);
    }
    if(problem) {
      generator()->log() << "Problem with position of " << **it << "\n"
			 << (*it)->vertex()/mm << "\n"
			 << (*it)->labVertex()/mm << "\n"
			 << (*it)->decayVertex()/mm << "\n"
			 << (*it)->labDecayVertex()/mm << "\n"
			 << (*it)->lifeLength()/mm << "\n"; 
      output=true;
    }
  }
  if(output) generator()->log() << *event;
}

void BasicConsistency::persistentOutput(PersistentOStream & os) const {
  os << _checkquark << _checkcharge << _checkcluster << _checkBR;
}

void BasicConsistency::persistentInput(PersistentIStream & is, int) {
  is >> _checkquark >> _checkcharge >> _checkcluster >> _checkBR;
}

ClassDescription<BasicConsistency> BasicConsistency::initBasicConsistency;
// Definition of the static class description member.

void BasicConsistency::Init() {

  static ClassDocumentation<BasicConsistency> documentation
    ("The BasicConsistency analysis handler checks for"
     " momentum and charge conservation.");

  static Switch<BasicConsistency,bool> interfaceCheckQuark
    ("CheckQuark",
     "Check whether there are quarks in the final state",
     &BasicConsistency::_checkquark, true, false, false);
  static SwitchOption interfaceCheckQuarkCheck
    (interfaceCheckQuark,
     "Yes",
     "Check for quarks",
     true);
  static SwitchOption interfaceCheckQuarkNoCheck
    (interfaceCheckQuark,
     "No",
     "Don't check for quarks",
     false);

  static Switch<BasicConsistency,bool> interfaceCheckCharge
    ("CheckCharge",
     "Check whether charge is conserved",
     &BasicConsistency::_checkcharge, true, false, false);
  static SwitchOption interfaceCheckChargeCheck
    (interfaceCheckCharge,
     "Yes",
     "Check charge conservation",
     true);
  static SwitchOption interfaceCheckChargeNoCheck
    (interfaceCheckCharge,
     "No",
     "Don't check charge conservation",
     false);

  static Switch<BasicConsistency,bool> interfaceCheckCluster
    ("CheckCluster",
     "Check whether there are clusters in the final state",
     &BasicConsistency::_checkcluster, true, false, false);
  static SwitchOption interfaceCheckClusterCheck
    (interfaceCheckCluster,
     "Yes",
     "Check for clusters",
     true);
  static SwitchOption interfaceCheckClusterNoCheck
    (interfaceCheckCluster,
     "No",
     "Don't check for clusters",
     false);

  static Switch<BasicConsistency,bool> interfaceCheckBranchingRatios
    ("CheckBranchingRatios",
     "Check whether the branching ratios of the particles add up to one.",
     &BasicConsistency::_checkBR, true, false, false);
  static SwitchOption interfaceCheckBranchingRatiosYes
    (interfaceCheckBranchingRatios,
     "Yes",
     "Perform the check",
     true);
  static SwitchOption interfaceCheckBranchingRatiosNo
    (interfaceCheckBranchingRatios,
     "No",
     "Don't perform the check",
     false);

}

void BasicConsistency::dofinish() {
  AnalysisHandler::dofinish();
  cout << "\nBasicConsistency: maximum 4-momentum violation: " 
       << _epsmom/MeV << " MeV\n";
}

void BasicConsistency::doinitrun() {
  AnalysisHandler::doinitrun();
  static double eps=1e-12;
  for(ParticleMap::const_iterator it=generator()->particles().begin();
      it!=generator()->particles().end();++it) {
    if(it->second->stable()) continue;
    double total(0.);
    for(DecaySet::const_iterator dit=it->second->decayModes().begin();
	dit!=it->second->decayModes().end();++dit) {
      if((**dit).on()) total +=(**dit).brat();
    }
    if(abs(total-1.)>eps) {
      cerr << "Warning: Total BR for " 
	   << it->second->PDGName() 
	   << " does not add up to 1. sum = " << total << "\n";
    }
  }
}
