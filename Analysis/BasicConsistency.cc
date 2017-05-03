// -*- C++ -*-
//
// BasicConsistency.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BasicConsistency class.
//

#include "BasicConsistency.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Utilities/ColourOutput.h"

using namespace Herwig;
using namespace ThePEG;

BasicConsistency::BasicConsistency() 
  : _epsmom(ZERO),_checkquark(true), _checkcharge(true),
    _checkcluster(true), _checkBR(true),
    _absolutemomentumtolerance(1*MeV), _relativemomentumtolerance(1e-5)
{}

IBPtr BasicConsistency::clone() const {
  return new_ptr(*this);
}

IBPtr BasicConsistency::fullclone() const {
  return new_ptr(*this);
}

void BasicConsistency::analyze(tEventPtr event, long, int, int) {
  bool writeEvent=false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));
    
  int charge(-event->incoming().first->dataPtr()->iCharge()
	     -event->incoming().second->dataPtr()->iCharge());

  Lorentz5Momentum 
    ptotal(-event->incoming().first->momentum()
	   -event->incoming().second->momentum());

  const Energy beamenergy = ptotal.m();

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if (_checkquark && (*it)->coloured()) {
      cerr << "Had quarks in final state in event " 
	   << event->number()  
	   << '\n';
      generator()->log() << "Had quarks in final state in event " 
			 << event->number() << '\n';
      writeEvent = true;
    }
    else if( _checkcluster && (**it).id()==ParticleID::Cluster) {
      cerr << "Had clusters in final state in event " 
	   << event->number()  
	   << '\n';
      generator()->log() << "Had clusters in final state in event " 
			 << event->number()  << '\n';
      writeEvent = true;
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
      problem |= ! ( isfinite(double(test.x()/mm)) && 
                     isfinite(double(test.y()/mm)) && 
                     isfinite(double(test.z()/mm)) && 
                     isfinite(double(test.t()/mm)) );
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
		       << event->number() << '\n';
    writeEvent = true;
  }

  Energy mag = ptotal.m();
  Energy ee  = ptotal.e();

  if (std::isnan(double(mag/MeV))) {
    cerr << "\nMomentum is 'nan'; " << ptotal/MeV 
	 << " MeV in event " << event->number() << '\n';
    generator()->log() <<"\nMomentum is 'nan'; " << ptotal/MeV 
		       << " MeV in event " << event->number() << '\n';
    writeEvent = true;
  }

  const Energy epsilonmax = max( _absolutemomentumtolerance,
				 _relativemomentumtolerance * beamenergy );

  if (abs(mag) > epsilonmax || abs(ee) > epsilonmax) {
    cerr << "\nMomentum imbalance by " << ptotal/MeV 
	 << " MeV in event " << event->number() << '\n';
    generator()->log() <<"\nMomentum imbalance by " << ptotal/MeV 
		       << " MeV in event " << event->number() << '\n';
    writeEvent = true;
  }

  if (abs(mag) > _epsmom)
    _epsmom = abs(mag);

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
      problem |= ( ! isfinite(double(test.m2()/mm/mm)) );
    }
    if(problem) {
      generator()->log() << "Problem with position of " << **it << "\n"
			 << (*it)->vertex()/mm << "\n"
			 << (*it)->labVertex()/mm << "\n"
			 << (*it)->decayVertex()/mm << "\n"
			 << (*it)->labDecayVertex()/mm << "\n"
			 << (*it)->lifeLength()/mm << "\n"; 
      writeEvent=true;
    }
  }
  if(writeEvent) generator()->log() << *event;
}

void BasicConsistency::persistentOutput(PersistentOStream & os) const {
  os << _checkquark << _checkcharge << _checkcluster << _checkBR
     << ounit(_absolutemomentumtolerance,MeV) << _relativemomentumtolerance;
}

void BasicConsistency::persistentInput(PersistentIStream & is, int) {
  is >> _checkquark >> _checkcharge >> _checkcluster >> _checkBR
     >> iunit(_absolutemomentumtolerance,MeV) >> _relativemomentumtolerance;
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

  static Parameter<BasicConsistency,Energy> interfaceAbsoluteMomentumTolerance
    ("AbsoluteMomentumTolerance",
     "The value of the momentum imbalance above which warnings are issued/MeV.\n"
     "Final tolerance is the larger of AbsoluteMomentumTolerance and\n"
     "RelativeMomentumTolerance*beam energy.",
     &BasicConsistency::_absolutemomentumtolerance, MeV, 1*MeV, ZERO, 1e10*GeV,
     false, false, true);

  static Parameter<BasicConsistency,double> interfaceRelativeMomentumTolerance
    ("RelativeMomentumTolerance",
     "The value of the momentum imbalance as a fraction of the beam energy\n"
     "above which warnings are issued.\n"
     "Final tolerance is the larger of AbsoluteMomentumTolerance and\n"
     "RelativeMomentumTolerance*beam energy.",
     &BasicConsistency::_relativemomentumtolerance, 1e-5, 0.0, 1.0,
     false, false, true);

}

void BasicConsistency::dofinish() {
  AnalysisHandler::dofinish();
  cout << "\nBasicConsistency: maximum 4-momentum violation: " 
       << ANSI::blue
       << _epsmom/MeV << " MeV\n"
       << ANSI::reset;
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
