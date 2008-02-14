// -*- C++ -*-
//
// SplittingGenerator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingGenerator class.
//

#include "SplittingGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/Repository.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _qcdinteractionMode << _qedinteractionMode << _ewkinteractionMode
     << _isr_Mode << _isr_qcdMode << _isr_qedMode << _isr_ewkMode     
     << _fsr_Mode << _fsr_qcdMode << _fsr_qedMode << _fsr_ewkMode
     << _bbranchings << _fbranchings;
}

void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_qcdinteractionMode >> _qedinteractionMode >> _ewkinteractionMode
     >>	_isr_Mode >> _isr_qcdMode >> _isr_qedMode >> _isr_ewkMode 
     >> _fsr_Mode >> _fsr_qcdMode >> _fsr_qedMode >> _fsr_ewkMode
     >> _bbranchings >> _fbranchings;
}

ClassDescription<SplittingGenerator> SplittingGenerator::initSplittingGenerator;
// Definition of the static class description member.

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

  static Switch<SplittingGenerator, bool> interfaceQCDinteractionMode
    ("QCDinteractions",
     "Should shower include QCD interactions",
     &SplittingGenerator::_qcdinteractionMode, 1, false, false);
  static SwitchOption interfaceQCDinteractionMode0
    (interfaceQCDinteractionMode,"No","QCD interaction is Off", 0);
  static SwitchOption interfaceQCDinteractionMode1
    (interfaceQCDinteractionMode,"Yes","QCD interaction is On", 1);

  static Switch<SplittingGenerator, bool> interfaceQEDinteractionMode
    ("QEDinteractions",
     "Should shower include QED interactions",
     &SplittingGenerator::_qedinteractionMode, 0, false, false);
  static SwitchOption interfaceQEDinteractionMode0
    (interfaceQEDinteractionMode,"No","QED interaction is Off", 0);
  static SwitchOption interfaceQEDinteractionMode1
    (interfaceQEDinteractionMode,"Yes","QED interaction is On", 1);

  static Switch<SplittingGenerator, bool> interfaceEWKinteractionMode
    ("EWKinteractions",
     "Should shower include EWK interactions",
     &SplittingGenerator::_ewkinteractionMode, 0, false, false);
  static SwitchOption interfaceEWKinteractionMode0
    (interfaceEWKinteractionMode,"No","EWK interaction is OFF", 0);
  static SwitchOption interfaceEWKinteractionMode1
    (interfaceEWKinteractionMode,"Yes","EWK interaction is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISRMode
    ("ISR",
     "Include initial-state radiation?",
     &SplittingGenerator::_isr_Mode, 1, false, false);
  static SwitchOption interfaceISRMode0
    (interfaceISRMode,"No","ISR (Initial State Radiation) is OFF", 0);
  static SwitchOption interfaceISRMode1
    (interfaceISRMode,"Yes","ISR (Initial State Radiation) is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_qcdMode
    ("ISR_QCD",
     "Include initial-state QCD radiation?",
     &SplittingGenerator::_isr_qcdMode, 1, false, false);
  static SwitchOption interfaceISR_qcdMode0
    (interfaceISR_qcdMode,"No","QCD ISR is OFF", 0);
  static SwitchOption interfaceISR_qcdMode1
    (interfaceISR_qcdMode,"Yes","QCD ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_qedMode
    ("ISR_QED",
     "Include initial-state QED radiation?",
     &SplittingGenerator::_isr_qedMode, 1, false, false);
  static SwitchOption interfaceISR_qedMode0
    (interfaceISR_qedMode,"No","QED ISR is OFF", 0);
  static SwitchOption interfaceISR_qedMode1
    (interfaceISR_qedMode,"Yes","QED ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_ewkMode
    ("ISR_EWK",
     "Include initial-state EWK radiation?",
     &SplittingGenerator::_isr_ewkMode, 1, false, false);
  static SwitchOption interfaceISR_ewkMode0
    (interfaceISR_ewkMode,"No","EWK ISR is OFF", 0);
  static SwitchOption interfaceISR_ewkMode1
    (interfaceISR_ewkMode,"Yes","EWK ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSRMode
    ("FSR",
     "Include final-state radiation?",
     &SplittingGenerator::_fsr_Mode, 1, false, false);
  static SwitchOption interfaceFSRMode0
    (interfaceFSRMode,"No","FSR (Final State Radiation) is OFF", 0);
  static SwitchOption interfaceFSRMode1
    (interfaceFSRMode,"Yes","FSR (Final State Radiation) is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_qcdMode
    ("FSR_QCD",
     "Include final-state QCD radiation?",
     &SplittingGenerator::_fsr_qcdMode, 1, false, false);
  static SwitchOption interfaceFSR_qcdMode0
    (interfaceFSR_qcdMode,"No","QCD FSR is OFF", 0);
  static SwitchOption interfaceFSR_qcdMode1
    (interfaceFSR_qcdMode,"Yes","QCD FSR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_qedMode
    ("FSR_QED",
     "Include final-state QED radiation?",
     &SplittingGenerator::_fsr_qedMode, 1, false, false);
  static SwitchOption interfaceFSR_qedMode0
    (interfaceFSR_qedMode,"No","QED FSR is OFF", 0);
  static SwitchOption interfaceFSR_qedMode1
    (interfaceFSR_qedMode,"Yes","QED FSR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_ewkMode
    ("FSR_EWK",
     "Include final-state EWK radiation?",
     &SplittingGenerator::_fsr_ewkMode, 1, false, false);

  static SwitchOption interfaceFSR_ewkMode0
    (interfaceFSR_ewkMode,"No","EWK FSR is OFF", 0);

  static SwitchOption interfaceFSR_ewkMode1
    (interfaceFSR_ewkMode,"Yes","EWK FSR is ON", 1);

  static Command<SplittingGenerator> interfaceAddSplitting
    ("AddFinalSplitting",
     "Adds another splitting to the list of splittings considered "
     "in the shower. Command is a->b,c; Sudakov",
     &SplittingGenerator::addFinalSplitting);
  static Command<SplittingGenerator> interfaceAddInitialSplitting
    ("AddInitialSplitting",
     "Adds another splitting to the list of initial splittings to consider "
     "in the shower. Command is a->b,c; Sudakov. Here the particle a is the "
     "particle that is PRODUCED by the splitting. b is the initial state "
     "particle that is splitting in the shower.",
     &SplittingGenerator::addInitialSplitting);
}

string SplittingGenerator::addSplitting(string arg, bool final) {
  string partons = StringUtils::car(arg);
  string sudakov = StringUtils::cdr(arg);
  vector<tPDPtr> products;
  string::size_type next = partons.find("->");
  if(next == string::npos) 
    return "Error: Invalid string for splitting " + arg;
  if(partons.find(';') == string::npos) 
    return "Error: Invalid string for splitting " + arg;
  tPDPtr parent = Repository::findParticle(partons.substr(0,next));
  partons = partons.substr(next+2);
  do {
    next = min(partons.find(','), partons.find(';'));
    tPDPtr pdp = Repository::findParticle(partons.substr(0,next));
    partons = partons.substr(next+1);
    if(pdp) products.push_back(pdp);
    else return "Error: Could not create splitting from " + arg;
  } while(partons[0] != ';' && partons.size());
  SudakovPtr s;
  s = dynamic_ptr_cast<SudakovPtr>(Repository::TraceObject(sudakov));
  if(!s) return "Error: Could not load Sudakov " + sudakov + '\n';
  IdList ids;
  ids.push_back(parent->id());
  for(vector<tPDPtr>::iterator it = products.begin(); it!=products.end(); ++it)
    ids.push_back((*it)->id());
  // check splitting can handle this
  if(!s->splittingFn()->accept(ids)) 
    return "Error: Sudakov " + sudakov + "can't handle particles\n";
  // add to map
  addToMap(ids,s,final);
  return "";
}

void SplittingGenerator::addToMap(const IdList &ids, const SudakovPtr &s, bool final) {
  if(isISRadiationON(s->splittingFn()->interactionType()) && !final) {
    _bbranchings.insert(BranchingInsert(ids[1],BranchingElement(s,ids)));
    s->addSplitting(ids);
  }
  if(isFSRadiationON(s->splittingFn()->interactionType()) && final) {
    _fbranchings.insert(BranchingInsert(ids[0],BranchingElement(s,ids)));
    s->addSplitting(ids);
  }
}

Branching SplittingGenerator::chooseForwardBranching(ShowerParticle &particle,
						     double enhance) const {
  Energy newQ = Energy();
  ShoKinPtr kinematics = ShoKinPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), IdList());
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) {
    ShowerIndex::InteractionType i = cit->second.first->interactionType();
    // check size of scales beforehand...
    ShoKinPtr newKin= cit->second.first->
      generateNextTimeBranching(particle.evolutionScales()[i], 
				cit->second.second,particle.id()!=cit->first,
				enhance);
    if(!newKin) continue;
    // select highest scale 
    if(newKin->scale() > newQ && 
       newKin->scale() <= particle.evolutionScales()[i]) {
      kinematics=newKin;
      newQ = newKin->scale();
      ids = cit->second.second;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList());
  // If a branching has been selected initialize it
  kinematics->initialize(particle,PPtr());
  // and return it
  return Branching(kinematics, ids);
}

Branching SplittingGenerator::chooseDecayBranching(ShowerParticle &particle,
						   vector<Energy> stoppingScale,
						   Energy minmass,
						   double enhance) const {
  Energy newQ = Constants::MaxEnergy;
  ShoKinPtr kinematics;
  IdList ids;
  // First, find the eventual branching, corresponding to the lowest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), IdList());
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit)  {
    ShowerIndex::InteractionType i = cit->second.first->interactionType();
    ShoKinPtr newKin;
    if(particle.evolutionScales()[i] < stoppingScale[i]) 
      newKin = cit->second.first->
	generateNextDecayBranching(particle.evolutionScales()[i],
				   stoppingScale[i],minmass,
				   cit->second.second,
				   particle.id()!=cit->first,enhance);
    if(!newKin) continue;
    if(newKin->scale() < newQ && newKin->scale() > particle.evolutionScales()[i]) {
      newQ = newKin->scale();
      ids = cit->second.second;
      kinematics=newKin;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList());
  // initialize the branching
  kinematics->initialize(particle,PPtr());
  // and return it
  return Branching(kinematics, ids);
}

Branching SplittingGenerator::
chooseBackwardBranching(ShowerParticle &particle,PPtr beamparticle,
			double enhance,
			Ptr<BeamParticleData>::transient_const_pointer beam) const {
  Energy newQ=Energy();
  ShoKinPtr kinematics=ShoKinPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.id());
  // if no possible branching return
  if(_bbranchings.find(index) == _bbranchings.end())
    return Branching(ShoKinPtr(), IdList());
  // select the branching
  for(BranchingList::const_iterator cit = _bbranchings.lower_bound(index); 
      cit != _bbranchings.upper_bound(index); ++cit ) {
    ShowerIndex::InteractionType i = cit->second.first->interactionType(); 
    ShoKinPtr newKin=cit->second.first->
      generateNextSpaceBranching(particle.evolutionScales()[i],
				 cit->second.second, particle.x(),
				 particle.id()!=cit->first,enhance,beam);
    if(!newKin) continue;
    if(newKin->scale() > newQ) {
      newQ = newKin->scale();
      kinematics=newKin;
      ids = cit->second.second;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics) return Branching(ShoKinPtr(), IdList());
  // initialize the ShowerKinematics 
  // and return it
  kinematics->initialize(particle,beamparticle);
  // return the answer
  return Branching(kinematics, ids);
}

