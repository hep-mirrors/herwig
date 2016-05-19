// -*- C++ -*-
//
// SplittingGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/Repository.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Utilities/Rebinder.h"
#include <cassert>
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<SplittingGenerator,Interfaced>
describeSplittingGenerator ("Herwig::SplittingGenerator","");

IBPtr SplittingGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr SplittingGenerator::fullclone() const {
  return new_ptr(*this);
}


void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _isr_Mode << _fsr_Mode << _bbranchings << _fbranchings << _deTuning;
}

void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_isr_Mode >> _fsr_Mode >> _bbranchings >> _fbranchings >> _deTuning;
}

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

  static Switch<SplittingGenerator, bool> interfaceISRMode
    ("ISR",
     "Include initial-state radiation?",
     &SplittingGenerator::_isr_Mode, 1, false, false);
  static SwitchOption interfaceISRMode0
    (interfaceISRMode,"No","ISR (Initial State Radiation) is OFF", 0);
  static SwitchOption interfaceISRMode1
    (interfaceISRMode,"Yes","ISR (Initial State Radiation) is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSRMode
    ("FSR",
     "Include final-state radiation?",
     &SplittingGenerator::_fsr_Mode, 1, false, false);
  static SwitchOption interfaceFSRMode0
    (interfaceFSRMode,"No","FSR (Final State Radiation) is OFF", 0);
  static SwitchOption interfaceFSRMode1
    (interfaceFSRMode,"Yes","FSR (Final State Radiation) is ON", 1);

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

  static Command<SplittingGenerator> interfaceDeleteSplitting
    ("DeleteFinalSplitting",
     "Deletes a splitting from the list of splittings considered "
     "in the shower. Command is a->b,c; Sudakov",
     &SplittingGenerator::deleteFinalSplitting);

  static Command<SplittingGenerator> interfaceDeleteInitialSplitting
    ("DeleteInitialSplitting",
     "Deletes a splitting from the list of initial splittings to consider "
     "in the shower. Command is a->b,c; Sudakov. Here the particle a is the "
     "particle that is PRODUCED by the splitting. b is the initial state "
     "particle that is splitting in the shower.",
     &SplittingGenerator::deleteInitialSplitting);

  static Parameter<SplittingGenerator,double> interfaceDetuning
    ("Detuning",
     "The Detuning parameter to make the veto algorithm less efficient to improve the weight variations",
     &SplittingGenerator::_deTuning, 1.0, 1.0, 10.0,
     false, false, Interface::limited);

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

string SplittingGenerator::deleteSplitting(string arg, bool final) {
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
  // delete from map
  deleteFromMap(ids,s,final);
  return "";
}

void SplittingGenerator::addToMap(const IdList &ids, const SudakovPtr &s, bool final) {
  if(isISRadiationON() && !final) {
    _bbranchings.insert(BranchingInsert(ids[1],BranchingElement(s,ids)));
    s->addSplitting(ids);
  }
  if(isFSRadiationON() &&  final) {
    _fbranchings.insert(BranchingInsert(ids[0],BranchingElement(s,ids)));
    s->addSplitting(ids);
  }
}

void SplittingGenerator::deleteFromMap(const IdList &ids, 
				       const SudakovPtr &s, bool final) {
  if(isISRadiationON() && !final) {
    pair<BranchingList::iterator,BranchingList::iterator> 
      range = _bbranchings.equal_range(ids[1]);
    for(BranchingList::iterator it=range.first;
	it!=range.second&&it!=_bbranchings.end()&&it->first==ids[1];++it) {
      if(it->second.first==s&&it->second.second==ids) {
	BranchingList::iterator it2=it;
	--it;
	_bbranchings.erase(it2);
      }
    }
    s->removeSplitting(ids);
  }
  if(isFSRadiationON() &&  final) {
    pair<BranchingList::iterator,BranchingList::iterator> 
      range = _fbranchings.equal_range(ids[0]);
    for(BranchingList::iterator it=range.first;
	it!=range.second&&it!=_fbranchings.end()&&it->first==ids[0];++it) {
      if(it->second.first==s&&it->second.second==ids) {
	BranchingList::iterator it2 = it;
	--it;
	_fbranchings.erase(it2);
      }
    }
    s->removeSplitting(ids);
  }
}

Branching SplittingGenerator::chooseForwardBranching(ShowerParticle &particle,
						     double enhance,
						     ShowerInteraction::Type type) const {
  Energy newQ = ZERO;
  ShoKinPtr kinematics = ShoKinPtr();
  ShowerPartnerType::Type partnerType(ShowerPartnerType::Undefined);
  SudakovPtr sudakov   = SudakovPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if( _fbranchings.find(index) == _fbranchings.end() ) 
    return Branching(ShoKinPtr(), IdList(),SudakovPtr(),ShowerPartnerType::Undefined);
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) {
    // check either right interaction or doing both
    if(type != cit->second.first->interactionType() &&
       type != ShowerInteraction::Both ) continue;
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.first->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    ShowerPartnerType::Type type;
    // work out which starting scale we need
    if(cit->second.first->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      newKin = cit->second.first->
    	generateNextTimeBranching(startingScale,cit->second.second,
				  particle.id()!=cit->first,enhance,_deTuning);
    }
    else if(cit->second.first->interactionType()==ShowerInteraction::QCD) {
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.first->splittingFn()->colourStructure()==OctetOctetOctet) {
    	  type = ShowerPartnerType::QCDColourLine;
	  Energy startingScale = angularOrdered ? particle.scales().QCD_c : particle.scales().QCD_c_noAO;
    	  newKin= cit->second.first->
    	    generateNextTimeBranching(startingScale,cit->second.second,
				      particle.id()!=cit->first,0.5*enhance,_deTuning);
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
    	  ShoKinPtr newKin2 = cit->second.first->
	    generateNextTimeBranching(startingScale,cit->second.second,
				      particle.id()!=cit->first,0.5*enhance,_deTuning);
    	  // pick the one with the highest scale
    	  if( ( newKin && newKin2 && newKin2->scale() > newKin->scale()) ||
    	      (!newKin && newKin2) ) {
    	    newKin = newKin2;
    	    type = ShowerPartnerType::QCDAntiColourLine;
    	  }
    	}
     	// other g -> q qbar
     	else {
     	  Energy startingScale = angularOrdered ? 
	    max(particle.scales().QCD_c     , particle.scales().QCD_ac     ) : 
	    max(particle.scales().QCD_c_noAO, particle.scales().QCD_ac_noAO);
    	  newKin= cit->second.first->
    	    generateNextTimeBranching(startingScale, cit->second.second,
				      particle.id()!=cit->first,enhance,_deTuning);
    	  type = UseRandom::rndbool() ? 
    	    ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine;
	}
      }
      // everything else q-> qg etc
      else {
	Energy startingScale;
	if(particle.hasColour()) {
	  type = ShowerPartnerType::QCDColourLine;
	  startingScale = angularOrdered ? particle.scales().QCD_c  : particle.scales().QCD_c_noAO;
	}
	else {
	  type = ShowerPartnerType::QCDAntiColourLine;
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	}
	newKin= cit->second.first->
	  generateNextTimeBranching(startingScale,cit->second.second,
				    particle.id()!=cit->first,enhance,_deTuning);
      }
    }
    // shouldn't be anything else
    else
      assert(false);
    // if no kinematics contine
    if(!newKin) continue;
    // select highest scale 
    if( newKin->scale() > newQ ) {
      kinematics  = newKin;
      newQ        = newKin->scale();
      ids         = cit->second.second;
      sudakov     = cit->second.first;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
				    ShowerPartnerType::Undefined);
  // If a branching has been selected initialize it
  kinematics->initialize(particle,PPtr());
  // and return it
  return Branching(kinematics, ids,sudakov,partnerType);
}

Branching SplittingGenerator::
chooseDecayBranching(ShowerParticle &particle,
		     const ShowerParticle::EvolutionScales & stoppingScales,
		     Energy minmass, double enhance,
		     ShowerInteraction::Type interaction) const {
  Energy newQ = Constants::MaxEnergy;
  ShoKinPtr kinematics;
  SudakovPtr sudakov;
  ShowerPartnerType::Type partnerType(ShowerPartnerType::Undefined);
  IdList ids;
  // First, find the eventual branching, corresponding to the lowest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), IdList(),SudakovPtr(),ShowerPartnerType::Undefined);
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit)  {
    // check interaction doesn't change flavour
    if(cit->second.second[1]!=index&&cit->second.second[2]!=index) continue;
    // check either right interaction or doing both
    if(interaction != cit->second.first->interactionType() &&
       interaction != ShowerInteraction::Both ) continue;
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.first->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    ShowerPartnerType::Type type;
    // work out which starting scale we need
    if(cit->second.first->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy stoppingScale = angularOrdered ? stoppingScales.QED    : stoppingScales.QED_noAO;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      if(startingScale < stoppingScale ) { 
    	newKin = cit->second.first->
    	  generateNextDecayBranching(startingScale,stoppingScale,minmass,cit->second.second,
    				     particle.id()!=cit->first,enhance,_deTuning);
      }
    }
    else if(cit->second.first->interactionType()==ShowerInteraction::QCD) {
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.first->splittingFn()->colourStructure()==OctetOctetOctet) {
	  Energy stoppingColour = angularOrdered ? stoppingScales.QCD_c     : stoppingScales.QCD_c_noAO;
	  Energy stoppingAnti   = angularOrdered ? stoppingScales.QCD_ac    : stoppingScales.QCD_ac_noAO;
	  Energy startingColour = angularOrdered ? particle.scales().QCD_c  : particle.scales().QCD_c_noAO;
	  Energy startingAnti   = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	  type = ShowerPartnerType::QCDColourLine;
	  if(startingColour<stoppingColour) {
	    newKin= cit->second.first->	
	      generateNextDecayBranching(startingColour,stoppingColour,minmass,
					 cit->second.second,
					 particle.id()!=cit->first,0.5*enhance,_deTuning);
	  }
	  ShoKinPtr newKin2; 
	  if(startingAnti<stoppingAnti) {
	    newKin2 = cit->second.first->
	      generateNextDecayBranching(startingAnti,stoppingAnti,minmass,
					 cit->second.second,
					 particle.id()!=cit->first,0.5*enhance,_deTuning);
	  }
	  // pick the one with the lowest scale
	  if( (newKin&&newKin2&&newKin2->scale()<newKin->scale()) ||
	      (!newKin&&newKin2) ) {
	    newKin = newKin2;
	    type = ShowerPartnerType::QCDAntiColourLine;
	  }
	}
	// other
	else {
	  assert(false);
	}
      }
      // everything else
      else {
	Energy startingScale,stoppingScale;
	if(particle.hasColour()) {
	  type = ShowerPartnerType::QCDColourLine;
	  stoppingScale = angularOrdered ? stoppingScales.QCD_c     : stoppingScales.QCD_c_noAO;
	  startingScale = angularOrdered ? particle.scales().QCD_c  : particle.scales().QCD_c_noAO;
	}
	else {
	  type = ShowerPartnerType::QCDAntiColourLine;
	  stoppingScale = angularOrdered ? stoppingScales.QCD_ac    : stoppingScales.QCD_ac_noAO;
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	}
	if(startingScale < stoppingScale ) { 
	  newKin = cit->second.first->
	    generateNextDecayBranching(startingScale,stoppingScale,minmass,cit->second.second,
				       particle.id()!=cit->first,enhance,_deTuning);
	}
      }
    }
    // shouldn't be anything else
    else
      assert(false);
    if(!newKin) continue;
    // select highest scale
    if(newKin->scale() < newQ ) {
      newQ = newKin->scale();
      ids = cit->second.second;
      kinematics=newKin;
      sudakov=cit->second.first;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
				    ShowerPartnerType::Undefined);
  // initialize the branching
  kinematics->initialize(particle,PPtr());
  // and generate phi
  kinematics->phi(sudakov->generatePhiDecay(particle,ids,kinematics));
  // and return it
  return Branching(kinematics, ids,sudakov,partnerType);
}

Branching SplittingGenerator::
chooseBackwardBranching(ShowerParticle &particle,PPtr beamparticle,
			double enhance,
			Ptr<BeamParticleData>::transient_const_pointer beam,
			ShowerInteraction::Type type,
			tcPDFPtr pdf, Energy freeze) const {
  Energy newQ=ZERO;
  ShoKinPtr kinematics=ShoKinPtr();
  ShowerPartnerType::Type partnerType(ShowerPartnerType::Undefined);
  SudakovPtr sudakov;
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.id());
  // if no possible branching return
  if(_bbranchings.find(index) == _bbranchings.end())
    return Branching(ShoKinPtr(), IdList(),SudakovPtr(),ShowerPartnerType::Undefined);
  // otherwise select branching
  for(BranchingList::const_iterator cit = _bbranchings.lower_bound(index); 
      cit != _bbranchings.upper_bound(index); ++cit ) {
    // check either right interaction or doing both
    if(type != cit->second.first->interactionType() &&
       type != ShowerInteraction::Both ) continue;
    // setup the PDF
    cit->second.first->setPDF(pdf,freeze);
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.first->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    ShowerPartnerType::Type type;
    if(cit->second.first->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      newKin=cit->second.first->
    	generateNextSpaceBranching(startingScale,cit->second.second, particle.x(),
    				   particle.id()!=cit->first,enhance,beam,_deTuning);
    }
    else if(cit->second.first->interactionType()==ShowerInteraction::QCD) { 
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.first->splittingFn()->colourStructure()==OctetOctetOctet) {
    	  type = ShowerPartnerType::QCDColourLine;
	  Energy startingScale = angularOrdered ? particle.scales().QCD_c : particle.scales().QCD_c_noAO;
	  newKin = cit->second.first->
	    generateNextSpaceBranching(startingScale,cit->second.second, particle.x(),
				       particle.id()!=cit->first,0.5*enhance,beam,_deTuning);
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	  ShoKinPtr newKin2 = cit->second.first->
	    generateNextSpaceBranching(startingScale,cit->second.second, particle.x(),
				       particle.id()!=cit->first,0.5*enhance,beam,_deTuning);
	  // pick the one with the highest scale
	  if( (newKin&&newKin2&&newKin2->scale()>newKin->scale()) ||
	      (!newKin&&newKin2) ) {
	    newKin = newKin2;
	    type = ShowerPartnerType::QCDAntiColourLine;
	  }
	}
	else {
     	  Energy startingScale = angularOrdered ? 
	    max(particle.scales().QCD_c     , particle.scales().QCD_ac    ) : 
	    max(particle.scales().QCD_c_noAO, particle.scales().QCD_ac_noAO);
	  type = UseRandom::rndbool() ? 
	    ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine;
	  newKin=cit->second.first->
	    generateNextSpaceBranching(startingScale,cit->second.second, particle.x(),
				       particle.id()!=cit->first,enhance,beam,_deTuning);
	}
      }
      // everything else
      else {
	Energy startingScale;
	if(particle.hasColour()) {
	  type = ShowerPartnerType::QCDColourLine;
	  startingScale = angularOrdered ? particle.scales().QCD_c  : particle.scales().QCD_c_noAO;
	}
	else {
	  type = ShowerPartnerType::QCDAntiColourLine;
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	}
    	newKin=cit->second.first->
    	  generateNextSpaceBranching(startingScale,cit->second.second, particle.x(),
    				     particle.id()!=cit->first,enhance,beam,_deTuning);
      }
    }
    // shouldn't be anything else
    else
      assert(false);
    // if no kinematics contine
    if(!newKin) continue;
    // select highest scale
    if(newKin->scale() > newQ) {
      newQ = newKin->scale();
      kinematics=newKin;
      ids = cit->second.second;
      sudakov=cit->second.first;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics) return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
				   ShowerPartnerType::Undefined);
  // initialize the ShowerKinematics 
  // and return it
  kinematics->initialize(particle,beamparticle);
  // and generate phi
  kinematics->phi(sudakov->generatePhiBackward(particle,ids,kinematics));
  // return the answer
  return Branching(kinematics, ids,sudakov,partnerType);
}

void SplittingGenerator::rebind(const TranslationMap & trans) {
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit)
    {(cit->second).first=trans.translate((cit->second).first);}
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit)
    {(cit->second).first=trans.translate((cit->second).first);}
  Interfaced::rebind(trans);
}

IVector SplittingGenerator::getReferences() {
  IVector ret = Interfaced::getReferences();
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit)
    {ret.push_back((cit->second).first);}
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit)
    {ret.push_back((cit->second).first);}
  return ret;
}

void SplittingGenerator::factorizationScaleFactor(double f) {
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit)
    {(cit->second).first->factorizationScaleFactor(f);}
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit)
    {(cit->second).first->factorizationScaleFactor(f);}
}

void SplittingGenerator::renormalizationScaleFactor(double f) {
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit)
    {(cit->second).first->renormalizationScaleFactor(f);}
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit)
    {(cit->second).first->renormalizationScaleFactor(f);}
}

