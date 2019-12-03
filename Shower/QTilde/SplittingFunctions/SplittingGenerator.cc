// -*- C++ -*-
//
// SplittingGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "ThePEG/Utilities/Rebinder.h"
#include <cassert>
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

namespace {

bool checkInteraction(ShowerInteraction allowed,
		      ShowerInteraction splitting) {
  if(allowed==ShowerInteraction::Both)
    return true;
  else if(allowed == splitting)
    return true;
  else
    return false;
}

}

DescribeClass<SplittingGenerator,Interfaced>
describeSplittingGenerator ("Herwig::SplittingGenerator","");

IBPtr SplittingGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr SplittingGenerator::fullclone() const {
  return new_ptr(*this);
}


void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _bbranchings << _fbranchings << _deTuning;
}

void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _bbranchings >> _fbranchings >> _deTuning;
}

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

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
  ids.push_back(parent);
  for(vector<tPDPtr>::iterator it = products.begin(); it!=products.end(); ++it)
    ids.push_back(*it);
  // check splitting can handle this
  if(!s->splittingFn()->accept(ids)) 
    return "Error: Sudakov " + sudakov + " SplittingFunction can't handle particles\n";
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
  ids.push_back(parent);
  for(vector<tPDPtr>::iterator it = products.begin(); it!=products.end(); ++it)
    ids.push_back(*it);
  // check splitting can handle this
  if(!s->splittingFn()->accept(ids)) 
    return "Error: Sudakov " + sudakov + " Splitting Function can't handle particles\n";
  // delete from map
  deleteFromMap(ids,s,final);
  return "";
}

void SplittingGenerator::addToMap(const IdList &ids, const SudakovPtr &s, bool final) {
  if(!final) {
      // search if the branching was already included.
    auto binsert =BranchingInsert(abs(ids[1]->id()),BranchingElement(s,ids));
      // get the range of already inserted splittings.
    auto eqrange=_bbranchings.equal_range(binsert.first);
    for(auto it = eqrange.first; it != eqrange.second; ++it){
      if((*it).second == binsert.second)
        throw Exception()<<"SplittingGenerator: Trying to insert existing splitting.\n"
        << Exception::setuperror;
    }
    _bbranchings.insert(binsert);
    s->addSplitting(ids);
  }
  else {
      // search if the branching was already included.
    auto binsert =BranchingInsert(abs(ids[0]->id()),BranchingElement(s,ids));
      // get the range of already inserted splittings.
    auto eqrange=_fbranchings.equal_range(binsert.first);
    for(auto it = eqrange.first; it != eqrange.second; ++it){
      if((*it).second ==binsert.second)
        throw Exception()<<"SplittingGenerator: Trying to insert existing splitting.\n"
        << Exception::setuperror;
    }
    
    _fbranchings.insert(binsert);
    s->addSplitting(ids);
  }
}

void SplittingGenerator::deleteFromMap(const IdList &ids, 
				       const SudakovPtr &s, bool final) {
  bool didRemove=false;
  if(!final) {
    pair<BranchingList::iterator,BranchingList::iterator> 
      range = _bbranchings.equal_range(abs(ids[1]->id()));
    for(BranchingList::iterator it=range.first;
	it!=range.second&&it!=_bbranchings.end();++it) {
      if(it->second.sudakov==s&&it->second.particles==ids) {
	BranchingList::iterator it2=it;
	--it;
	_bbranchings.erase(it2);
    didRemove=true;
      }
    }
    s->removeSplitting(ids);
  }
  else {
    pair<BranchingList::iterator,BranchingList::iterator> 
      range = _fbranchings.equal_range(abs(ids[0]->id()));
    for(BranchingList::iterator it=range.first;
	it!=range.second&&it!=_fbranchings.end();++it) {
      if(it->second.sudakov==s&&it->second.particles==ids) {
	BranchingList::iterator it2 = it;
	--it;
	_fbranchings.erase(it2);
    didRemove=true;
      }
    }
    s->removeSplitting(ids);
  }
  if (!didRemove)
    throw Exception()<<"SplittingGenerator: Try to remove non existing splitting.\n"
                      << Exception::setuperror;
}

Branching SplittingGenerator::chooseForwardBranching(ShowerParticle &particle,
						     double enhance,
						     ShowerInteraction type) const {
  RhoDMatrix rho;
  bool rhoCalc(false);
  Energy newQ = ZERO;
  ShoKinPtr kinematics = ShoKinPtr();
  ShowerPartnerType partnerType(ShowerPartnerType::Undefined);
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
    if(!checkInteraction(type,cit->second.sudakov->interactionType())) continue;
    if(!rhoCalc) {
      rho = particle.extractRhoMatrix(true);
      rhoCalc = true;
    }
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.sudakov->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    ShowerPartnerType type;
    IdList particles = particle.id()!=cit->first ? cit->second.conjugateParticles : cit->second.particles;
    // work out which starting scale we need
    if(cit->second.sudakov->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      newKin = cit->second.sudakov->
    	generateNextTimeBranching(startingScale,particles,rho,enhance,
				  _deTuning);
    }
    else if(cit->second.sudakov->interactionType()==ShowerInteraction::QCD) {
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.sudakov->splittingFn()->colourStructure()==OctetOctetOctet) {
    	  type = ShowerPartnerType::QCDColourLine;
	  Energy startingScale = angularOrdered ? particle.scales().QCD_c : particle.scales().QCD_c_noAO;
    	  newKin= cit->second.sudakov->
    	    generateNextTimeBranching(startingScale,particles,rho,0.5*enhance,
				      _deTuning);
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
    	  ShoKinPtr newKin2 = cit->second.sudakov->
	    generateNextTimeBranching(startingScale,particles,rho,0.5*enhance,
				      _deTuning);
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
    	  newKin= cit->second.sudakov->
    	    generateNextTimeBranching(startingScale,particles,rho,enhance,
				      _deTuning);
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
	newKin= cit->second.sudakov->
	  generateNextTimeBranching(startingScale,particles,rho,enhance,
				    _deTuning);
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
      ids         = particles;
      sudakov     = cit->second.sudakov;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics) {
    if ( particle.spinInfo() ) particle.spinInfo()->undecay();
    return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
		     ShowerPartnerType::Undefined);
  }
  // if not hard generate phi
  kinematics->phi(sudakov->generatePhiForward(particle,ids,kinematics,rho));
  // and return it
  return Branching(kinematics, ids,sudakov,partnerType);
}

Branching SplittingGenerator::
chooseDecayBranching(ShowerParticle &particle,
		     const ShowerParticle::EvolutionScales & stoppingScales,
		     Energy minmass, double enhance,
		     ShowerInteraction interaction) const {
  RhoDMatrix rho(particle.dataPtr()->iSpin());
  Energy newQ = Constants::MaxEnergy;
  ShoKinPtr kinematics;
  SudakovPtr sudakov;
  ShowerPartnerType partnerType(ShowerPartnerType::Undefined);
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
    if(cit->second.particles[1]->id()!=index&&cit->second.particles[2]->id()!=index) continue;
    // check either right interaction or doing both
    if(!checkInteraction(interaction,cit->second.sudakov->interactionType())) continue;
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.sudakov->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    IdList particles = particle.id()!=cit->first ? cit->second.conjugateParticles : cit->second.particles;
    ShowerPartnerType type;
    // work out which starting scale we need
    if(cit->second.sudakov->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy stoppingScale = angularOrdered ? stoppingScales.QED    : stoppingScales.QED_noAO;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      if(startingScale < stoppingScale ) { 
    	newKin = cit->second.sudakov->
    	  generateNextDecayBranching(startingScale,stoppingScale,minmass,particles,rho,enhance,_deTuning);
      }
    }
    else if(cit->second.sudakov->interactionType()==ShowerInteraction::QCD) {
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.sudakov->splittingFn()->colourStructure()==OctetOctetOctet) {
	  Energy stoppingColour = angularOrdered ? stoppingScales.QCD_c     : stoppingScales.QCD_c_noAO;
	  Energy stoppingAnti   = angularOrdered ? stoppingScales.QCD_ac    : stoppingScales.QCD_ac_noAO;
	  Energy startingColour = angularOrdered ? particle.scales().QCD_c  : particle.scales().QCD_c_noAO;
	  Energy startingAnti   = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	  type = ShowerPartnerType::QCDColourLine;
	  if(startingColour<stoppingColour) {
	    newKin= cit->second.sudakov->	
	      generateNextDecayBranching(startingColour,stoppingColour,minmass,
					 particles,rho,0.5*enhance,_deTuning);
	  }
	  ShoKinPtr newKin2; 
	  if(startingAnti<stoppingAnti) {
	    newKin2 = cit->second.sudakov->
	      generateNextDecayBranching(startingAnti,stoppingAnti,minmass,particles,rho,0.5*enhance,_deTuning);
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
	  newKin = cit->second.sudakov->
	    generateNextDecayBranching(startingScale,stoppingScale,minmass,particles,rho,enhance,_deTuning);
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
      ids = particles;
      kinematics=newKin;
      sudakov=cit->second.sudakov;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
				    ShowerPartnerType::Undefined);
  // and generate phi
  kinematics->phi(sudakov->generatePhiDecay(particle,ids,kinematics,rho));
  // and return it
  return Branching(kinematics, ids,sudakov,partnerType);
}

Branching SplittingGenerator::
chooseBackwardBranching(ShowerParticle &particle,PPtr ,
			double enhance,
			Ptr<BeamParticleData>::transient_const_pointer beam,
			ShowerInteraction type,
			tcPDFPtr pdf, Energy freeze) const {
  RhoDMatrix rho;
  bool rhoCalc(false);
  Energy newQ=ZERO;
  ShoKinPtr kinematics=ShoKinPtr();
  ShowerPartnerType partnerType(ShowerPartnerType::Undefined);
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
    if(!checkInteraction(type,cit->second.sudakov->interactionType())) continue;
    // setup the PDF
    cit->second.sudakov->setPDF(pdf,freeze);
    //calc rho as needed
    if(!rhoCalc) {
      rho = particle.extractRhoMatrix(false);
      rhoCalc = true;
    }
    // whether or not this interaction should be angular ordered
    bool angularOrdered = cit->second.sudakov->splittingFn()->angularOrdered();
    ShoKinPtr newKin;
    IdList particles = particle.id()!=cit->first ? cit->second.conjugateParticles : cit->second.particles;
    ShowerPartnerType type;
    if(cit->second.sudakov->interactionType()==ShowerInteraction::QED) {
      type = ShowerPartnerType::QED;
      Energy startingScale = angularOrdered ? particle.scales().QED : particle.scales().QED_noAO;
      newKin=cit->second.sudakov->
    	generateNextSpaceBranching(startingScale,particles,particle.x(),rho,enhance,beam,_deTuning);
    }
    else if(cit->second.sudakov->interactionType()==ShowerInteraction::QCD) { 
      // special for octets
      if(particle.dataPtr()->iColour()==PDT::Colour8) {
	// octet -> octet octet
	if(cit->second.sudakov->splittingFn()->colourStructure()==OctetOctetOctet) {
    	  type = ShowerPartnerType::QCDColourLine;
	  Energy startingScale = angularOrdered ? particle.scales().QCD_c : particle.scales().QCD_c_noAO;
	  newKin = cit->second.sudakov->
	    generateNextSpaceBranching(startingScale,particles, particle.x(),rho,0.5*enhance,beam,_deTuning);
	  startingScale = angularOrdered ? particle.scales().QCD_ac : particle.scales().QCD_ac_noAO;
	  ShoKinPtr newKin2 = cit->second.sudakov->
	    generateNextSpaceBranching(startingScale,particles, particle.x(),rho,0.5*enhance,beam,_deTuning);
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
	  newKin=cit->second.sudakov->
	    generateNextSpaceBranching(startingScale,particles, particle.x(),rho,enhance,beam,_deTuning);
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
    	newKin=cit->second.sudakov->
    	  generateNextSpaceBranching(startingScale,particles, particle.x(),rho,enhance,beam,_deTuning);
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
      ids = particles;
      sudakov=cit->second.sudakov;
      partnerType = type;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics) {
    if ( particle.spinInfo() ) particle.spinInfo()->undecay();
    return Branching(ShoKinPtr(), IdList(),SudakovPtr(),
		     ShowerPartnerType::Undefined);
  }
  // initialize the ShowerKinematics 
  // and generate phi
  kinematics->phi(sudakov->generatePhiBackward(particle,ids,kinematics,rho));
  // return the answer
  return Branching(kinematics, ids,sudakov,partnerType);
}

void SplittingGenerator::rebind(const TranslationMap & trans) {
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit) {
    (cit->second).sudakov=trans.translate((cit->second).sudakov);
    for(unsigned int ix=0;ix<(cit->second).particles.size();++ix) {
      (cit->second).particles[ix]=trans.translate((cit->second).particles[ix]);
    }
    for(unsigned int ix=0;ix<(cit->second).conjugateParticles.size();++ix) {
      (cit->second).conjugateParticles[ix]=trans.translate((cit->second).conjugateParticles[ix]);
    }
  }
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit) {
    (cit->second).sudakov=trans.translate((cit->second).sudakov);
    for(unsigned int ix=0;ix<(cit->second).particles.size();++ix) {
      (cit->second).particles[ix]=trans.translate((cit->second).particles[ix]);
    }
    for(unsigned int ix=0;ix<(cit->second).conjugateParticles.size();++ix) {
      (cit->second).conjugateParticles[ix]=trans.translate((cit->second).conjugateParticles[ix]);
    }
  }
  Interfaced::rebind(trans);
}

IVector SplittingGenerator::getReferences() {
  IVector ret = Interfaced::getReferences();
  BranchingList::iterator cit;
  for(cit=_fbranchings.begin();cit!=_fbranchings.end();++cit) {
    ret.push_back((cit->second).sudakov);
    for(unsigned int ix=0;ix<(cit->second).particles.size();++ix) 
      ret.push_back(const_ptr_cast<tPDPtr>((cit->second).particles[ix]));
    for(unsigned int ix=0;ix<(cit->second).conjugateParticles.size();++ix) 
      ret.push_back(const_ptr_cast<tPDPtr>((cit->second).conjugateParticles[ix]));
  }
  for(cit=_bbranchings.begin();cit!=_bbranchings.end();++cit) {
    ret.push_back((cit->second).sudakov);
    for(unsigned int ix=0;ix<(cit->second).particles.size();++ix) 
      ret.push_back(const_ptr_cast<tPDPtr>((cit->second).particles[ix]));
    for(unsigned int ix=0;ix<(cit->second).conjugateParticles.size();++ix) 
      ret.push_back(const_ptr_cast<tPDPtr>((cit->second).conjugateParticles[ix]));
  }
  return ret;
}


