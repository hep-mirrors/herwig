// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingGenerator class.
//

#include "SplittingGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SplittingGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/Repository.h"
#include "Herwig++/Shower/Kinematics/FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Kinematics/Decay_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Kinematics/IS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Kinematics/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _qcdinteractionMode << _qedinteractionMode << _ewkinteractionMode
     << _isr_Mode << _isr_qcdMode << _isr_qedMode << _isr_ewkMode     
     << _fsr_Mode << _fsr_qcdMode << _fsr_qedMode << _fsr_ewkMode
     << _showerVariables << _bbranchings << _fbranchings;
}

void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_qcdinteractionMode >> _qedinteractionMode >> _ewkinteractionMode
     >>	_isr_Mode >> _isr_qcdMode >> _isr_qedMode >> _isr_ewkMode 
     >> _fsr_Mode >> _fsr_qcdMode >> _fsr_qedMode >> _fsr_ewkMode
     >> _showerVariables >> _bbranchings >> _fbranchings;
}

ClassDescription<SplittingGenerator> SplittingGenerator::initSplittingGenerator;
// Definition of the static class description member.

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

  static Switch<SplittingGenerator, bool> interfaceQCDinteractionMode
    ("OnOffQCDinteractionMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_qcdinteractionMode, 1, false, false);
  static SwitchOption interfaceQCDinteractionMode0
    (interfaceQCDinteractionMode,"QCDinteraction-OFF","QCD interaction is OFF", 0);
  static SwitchOption interfaceQCDinteractionMode1
    (interfaceQCDinteractionMode,"QCDinteraction-ON","QCD interaction is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceQEDinteractionMode
    ("OnOffQEDinteractionMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_qedinteractionMode, 0, false, false);
  static SwitchOption interfaceQEDinteractionMode0
    (interfaceQEDinteractionMode,"QEDinteraction-OFF","QED interaction is OFF", 0);
  static SwitchOption interfaceQEDinteractionMode1
    (interfaceQEDinteractionMode,"QEDinteraction-ON","QED interaction is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceEWKinteractionMode
    ("OnOffEWKinteractionMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_ewkinteractionMode, 0, false, false);
  static SwitchOption interfaceEWKinteractionMode0
    (interfaceEWKinteractionMode,"EWKinteraction-OFF","EWK interaction is OFF", 0);
  static SwitchOption interfaceEWKinteractionMode1
    (interfaceEWKinteractionMode,"EWKinteraction-ON","EWK interaction is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISRMode
    ("OnOffISRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_isr_Mode, 1, false, false);
  static SwitchOption interfaceISRMode0
    (interfaceISRMode,"ISR-OFF","ISR (Initial State Radiation) is OFF", 0);
  static SwitchOption interfaceISRMode1
    (interfaceISRMode,"ISR-ON","ISR (Initial State Radiation) is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_qcdMode
    ("OnOffISR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_isr_qcdMode, 1, false, false);
  static SwitchOption interfaceISR_qcdMode0
    (interfaceISR_qcdMode,"ISR_QCD-OFF","QCD ISR is OFF", 0);
  static SwitchOption interfaceISR_qcdMode1
    (interfaceISR_qcdMode,"ISR_QCD-ON","QCD ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_qedMode
    ("OnOffISR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_isr_qedMode, 1, false, false);
  static SwitchOption interfaceISR_qedMode0
    (interfaceISR_qedMode,"ISR_QED-OFF","QED ISR is OFF", 0);
  static SwitchOption interfaceISR_qedMode1
    (interfaceISR_qedMode,"ISR_QED-ON","QED ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceISR_ewkMode
    ("OnOffISR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_isr_ewkMode, 1, false, false);
  static SwitchOption interfaceISR_ewkMode0
    (interfaceISR_ewkMode,"ISR_EWK-OFF","EWK ISR is OFF", 0);
  static SwitchOption interfaceISR_ewkMode1
    (interfaceISR_ewkMode,"ISR_EWK-ON","EWK ISR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSRMode
    ("OnOffFSRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_fsr_Mode, 1, false, false);
  static SwitchOption interfaceFSRMode0
    (interfaceFSRMode,"FSR-OFF","FSR (Final State Radiation) is OFF", 0);
  static SwitchOption interfaceFSRMode1
    (interfaceFSRMode,"FSR-ON","FSR (Final State Radiation) is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_qcdMode
    ("OnOffFSR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_fsr_qcdMode, 1, false, false);
  static SwitchOption interfaceFSR_qcdMode0
    (interfaceFSR_qcdMode,"FSR_QCD-OFF","QCD FSR is OFF", 0);
  static SwitchOption interfaceFSR_qcdMode1
    (interfaceFSR_qcdMode,"FSR_QCD-ON","QCD FSR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_qedMode
    ("OnOffFSR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_fsr_qedMode, 1, false, false);
  static SwitchOption interfaceFSR_qedMode0
    (interfaceFSR_qedMode,"FSR_QED-OFF","QED FSR is OFF", 0);
  static SwitchOption interfaceFSR_qedMode1
    (interfaceFSR_qedMode,"FSR_QED-ON","QED FSR is ON", 1);

  static Switch<SplittingGenerator, bool> interfaceFSR_ewkMode
    ("OnOffFSR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_fsr_ewkMode, 1, false, false);

  static SwitchOption interfaceFSR_ewkMode0
    (interfaceFSR_ewkMode,"FSR_EWK-OFF","EWK FSR is OFF", 0);

  static SwitchOption interfaceFSR_ewkMode1
    (interfaceFSR_ewkMode,"FSR_EWK-ON","EWK FSR is ON", 1);

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

Branching SplittingGenerator::chooseForwardBranching(ShowerParticle &particle) const {
  Timer<1200> timer("SplittingGenerator::chooseForwardBranching");
  Energy newQ = Energy(), newpT=Energy();
  double newZ(0.0),newPhi(0.0);
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), SudakovPtr(), IdList());
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) 
    {
      tSudakovPtr candidateSudakov = cit->second.first;
      assert(candidateSudakov);
      ShowerIndex::InteractionType i = candidateSudakov->interactionType();
      // check size of scales beforehand...
      Energy candidateNewQ = Energy();
      if(particle.evolutionScales()[i] > _showerVariables->cutoffQScale(i)) 
	candidateNewQ = candidateSudakov->
	  generateNextTimeBranching(particle.evolutionScales()[i], 
				    cit->second.second,particle.id()!=cit->first);
      // select highest scale 
      if(candidateNewQ > newQ && candidateNewQ <= particle.evolutionScales()[i]) 
	{
	  newQ = candidateNewQ;
	  sudakov = candidateSudakov;
	  ids = cit->second.second;
	  newZ = sudakov->z();
	  newPhi = sudakov->phi();
	  newpT = sudakov->pT();
	}
    }
  // return empty branching if nothing happened
  if(!sudakov)  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());
  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  // For the time being we are considering only 1->2 branching
  tSplittingFnPtr splitFun = sudakov->splittingFn();
  assert(splitFun);
  // get the reference vectors
  Lorentz5Momentum p, n;
  finalStateBasisVectors(particle,sudakov->splittingFn()->interactionType(),
			 p,n);
  // construct the shower kinematics object
  FS_QtildaShowerKinematics1to2Ptr showerKin = 
    new_ptr(FS_QtildaShowerKinematics1to2(p, n,showerVariables()));
  showerKin->qtilde(newQ);
  showerKin->setResScale(sudakov->resScale());
  showerKin->setKinScale(sudakov->kinScale()); 
  showerKin->z(newZ);
  showerKin->phi(newPhi);
  showerKin->pT(newpT);
  return Branching(showerKin, sudakov, ids);
}


Branching SplittingGenerator::chooseDecayBranching(ShowerParticle &particle,
						   vector<Energy> stoppingScale,
						   Energy minmass) const 
{
  Energy newQ = ShowerVariables::HUGEMASS,newpT(0.);
  double newZ(0.0),newPhi(0.0);
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the lowest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), SudakovPtr(), IdList());
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) 
    {
      tSudakovPtr candidateSudakov = cit->second.first;
      assert(candidateSudakov);
      ShowerIndex::InteractionType i = candidateSudakov->interactionType();
      Energy candidateNewQ = ShowerVariables::HUGEMASS;
      if(particle.evolutionScales()[i] < stoppingScale[i]) 
  	candidateNewQ = candidateSudakov->
 	  generateNextDecayBranching(particle.evolutionScales()[i],
				     stoppingScale[i],minmass,
 				     cit->second.second,
				     particle.id()!=cit->first);
      if(candidateNewQ < newQ && candidateNewQ > particle.evolutionScales()[i])
	{
	  newQ = candidateNewQ;
	  sudakov = candidateSudakov;
	  ids = cit->second.second;
	  newZ = sudakov->z();
	  newPhi = sudakov->phi();
	  newpT = sudakov->pT();
	}
    }
  // return empty branching if nothing happened
  if(!sudakov)  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());
  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  // For the time being we are considering only 1->2 branching
  tSplittingFnPtr splitFun = sudakov->splittingFn();
  assert(splitFun);
  Lorentz5Momentum p, n, ppartner, pcm;
  assert(particle.perturbative()!=1);
  // this is for the initial decaying particle
  if(particle.perturbative()==2)
     {
       p = particle.momentum();
       ShowerParticlePtr partner=particle.partners()[sudakov->splittingFn()->
						     interactionType()];
       Lorentz5Momentum ppartner(partner->momentum());
       if(partner->getThePEGBase()) ppartner=partner->getThePEGBase()->momentum();
       pcm=ppartner;
       Hep3Vector boost(p.findBoostToCM());
       pcm.boost(boost);
       n = Lorentz5Momentum( 0.0,0.5*p.mass()*pcm.vect().unit()); 
       n.boost( -boost);
     }
   else 
     {
       tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	 ->showerKinematics();
       p = kin->getBasis()[0];
       n = kin->getBasis()[1];
     }
   Decay_QtildaShowerKinematics1to2Ptr showerKin = 
     new_ptr(Decay_QtildaShowerKinematics1to2(p, n,showerVariables()));
   showerKin->qtilde(newQ);
   showerKin->setResScale(sudakov->resScale());
   showerKin->setKinScale(sudakov->kinScale()); 
   showerKin->z(newZ);
   showerKin->phi(newPhi);
   showerKin->pT(newpT);
   return Branching(showerKin, sudakov, ids);
}

Branching SplittingGenerator::chooseBackwardBranching(ShowerParticle &particle) const {
  Energy newQ=Energy(),newpT=Energy();
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;
  double newZ(0.0), newPhi(0.0);  
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.id());
  // if no possible branching return
  if(_bbranchings.find(index) == _bbranchings.end()) 
    return Branching(ShoKinPtr(), tSudakovPtr(), IdList());
  // select the branching
  for(BranchingList::const_iterator cit = _bbranchings.lower_bound(index); 
      cit != _bbranchings.upper_bound(index); ++cit ) 
    {
      tSudakovPtr candidateSudakov = cit->second.first;
      assert(candidateSudakov);
      Energy candidateNewQ = Energy();
      ShowerIndex::InteractionType i = candidateSudakov->interactionType(); 
      candidateNewQ = candidateSudakov->
	generateNextSpaceBranching(particle.evolutionScales()[i],
				   cit->second.second, particle.x(),
				   particle.id()!=cit->first);
      if(candidateNewQ > newQ) 
	{
	  newQ = candidateNewQ;
	  sudakov = candidateSudakov;
	  ids = cit->second.second;
	  newZ = sudakov->z();
	  newPhi = sudakov->phi();
	  newpT  = sudakov->pT();
	}
    } 
  // return empty branching if nothing happened
  if(!sudakov) return Branching(ShoKinPtr(), tSudakovPtr(), IdList());
  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  tSplittingFnPtr splitFun = sudakov->splittingFn();
  assert(splitFun);
  // For the time being we are considering only 1->2 branching
  Lorentz5Momentum p, n, pthis, ppartner, pcm;
  assert(particle.perturbative()!=2);
  if(particle.perturbative()==1) 
    {
      pcm = particle.parents()[0]->momentum();
      p = Lorentz5Momentum(0.0, pcm.vect());
      n = Lorentz5Momentum(0.0, -pcm.vect());
    } 
  else 
    {
      p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
	->showerKinematics()->getBasis()[0];
      n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
	->showerKinematics()->getBasis()[1];
    }
  // create the shower kinematics
  Ptr<IS_QtildaShowerKinematics1to2>::pointer 
    showerKin = new_ptr(IS_QtildaShowerKinematics1to2(p, n,showerVariables()));
  showerKin->qtilde(newQ);
  showerKin->setResScale(sudakov->resScale());
  showerKin->setKinScale(sudakov->kinScale()); 
  showerKin->z(newZ);
  showerKin->phi(newPhi);
  showerKin->pT(newpT);
  // return the answer
  return Branching(showerKin, sudakov, ids);
}

void SplittingGenerator::finalStateBasisVectors(ShowerParticle particle,
						ShowerIndex::InteractionType type,
						Lorentz5Momentum & p,
						Lorentz5Momentum & n) const
{
  if(particle.perturbative()!=0)
    {
      // find the partner and its momentum
      ShowerParticlePtr partner=particle.partners()[type];
      Lorentz5Momentum ppartner(partner->momentum());
      if(partner->getThePEGBase()) ppartner=partner->getThePEGBase()->momentum();
      // momentum of the emitting particle
      p = particle.momentum();
      Lorentz5Momentum pcm;
      // if the partner is a final-state particle then the reference
      // vector is along the partner in the rest frame of the pair
      if(partner->isFinalState())
	{
	  Hep3Vector boost=(p + ppartner).findBoostToCM();
	  pcm = ppartner;
	  pcm.boost(boost);
	  n = Lorentz5Momentum(0.,pcm.vect());
	  n.boost( -boost);
	}
      else if(!partner->isFinalState())
	{
	  // if the partner is an initial-state particle then the reference
	  // vector is along the partner which should be massless
	  if(particle.perturbative()==1)
	    {n = Lorentz5Momentum(0.,ppartner.vect());}
	  // if the partner is an initial-state decaying particle then the reference
	  // vector is along the backwards direction in rest frame of decaying particle
	  else
	    {
	      Hep3Vector boost=ppartner.findBoostToCM();
	      pcm = p;
	      pcm.boost(boost);
	      n = Lorentz5Momentum( 0.0, -pcm.vect()); 
	      n.boost( -boost);
	    } 
	} 
    }
  else if(particle.initiatesTLS())
    {
      tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>
	(particle.parents()[0]->children()[0])->showerKinematics();
      p = kin->getBasis()[0];
      n = kin->getBasis()[1];
    }
  else 
    {
      tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	->showerKinematics();
      p = kin->getBasis()[0];
      n = kin->getBasis()[1];
    }
}
