// -*- C++ -*-
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
#include "ThePEG/PDT/EnumParticles.h"
#include "ShowerParticle.h"
#include "ThePEG/Handlers/PartialCollisionHandler.h"
#include "ThePEG/Repository/FullEventGenerator.h"
#include "ShowerAlpha.h"
#include "SplittingFunction.h"
#include "SudakovFormFactor.h"
#include "IS_QtildaShowerKinematics1to2.h"
#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConfig.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Repository/Repository.h"

using namespace Herwig;


SplittingGenerator::~SplittingGenerator() {}


void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _QCDinteractionMode << _QEDinteractionMode << _EWKinteractionMode
     << _ISR_Mode << _ISR_QCDMode << _ISR_QEDMode << _ISR_EWKMode
     << _FSR_Mode << _FSR_QCDMode << _FSR_QEDMode << _FSR_EWKMode
     << _showerAlphaQCD << _showerAlphaQED << _showerVariables << _branchings;
}


void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_QCDinteractionMode >> _QEDinteractionMode >> _EWKinteractionMode
     >>	_ISR_Mode >> _ISR_QCDMode >> _ISR_QEDMode >> _ISR_EWKMode 
     >> _FSR_Mode >> _FSR_QCDMode >> _FSR_QEDMode >> _FSR_EWKMode
     >> _showerAlphaQCD >> _showerAlphaQED >> _showerVariables >> _branchings;
  setSVtoAlpha(_showerVariables);
}


ClassDescription<SplittingGenerator> SplittingGenerator::initSplittingGenerator;
// Definition of the static class description member.

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

  static Switch<SplittingGenerator, int> interfaceQCDinteractionMode
    ("OnOffQCDinteractionMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_QCDinteractionMode, 1, false, false);
  static SwitchOption interfaceQCDinteractionMode0
    (interfaceQCDinteractionMode,"QCDinteraction-OFF","QCD interaction is OFF", 0);
  static SwitchOption interfaceQCDinteractionMode1
    (interfaceQCDinteractionMode,"QCDinteraction-ON","QCD interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceQEDinteractionMode
    ("OnOffQEDinteractionMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_QEDinteractionMode, 0, false, false);
  static SwitchOption interfaceQEDinteractionMode0
    (interfaceQEDinteractionMode,"QEDinteraction-OFF","QED interaction is OFF", 0);
  static SwitchOption interfaceQEDinteractionMode1
    (interfaceQEDinteractionMode,"QEDinteraction-ON","QED interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceEWKinteractionMode
    ("OnOffEWKinteractionMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_EWKinteractionMode, 0, false, false);
  static SwitchOption interfaceEWKinteractionMode0
    (interfaceEWKinteractionMode,"EWKinteraction-OFF","EWK interaction is OFF", 0);
  static SwitchOption interfaceEWKinteractionMode1
    (interfaceEWKinteractionMode,"EWKinteraction-ON","EWK interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISRMode
    ("OnOffISRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_ISR_Mode, 0, false, false);
  static SwitchOption interfaceISRMode0
    (interfaceISRMode,"ISR-OFF","ISR (Initial State Radiation) is OFF", 0);
  static SwitchOption interfaceISRMode1
    (interfaceISRMode,"ISR-ON","ISR (Initial State Radiation) is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_QCDMode
    ("OnOffISR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_ISR_QCDMode, 1, false, false);
  static SwitchOption interfaceISR_QCDMode0
    (interfaceISR_QCDMode,"ISR_QCD-OFF","QCD ISR is OFF", 0);
  static SwitchOption interfaceISR_QCDMode1
    (interfaceISR_QCDMode,"ISR_QCD-ON","QCD ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_QEDMode
    ("OnOffISR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_ISR_QEDMode, 1, false, false);
  static SwitchOption interfaceISR_QEDMode0
    (interfaceISR_QEDMode,"ISR_QED-OFF","QED ISR is OFF", 0);
  static SwitchOption interfaceISR_QEDMode1
    (interfaceISR_QEDMode,"ISR_QED-ON","QED ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_EWKMode
    ("OnOffISR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_ISR_EWKMode, 1, false, false);
  static SwitchOption interfaceISR_EWKMode0
    (interfaceISR_EWKMode,"ISR_EWK-OFF","EWK ISR is OFF", 0);
  static SwitchOption interfaceISR_EWKMode1
    (interfaceISR_EWKMode,"ISR_EWK-ON","EWK ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSRMode
    ("OnOffFSRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_FSR_Mode, 1, false, false);
  static SwitchOption interfaceFSRMode0
    (interfaceFSRMode,"FSR-OFF","FSR (Final State Radiation) is OFF", 0);
  static SwitchOption interfaceFSRMode1
    (interfaceFSRMode,"FSR-ON","FSR (Final State Radiation) is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_QCDMode
    ("OnOffFSR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_FSR_QCDMode, 1, false, false);
  static SwitchOption interfaceFSR_QCDMode0
    (interfaceFSR_QCDMode,"FSR_QCD-OFF","QCD FSR is OFF", 0);
  static SwitchOption interfaceFSR_QCDMode1
    (interfaceFSR_QCDMode,"FSR_QCD-ON","QCD FSR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_QEDMode
    ("OnOffFSR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_FSR_QEDMode, 1, false, false);
  static SwitchOption interfaceFSR_QEDMode0
    (interfaceFSR_QEDMode,"FSR_QED-OFF","QED FSR is OFF", 0);
  static SwitchOption interfaceFSR_QEDMode1
    (interfaceFSR_QEDMode,"FSR_QED-ON","QED FSR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_EWKMode
    ("OnOffFSR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_FSR_EWKMode, 1, false, false);

  static SwitchOption interfaceFSR_EWKMode0
    (interfaceFSR_EWKMode,"FSR_EWK-OFF","EWK FSR is OFF", 0);

  static SwitchOption interfaceFSR_EWKMode1
    (interfaceFSR_EWKMode,"FSR_EWK-ON","EWK FSR is ON", 1);

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQED
    ("ShowerAlphaQED", "A reference to the IS_ShowerAlphaQED object", 
     &Herwig::SplittingGenerator::_showerAlphaQED, false, false, true, false, 
     &Herwig::SplittingGenerator::setQED, 0, 0);
  // up to here.

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQCD
    ("ShowerAlphaQCD", "A reference to the IS_ShowerAlphaQCD object", 
     &Herwig::SplittingGenerator::_showerAlphaQCD, false, false, true, false,
     &Herwig::SplittingGenerator::setQCD, 0, 0);

  static Reference<SplittingGenerator,ShowerVariables> interfaceShowerVariables
    ("ShowerVariables", "A reference to the ShowerVariables object", 
     &Herwig::SplittingGenerator::_showerVariables, false, false, true, false, 
     &Herwig::SplittingGenerator::setSVtoAlpha, 0, 0);

  static Command<SplittingGenerator> interfaceAddSplitting
    ("AddSplitting",
     "Adds another splitting to the list of splittings considered "
     "in the shower. Command is a->b,c; SplittingFn, where Sudakov "
     "is the name of a created Sudakov object.",
     &SplittingGenerator::addSplitting);
}

string SplittingGenerator::addSplitting(string arg) {
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
  if(!s) cerr << "Error: Could not load Sudakov " << sudakov << endl;
  IdList ids;
  ids.push_back(parent->id());
  for(vector<tPDPtr>::iterator it = products.begin(); it!=products.end(); ++it)
    ids.push_back((*it)->id());
  long i = parent->id();
  addToMap(i,ids,s);
  return "";
}

//--------------------------------------------------------------------------------

void SplittingGenerator::setSVtoAlpha(ShowerVarsPtr p) {
  _showerVariables = p; 
  if (_showerAlphaQCD) _showerAlphaQCD->setSV(p);
  if (_showerAlphaQED) _showerAlphaQED->setSV(p);
}

void SplittingGenerator::setQED(ShowerAlphaPtr p) {
  _showerAlphaQED = p;
  if(_showerVariables) p->setSV(_showerVariables);
}

void SplittingGenerator::setQCD(ShowerAlphaPtr p) {
  _showerAlphaQCD = p;
  if(_showerVariables) p->setSV(_showerVariables);
}

bool SplittingGenerator::isInteractionON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  switch ( interaction ) {
  case ShowerIndex::QCD : mode = _QCDinteractionMode; break; 
  case ShowerIndex::QED : mode = _QEDinteractionMode; break; 
  case ShowerIndex::EWK : mode = _EWKinteractionMode; break; 
  default: break;
  }
  return mode;
}

bool SplittingGenerator::isISRadiationON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  if ( isInteractionON(interaction) && isISRadiationON() ) { 
    switch ( interaction ) {
    case ShowerIndex::QCD : mode = _ISR_QCDMode; break; 
    case ShowerIndex::QED : mode = _ISR_QEDMode; break; 
    case ShowerIndex::EWK : mode = _ISR_EWKMode; break; 
    default: break;
    }
  }
  return mode;
}  

bool SplittingGenerator::isFSRadiationON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  if ( isInteractionON(interaction) && isFSRadiationON() ) { 
    switch ( interaction ) {
    case ShowerIndex::QCD : mode = _FSR_QCDMode; break; 
    case ShowerIndex::QED : mode = _FSR_QEDMode; break; 
    case ShowerIndex::EWK : mode = _FSR_EWKMode; break; 
    default: break;
    }
  }
  return mode;
}


Branching SplittingGenerator::chooseForwardBranching(tPartCollHdlPtr ch,
		                                     ShowerParticle &particle,
						     const bool reverseAngOrd) const {
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "SplittingGenerator::chooseForwardBranching(): full _____________________________" << endl; 
  }
  
  Energy newQ = Energy();
  double newZ = 0.0; 
  double newPhi = 0.0;
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;

  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  if(_branchings.find(index) == _branchings.end()) 
    return Branching(ShoKinPtr(), SudakovPtr(), IdList());
  for(BranchingList::const_iterator cit = _branchings.lower_bound(index); 
      cit != _branchings.upper_bound(index); ++cit) {
    //if(cit->second.second[0] != index) continue;
    tSudakovPtr candidateSudakov = cit->second.first;
    Energy candidateNewQ = 0*MeV;
    ShowerIndex::InteractionType i = candidateSudakov->interactionType();
    if(candidateSudakov) {
	  // check size of scales beforehand...
      if(particle.evolutionScales()[i] > _showerVariables->cutoffQScale(i)) {
	  // use this condition for roughly only one gluon per quark
	  // in the asymmetric case...
// 	  if ( (particle.evolutionScales()[i] > 172.0*GeV
// 		&& particle.evolutionScales()[i] < 175.0*GeV) 
// 	       || (particle.evolutionScales()[i] > 47.7*GeV
// 		   && particle.evolutionScales()[i] < 48.05*GeV) ) {
	  // ... and this one in the symmetric case
//	  if (particle.evolutionScales()[i] > 90.0*GeV) {
        candidateNewQ = candidateSudakov->
                           generateNextBranching(particle.evolutionScales()[i], 
		    	               cit->second.second, reverseAngOrd);
		//	    }
      } else candidateNewQ = 0*GeV; 
      if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {
	    generator()->log() << "  trying "
			       << particle.data().id() << " -> ("
	                       << cit->second.second[1] << ", "
	                       << cit->second.second[2] << ") in "
			       << candidateSudakov->interactionType() << "; "
			       << "(Q0 -> Q) = (" 
			       << particle.evolutionScales()[i]/MeV
			       << " -> " 
			       << (candidateNewQ/MeV > 0 ? candidateNewQ/MeV : 0)
			       << ") MeV"  << endl;
      }	  
      if((!reverseAngOrd && candidateNewQ > newQ && 
	      candidateNewQ < particle.evolutionScales()[i]) || 
	     (reverseAngOrd && candidateNewQ < newQ && 
	      candidateNewQ > particle.evolutionScales()[i])) { 
	    newQ = candidateNewQ;
	    sudakov = candidateSudakov;
	    ids = cit->second.second;
	    newZ = sudakov->z();
	    newPhi = sudakov->phi();
      } 
    }
  }
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  " << particle.data().PDGName() << ": ";
      //		       << " [" << particle.number() << "]: ";
    if(sudakov) {
      generator()->log() << particle.data().id() << " -> ("
	                 << ids[1] << ", "
	                 << ids[2] << ") in "
			 << sudakov->interactionType() << "; "
			 << "(Q0 -> Q) = (" 
			 << particle.evolutionScales()[sudakov->interactionType()]/MeV 
			 << " -> " 
			 << (newQ/MeV > 0 ? newQ/MeV : 0)
			 << ") MeV"  << endl;    
    } else 
      generator()->log() << " won't branch." << endl; 
  }

  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  if(newQ && sudakov) {
    if(sudakov->splittingFn()) {

      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitFun = sudakov->splittingFn();
      if(splitFun) {	  
	Lorentz5Momentum p, n, ppartner, pcm;
	double th;
	Energy nnorm;
	if(particle.isFromHardSubprocess()) {
	  p = particle.momentum();
	  //***LOOKHERE*** We choose  n  naively for the time being.  
	  // Lorentz5Momentum n = Lorentz5Momentum( 0.0, - p.vect() ); 
	  // for test purposes: n points into the negative z-direction:
	  ppartner = particle.partners()[sudakov->splittingFn()->
					 interactionType()]->momentum();
	  pcm = p; 
	  pcm.boost((p + ppartner).findBoostToCM());	  
	  th = (pi);
	  nnorm = pcm.vect().mag(); 
	  nnorm*=1.;
	  Vector3 nv = cos(th)*pcm.vect().unit() + 
	    sin(th)*pcm.vect().orthogonal().unit();
	  n = Lorentz5Momentum( 0.0, nnorm*nv ); 
	  n.boost( -(p + ppartner).findBoostToCM() );
	  // n = Lorentz5Momentum( 2000.0, 3000.0, -1000.0, sqrt(14.)*1000.0); 
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	    generator()->log() << "  chosen normalization = " << nnorm 
			       << ", theta = " << th << endl 
			       << "  with p_partner = " << ppartner << endl
			       << "  and    p in cm = " << pcm << endl;
	  }
	} else {
	  p = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[0];
	  n = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[1];
	} 

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  create ShowerKinematics with " 
			     << endl 
			     << "  p = " << p << endl 
			     << "  n = " << n << endl;
	}

	Ptr<FS_QtildaShowerKinematics1to2>::pointer showerKin = 
	  new_ptr(FS_QtildaShowerKinematics1to2(p, n));

        showerKin->qtilde(newQ);
        showerKin->setResScale(sudakov->resScale());
	showerKin->setKinScale(sudakov->kinScale()); 
	showerKin->z(newZ);
	showerKin->phi(newPhi);

// 	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
// 	  generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________" << endl; 
// 	}

	return Branching(showerKin, sudakov, ids);
      }
    }
  }
  
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________" << endl; 
//   }
  
  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());

}
 

Branching SplittingGenerator::chooseBackwardBranching(tPartCollHdlPtr ch,
			                              ShowerParticle &particle) const {
  Energy newQ = Energy();
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;
  
  // First, find the eventual branching, corresponding to the highest scale.
  for(int i = 0; i < ShowerIndex::NumInteractionTypes; ++i) {
    //ShowerIndex index;
    //index.id = abs( particle.data().id() ); 
    //index.interaction = ShowerIndex::int2Interaction( i );
    //index.timeFlag = ShowerIndex::IS; 
    long index = abs(particle.data().id());
    if(_branchings.find(index) != _branchings.end()) {
      for(BranchingList::const_iterator cit = _branchings.lower_bound(index); 
	    cit != _branchings.upper_bound(index); ++cit ) {
	tSudakovPtr candidateSudakov = cit->second.first;
	Energy candidateNewQ = Energy();
	if(candidateSudakov) {
	  candidateNewQ = candidateSudakov->
		  generateNextBranching(particle.evolutionScales()[i],
				        cit->second.second);
	  if(candidateNewQ > newQ) {
	    newQ = candidateNewQ;
	    sudakov = candidateSudakov;
	    ids = cit->second.second;
	  } 
	}
      } 
    }
  }

  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  if(newQ && sudakov) {
    if(sudakov->splittingFn()) {

      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitFun = sudakov->splittingFn();
      if(splitFun) {	  

        //***LOOKHERE*** Do something similar as in chooseForwardBranching
        //               but use IS_QtildaShowerKinematics1to2 instead.

      }
    }
  }

  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());

}


void SplittingGenerator::generateBranchingKinematics(tPartCollHdlPtr ch,
	                                             ShowerParticle &particle,
					             Branching &b) const {
  //***LOOKHERE*** Complete the kinematics of the branching by filling
  //               the eventual missing bits of the ShowerKinematics
  //               object created, and already at least partially filled,
  //               by  chooseForwardBranching  or  chooseBackwardBranching.
  //               Notice that this part could remain empty if such 
  //               ShowerKinematics object is already completely filled.

}
 
void SplittingGenerator::addToMap(long &i, IdList &ids, SudakovPtr &s) {
   if(isISRadiationON(s->splittingFn()->interactionType()) || 
      isFSRadiationON(s->splittingFn()->interactionType())) {
     //i.timeFlag = ShowerIndex::IS;
      _branchings.insert(BranchingInsert(i,BranchingElement(s,ids)));
   }
}

void SplittingGenerator::debuggingInfo() {

  generator()->log() << "SplittingGenerator::debuggingInfo() begin"
		     << " ______________________________________" << endl; 
  generator()->log() << "  no of initialized Sudakov FF's = " 
		     << _branchings.size() << endl
		     << "  id\t" << "int\t" << "timeFlag" << endl;
  for(BranchingList::const_iterator cit = _branchings.begin();
      cit != _branchings.end(); ++cit) {
    generator()->log() << "  " << cit->first
		       << "\t" << cit->second.first->splittingFn()->interactionType()
      //	       << "\t" << cit->first.timeFlag 
		       << endl;
  }

  //  generator()->log() << "SplittingGenerator::debuggingInfo() end ________________________________________" << endl; 
}
