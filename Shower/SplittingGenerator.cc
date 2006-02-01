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
#include "ThePEG/PDT/EnumParticles.h"
#include "ShowerParticle.h"
#include "ShowerAlpha.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "IS_QtildaShowerKinematics1to2.h"
#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConfig.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDF/BeamParticleData.h"

using namespace Herwig;

SplittingGenerator::~SplittingGenerator() {}

void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _qcdinteractionMode << _qedinteractionMode << _ewkinteractionMode
     << _isr_Mode << _isr_qcdMode << _isr_qedMode << _isr_ewkMode
     << _fsr_Mode << _fsr_qcdMode << _fsr_qedMode << _fsr_ewkMode
     << _showerAlphaQCD << _showerAlphaQED << _showerVariables << _bbranchings
     << _fbranchings;
}

void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_qcdinteractionMode >> _qedinteractionMode >> _ewkinteractionMode
     >>	_isr_Mode >> _isr_qcdMode >> _isr_qedMode >> _isr_ewkMode 
     >> _fsr_Mode >> _fsr_qcdMode >> _fsr_qedMode >> _fsr_ewkMode
     >> _showerAlphaQCD >> _showerAlphaQED >> _showerVariables >> _bbranchings
     >> _fbranchings;
  setSVtoAlpha(_showerVariables);
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
     &SplittingGenerator::_isr_Mode, 0, false, false);
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

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQED
    ("ShowerAlphaQED", "A reference to the IS_ShowerAlphaQED object", 
     &Herwig::SplittingGenerator::_showerAlphaQED, false, false, true, false, 
     &Herwig::SplittingGenerator::setQED, 0, 0);

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQCD
    ("ShowerAlphaQCD", "A reference to the IS_ShowerAlphaQCD object", 
     &Herwig::SplittingGenerator::_showerAlphaQCD, false, false, true, false,
     &Herwig::SplittingGenerator::setQCD, 0, 0);

  static Reference<SplittingGenerator,ShowerVariables> interfaceShowerVariables
    ("ShowerVariables", "A reference to the ShowerVariables object", 
     &Herwig::SplittingGenerator::_showerVariables, false, false, true, false, 
     &Herwig::SplittingGenerator::setSVtoAlpha, 0, 0);

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
  if(!s) cerr << "Error: Could not load Sudakov " << sudakov << '\n';
  IdList ids;
  ids.push_back(parent->id());
  for(vector<tPDPtr>::iterator it = products.begin(); it!=products.end(); ++it)
    ids.push_back((*it)->id());
  addToMap(ids,s,final);
  return "";
}

Branching SplittingGenerator::chooseForwardBranching(tEHPtr ch,
		                                     ShowerParticle &particle,
						     const bool reverseAngOrd) const {
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "SplittingGenerator::chooseForwardBranching(): full _____________________________\n";
  }
  
  Energy newQ = Energy();
  double newZ = 0.0; 
  double newPhi = 0.0;
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;

  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), SudakovPtr(), IdList());
  // otherwise select branching
  for(BranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) {
    tSudakovPtr candidateSudakov = cit->second.first;
    if(candidateSudakov) {
      ShowerIndex::InteractionType i = candidateSudakov->interactionType();
      // check size of scales beforehand...
      Energy candidateNewQ(0.);
      if(particle.evolutionScales()[i] > _showerVariables->cutoffQScale(i)) {
        candidateNewQ = candidateSudakov->generateNextTimeBranching(
				       particle.evolutionScales()[i], 
		    	               cit->second.second, reverseAngOrd);
      }
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
			       << ") MeV\n";
      }	  
      // won't work for reverse angluar ordering as newQ=0. at initialization
      // NOTE: PJS - must check with stefan to find what newQ should be
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
			 << ") MeV\n";
    } else 
      generator()->log() << " won't branch.\n";
  }


  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  // NOTE: PJS - modified by durham group, removed newQ from condition
  //if(newQ && sudakov) {
  if(sudakov) {
    if(sudakov->splittingFn()) {

      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitFun = sudakov->splittingFn();

      // find out whether we had a spacelike branching before. 
      
      bool fromHard = particle.isFromHardSubprocess();
//       cerr << particle.parents()[0]
// 	   << '\n'; 
//       cerr << dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
// 	   << '\n'; 
      
      ShowerParticlePtr partest = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0]); 

      bool parentHas = (partest && partest->showerKinematics());
      partest = dynamic_ptr_cast<ShowerParticlePtr>(*(particle.siblings().begin()));
      bool siblingHas = (partest && partest->showerKinematics());
      bool initiatesTLS = (particle.initiatesTLS());

      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
       cout << "fromHard, parentHas, siblingHas, initiatesTLS = " 
 	   << (fromHard ? "y":"n") << ", " 
 	   << (parentHas ? "y":"n") << ", "
 	   << (siblingHas ? "y":"n") << ", " 
 	   << (initiatesTLS ? "y":"n") << '\n';
      }
//       if (!(parentHas || fromHard) && siblingHas) {
// 	cout << "Partner...\n";
// 	cout << particle.partners()[splitFun->interactionType()] << ", "
// 	     << particle.partners()[splitFun->interactionType()]->id()
// 	     << particle.partners()[splitFun->interactionType()]->momentum()
// 	     << '\n'; 
//    }

      if(splitFun) {	  
	
	//	cerr << "B\n";
	Lorentz5Momentum p, n, ppartner, pcm;
	double th;
	Energy nnorm;
	//	if(particle.isFromHardSubprocess()) {
	//	if (fromHard || siblingHas && !parentHas) {
	if (fromHard || initiatesTLS) {
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
			       << ", theta = " << th 
			       << "\n  with p_partner = " << ppartner
			       << "\n  and    p in cm = " << pcm << '\n';
	  }
	} else {

	  //	  cerr << "A\n";
	  //for(tParticleSet::const_iterator cit = particle.siblings().begin(); 
	  //    cit != particle.siblings().end(); cit++) {
	  //  cout << (*cit)->id() << ", ";
	  //}
	  p = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[0];
	  n = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[1];
	} 


	//	cerr << "C\n";
  

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  create ShowerKinematics with\n" 
			     << "  p = " << p
			     << "\n  n = " << n << '\n';
	}

	Ptr<FS_QtildaShowerKinematics1to2>::pointer showerKin = 
	  new_ptr(FS_QtildaShowerKinematics1to2(p, n));

        showerKin->qtilde(newQ);
        showerKin->setResScale(sudakov->resScale());
	showerKin->setKinScale(sudakov->kinScale()); 
	showerKin->z(newZ);
	showerKin->phi(newPhi);

// 	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
// 	  generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________\n";
// 	}

	return Branching(showerKin, sudakov, ids);
      }
    }
  }
  
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________\n";
//   }
  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());

}


Branching SplittingGenerator::chooseBackwardBranching(tEHPtr ch,
			                              ShowerParticle &particle) const {
  Energy newQ = Energy();
  tSudakovPtr sudakov = tSudakovPtr();
  IdList ids;
  double newZ(0.0), newPhi(0.0);  

  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.id());
  double x = particle.x();
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Called cBB with " << &particle << ", id = " 
       << particle.id() << " and x = " << x << '\n';
  }
  if(_bbranchings.find(index) == _bbranchings.end()) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  no branchings...\n";
    }
    return Branching(ShoKinPtr(), tSudakovPtr(), IdList());
  }
  for(BranchingList::const_iterator cit = _bbranchings.lower_bound(index); 
      cit != _bbranchings.upper_bound(index); ++cit ) {
    tSudakovPtr candidateSudakov = cit->second.first;
    Energy candidateNewQ = Energy();
    if(candidateSudakov) {
      ShowerIndex::InteractionType i = candidateSudakov->interactionType();
      // The parton we are IS radiating should have one parent which is
      // the beam particle. If this isn't the case, we must be
      // backward evolving for a decay shower. In that case give a null
      // pdf pointer
//       cout << "A ";
//       cout << "particle.parents().size() = " << particle.parents().size()
// 	   << endl;
//       cout << ", particle.parents()[0] = " 
// 	   << particle.parents()[0] << endl;
//      cerr << "A\n";
//      cerr << "particle.parents().size() " 
//           << particle.parents().size() << '\n';
   
      // at least should be outside loop perhaps even outside function
      // and in backward evolver
   
      Ptr<BeamParticleData>::tcp p = 
	dynamic_ptr_cast<Ptr<BeamParticleData>::tcp>(particle.parents()[0]
						     ->dataPtr());
//      cerr << "B" << endl;;
      tcPDFPtr pdf;
      if(p) pdf = p->pdf();
      else pdf = tcPDFPtr();
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Looking at splitting " 
			   << cit->second.second[0] << "->" 
			   << cit->second.second[1] << ","
			   << cit->second.second[2] << ";"
			   << " Q = " << particle.evolutionScales()[i]/GeV 
			   << '\n';
      }
      candidateNewQ = candidateSudakov->
        generateNextSpaceBranching(particle.evolutionScales()[i],
				   cit->second.second, p, pdf, particle.x());
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  found Splitting with q = " 
			   << candidateNewQ/GeV 
			   << " [GeV]\n";
      }
      if(candidateNewQ > newQ) {
	newQ = candidateNewQ;
	sudakov = candidateSudakov;
	ids = cit->second.second;
	newZ = sudakov->z();
	newPhi = sudakov->phi();
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
	Lorentz5Momentum p, n, pthis, ppartner, pcm;
	if(particle.isFromHardSubprocess()) {
// 	  cout << particle.parents().size() 
// 	       << " parents. parents()[0]->id() = " 
// 	       << particle.parents()[0]->id() << endl;

// 	  pthis = particle.momentum();
// 	  ppartner = particle.partners()[sudakov->splittingFn()->
// 					 interactionType()]->momentum();
// 	  pcm = pthis; 
// 	  pcm.boost((pthis + ppartner).findBoostToCM());	  
	  pcm = particle.parents()[0]->momentum();
	  p = Lorentz5Momentum(0.0, pcm.vect());
	  n = Lorentz5Momentum(0.0, -pcm.vect()); 
// 	  p.boost( -(pthis + ppartner).findBoostToCM() );
// 	  n.boost( -(pthis + ppartner).findBoostToCM() );
	} else {
	  p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
	    ->showerKinematics()->getBasis()[0];
	  n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
	    ->showerKinematics()->getBasis()[1];
	} 
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  create ShowerKinematics with\n" 
			     << "  p = " << p
			     << "\n  n = " << n << '\n';
	}

	Ptr<IS_QtildaShowerKinematics1to2>::pointer 
	  showerKin = new_ptr(IS_QtildaShowerKinematics1to2(p, n));

        showerKin->qtilde(newQ);
        showerKin->setResScale(sudakov->resScale());
	showerKin->setKinScale(sudakov->kinScale()); 
	showerKin->z(newZ);
	showerKin->phi(newPhi);

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "cBB found branching, shoKin = "
			     << showerKin
			     << "\n  Q = " << newQ/GeV 
			     << ", z = " << newZ 
			     << ", phi = " << newPhi << '\n';
	}
	return Branching(showerKin, sudakov, ids);

      }
    }
  }

  return Branching(ShoKinPtr(), tSudakovPtr(), IdList());

}


void SplittingGenerator::generateBranchingKinematics(tEHPtr ch,
	                                             ShowerParticle &particle,
					             Branching &b) const {
  //***LOOKHERE*** Complete the kinematics of the branching by filling
  //               the eventual missing bits of the ShowerKinematics
  //               object created, and already at least partially filled,
  //               by  chooseForwardBranching  or  chooseBackwardBranching.
  //               Notice that this part could remain empty if such 
  //               ShowerKinematics object is already completely filled.

}

void SplittingGenerator::addToMap(IdList &ids, SudakovPtr &s, bool final) {
   if(isISRadiationON(s->splittingFn()->interactionType()) && !final) {
  
      _bbranchings.insert(BranchingInsert(ids[1],BranchingElement(s,ids)));
   } 
   if(isFSRadiationON(s->splittingFn()->interactionType()) && final) {
     //i.timeFlag = ShowerIndex::IS;
      _fbranchings.insert(BranchingInsert(ids[0],BranchingElement(s,ids)));
   }
}

tSplittingFnPtr SplittingGenerator::getSplittingFunction(long id1, long id2, 
							 bool init) {
  BranchingList *bl;
  BranchingList::const_iterator cit;
  int idx;
  if(init) { bl = &_bbranchings; idx = 0; }
  else { bl = &_fbranchings; idx = 1; }
  long index = abs(id1);
  if(bl->find(index) == bl->end()) return tSplittingFnPtr();
  for(cit = bl->lower_bound(index); cit != bl->upper_bound(index); ++cit) {
    if(cit->second.second[idx] == abs(id2)) 
      return cit->second.first->splittingFn();
  }
  return tSplittingFnPtr();
}

void SplittingGenerator::debuggingInfo() {

  generator()->log() << "SplittingGenerator::debuggingInfo() begin"
		     << " ______________________________________\n";
  generator()->log() << "  no of initialized Sudakov FF's = " 
		     << _fbranchings.size() << endl
		     << "  id\t" << "int\t" << "timeFlag\n";
  for(BranchingList::const_iterator cit = _fbranchings.begin();
      cit != _fbranchings.end(); ++cit) {
    generator()->log() << "  " << cit->first
		       << "\t" << cit->second.first->splittingFn()->interactionType()
      //	       << "\t" << cit->first.timeFlag 
		       << '\n';
  }

  //  generator()->log() << "SplittingGenerator::debuggingInfo() end ________________________________________" << endl; 
}

