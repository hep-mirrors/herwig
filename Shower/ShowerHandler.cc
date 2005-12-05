// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
//#include "ThePEG/Handlers/CollisionHandler.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDF/PDF.h"
#include "ShowerConfig.h"
#include "ShowerIndex.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/PDT/ParticleData.h"
// #include "ThePEG/MatrixElement/PhaseSpaceBase.h"
#include "ThePEG/Handlers/XComb.h"
#include "ShowerParticle.h"
#include "ThePEG/EventRecord/ColourLine.h"
//#include "ShowerColourLine.h"
#include "ThePEG/CLHEPWrap/Matrix.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/EventRecord/Step.h"
#include "QQbarG.h"

using namespace Herwig;


ShowerHandler::~ShowerHandler() {}


void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _globalParameters 
    //<< _MECorrections 
     << _showerVariables
     << _evolver;
}


void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _globalParameters 
    //>> _MECorrections 
     >> _showerVariables
     >> _evolver;
}


ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.


void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.");

  static Reference<ShowerHandler,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
			      "A reference to the GlobalParameters object", 
			      &Herwig::ShowerHandler::_globalParameters,
			      false, false, true, false);
  /*static Reference<ShowerHandler,MECorrections> 
    interfaceMECorrections("MECorrections", 
			   "A reference to the MECorrections object", 
			   &Herwig::ShowerHandler::_MECorrections,
			   false, false, true, false);*/
  static Reference<ShowerHandler,ShowerVariables> 
    interfaceShowerVariables("ShowerVariables", 
			       "A reference to the ShowerVariables object", 
			       &Herwig::ShowerHandler::_showerVariables,
			       false, false, true, false);
  static Reference<ShowerHandler,Evolver> 
    interfaceEvolver("Evolver", 
		     "A reference to the Evolver object", 
		     &Herwig::ShowerHandler::_evolver,
		     false, false, true, false);

}


void ShowerHandler::cascade() {


//   cerr << "=========================================" 
//        << "======================================\r"
//        << "ShowerHandler::cascade() started for event no " 
//        << generator()->currentEventNumber() << endl;


  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "ShowerHandler::debuggingInfo() Evt #" 
		       << generator()->currentEventNumber() 
		       << "\t _______________________________________" << endl;
  }

  tEHPtr ch = eventHandler();

  //Particle::PrintParticles(cout, ch->currentStep()->particles().begin(),
  //		   ch->currentStep()->particles().end());
  // const Hint & theHint = hint();	    // OK: COMMENTED TO AVOID COMPILATIO+N WARNINGS
  // const PDF & theFirstPDF = firstPDF();  // OK: COMMENTED TO AVOID COMPILATION+ WARNINGS
  // const PDF & theSecondPDF = secondPDF();// OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const pair<PDF,PDF> & thePdfs = pdfs();// OK: COMMENTED TO AVOID COMPILATION WARNINGS

  // From the ThePEG particles entering the hard subprocess, create
  // the corresponding starting ShowerParticle objects and put them
  // in the vector hardProcessParticles (not directly in _particles
  // because we could throw away the all showering and restart it).
  // ShowerParticleVector hardProcessParticles;
  //  createShowerParticlesFromP7Particles( ch, hardProcessParticles );

  static const int maxNumFailures = 100;
  int countFailures = 0;
  bool keepTrying = true;
  
  // tie in hard corrections here
  if(_showerVariables->hardMEC()) hardMEC(ch);
  bool evolverOK = true;

  // even more cleaning... ***ACHTUNG*** find something much, much
  // more sophisticated to save old properties. In particular for
  // the multiscaleshower...
  ShowerParticleVector hardProcessParticles;
  hardProcessParticles.clear();
  convertToShowerParticles( ch, hardProcessParticles );
  
  while(keepTrying && countFailures < maxNumFailures) {       
    try {

      // printing the current status of _particles
      //      cerr << "Beginning of Shower loop" << endl;
//       for(ShowerParticleVector::const_iterator cit = _particles.begin();
// 	  cit != _particles.end(); ++cit) {
// 	cerr << (*cit)->id() << "(" << (*cit) << ")" << " p: ";
// 	if((*cit)->parents().size()) {
// 	  cerr << (*cit)->parents()[0]->id() 
// 	       << "(" << (*cit)->parents()[0] << ")" << " c: ";
// 	}
// 	for (unsigned int j=0; j < (*cit)->children().size(); j++) {
// 	  cerr << (*cit)->children()[j]->id() 
// 	       << "(" << (*cit)->children()[j] << ")" << " ";
// 	} 
// 	cerr << endl;
//       }


      // Cleaning and initializing
      _particles.clear();

//       cerr << "after _particles.clear() " << endl;
//       for(ShowerParticleVector::const_iterator cit = _particles.begin();
// 	  cit != _particles.end(); ++cit) {
// 	cerr << (*cit)->id() << "(" << (*cit) << ")" << " p: ";
// 	if((*cit)->parents().size()) {
// 	  cerr << (*cit)->parents()[0]->id() 
// 	       << "(" << (*cit)->parents()[0] << ")" << " c: ";
// 	}
// 	for (unsigned int j=0; j < (*cit)->children().size(); j++) {
// 	  cerr << (*cit)->children()[j]->id() 
// 	       << "(" << (*cit)->children()[j] << ")" << " ";
// 	} 
// 	cerr << endl;
//       }
      _showerVariables->reset();     
      _evolver->clear();
      
      // Fill _particles with the particles from the hard subprocess
      ShowerParticleVector::const_iterator cit;
      for(cit = hardProcessParticles.begin(); 
	  cit != hardProcessParticles.end(); ++cit ) {
	if(*cit) { 
	  _particles.push_back( *cit );      
	}
      }

//       cerr << "after _particles refill " << endl;
//       for(ShowerParticleVector::const_iterator cit = _particles.begin();
// 	  cit != _particles.end(); ++cit) {
// 	cerr << (*cit)->id() << "(" << (*cit) << ")" << " p: ";
// 	if((*cit)->parents().size()) {
// 	  cerr << (*cit)->parents()[0]->id() 
// 	       << "(" << (*cit)->parents()[0] << ")" << " c: ";
// 	}
// 	for (unsigned int j=0; j < (*cit)->children().size(); j++) {
// 	  cerr << (*cit)->children()[j]->id() 
// 	       << "(" << (*cit)->children()[j] << ")" << " ";
// 	} 
// 	cerr << endl;
//       }

      // clean up children from previous trial showers:
      //      for ( ShowerParticleVector::iterator it = _particles.begin(); 
      // 	    it != _particles.end(); ++it) {
      //	(*it)->removeChildren(); 
      //}
      
      /*bool isMECorrectionApplied = false;
      tMECorrectionPtr softMECorrection = tMECorrectionPtr();
      if ( _MECorrections->isMECorrectionsON() &&
	   _MECorrections->getMECorrection( ch->lastME() ) ) { 
        isMECorrectionApplied = true;
	if ( _MECorrections->getMECorrection( ch->lastME() )->hardProcessPlusJetME() ) {
	  _MECorrections->getMECorrection( ch->lastME() )->hardMECorrection();
	} else { 
	  softMECorrection = _MECorrections->getMECorrection( ch->lastME() );
	}
	}*/

      if(_showerVariables->isDecayBeforeShowerON()) {      
        // Loop over _particles and for those that satisfy the condition
        // if ( _showerVariables->hasToDecayBeforeShower( id ) ) 
        // then decay them. Apply the same procedure iteratively for
        // all decay products, until none particles in _particles
        // satisfies the above condition.
        // ***LOOKHERE*** --- code here --- 
      }
      Energy largestWidth = Energy();
      if(_showerVariables->isMultiScaleShowerON()) {      
        for(ShowerParticleVector::const_iterator cit = _particles.begin();
	    cit != _particles.end(); ++cit) {
	  if((*cit)->children().size()) {
	    if((*cit)->data().width() > largestWidth) { 
	      largestWidth = (*cit)->data().width();
	    }
	  }
	}
	if(largestWidth < _globalParameters->hadronizationScale()) {  
	  largestWidth = Energy();
	}
      }

      // In the case there are no more shower steps, we postpone 
      // the final kinematics reconstruction after we have decided
      // whether the gluons should be on their physical massless 
      // shell or on the effective mass shell. This allows us to
      // avoid to do twice the same kinematics reconstruction, 
      // the first time with massless gluons, and the second one 
      // with gluons on their effective mass shell.
      bool skipKinReco = false;
      if(largestWidth < _globalParameters->hadronizationScale()) {  
	skipKinReco = true;
      }      
      
      _showerVariables->stopShowerAtMassScale(largestWidth);
      // ***MAYBE NOT NEEDED***
      _showerVariables->vetoBelowPtScale(largestWidth);  
      //      cerr << *(ch->currentStep());

      evolverOK = _evolver
	->showerNormally(ch, _showerVariables, _particles, skipKinReco);
      if (!evolverOK) {
	// 	cerr << "ShowerHandler: evolver not ok, jumping out of big while loop." << endl;
	
	keepTrying = true; 
	continue;
      }
      // The following loop is done only for multi-scale showering.
      while(largestWidth > _globalParameters->hadronizationScale()) {  
	
        // Because is forbidden to add elements to a STL container while
	// looping over it, we have to create first a new container which
        // holds the particles that need to be considered, that is the 
        // current initial and final states; then we can loop over the
        // new container, and adding to the initial one the new particles
        // that are produced by the showering.
	
        ShowerParticleVector currentFinalOrInitialParticles;
        for(cit = _particles.begin(); cit != _particles.end(); ++cit) {
	  if((*cit)->children().size() == 0) {
	    currentFinalOrInitialParticles.push_back(*cit);
	  }
	}
	ShowerParticleVector::iterator it;
        for(it = currentFinalOrInitialParticles.begin();
	    it != currentFinalOrInitialParticles.end(); ++it) {
	  Energy particleWidth = (*it)->data().width();
	  if(particleWidth >= largestWidth) {
            ShowerParticleVector decayParticles;
	    // Decay the particle *it (need pointer to Decayer)
            // and put the decaying ShowerParticle object and all its
            // decay products (for each of them we need to create a
            // ShowerParticle object) in the collection decayParticles.
            // ***LOOKHERE*** --- code here ---
	    
	    /* tMECorrectionPtr softMEDecayCorrection = tMECorrectionPtr();
	       if ( _MECorrections->isMECorrectionsON()  &&
	       ( ( ! isMECorrectionApplied ) || 
	       ( isMECorrectionApplied  &&  _MECorrections->isComposeMECorrectionsON() )
	       ) ) {
	       // if ( _MECorrections->getMECorrection( "decayer pointer" ) {
	       //    isMECorrectionApplied = true;
	       //    if ( _MECorrections->
	       //         getMECorrection( "decayer pointer" )->decayProcessPlusJetME() ) {
	       //      _MECorrections->getMECorrection( ch->lastME() )->hardMECorrection();
	       //    } else { 
	       //      softMEDecayCorrection = _MECorrections->
	       //                              getMECorrection( "decayer pointer" );
	      //    }
              // }
	      }*/
	    
            // In order to proper treat exceptional cases in which 
            // a decay product has larger width than the decaying parent, 
            // we have to keep track of the largest width for each
            // decaying system (including the decaying parent), and
            // also of the global largest width for all possible
            // decaying systems.  
            Energy largestWidthDecayingSystem = Energy();
	    for(cit = decayParticles.begin();
		cit != decayParticles.end(); ++cit) {
	      if((*cit)->data().width() > largestWidthDecayingSystem) { 
		largestWidthDecayingSystem = (*cit)->data().width();
	      }
	    }
	    if(largestWidthDecayingSystem > largestWidth) {
	      largestWidth = largestWidthDecayingSystem;
	    }
	    
	    _showerVariables->stopShowerAtMassScale(largestWidthDecayingSystem); 
	    //***MAYBE NOT NEEDED***
	    _showerVariables->vetoBelowPtScale(largestWidthDecayingSystem);  
   
	    _evolver->showerDecay(ch, _showerVariables, decayParticles);
            
	    // Add the particles in  decayParticles  in  _particles .
            // Notice that the first particle in  decayParticles  must be
            // skipped because it corresponds to the decaying particle
            // which was already in  _particles .
	    ShowerParticleVector::const_iterator cit = decayParticles.begin();
            do {
	      ++cit;
	      _particles.insert( _particles.end(), *cit );
	    } while ( cit != decayParticles.end() );
	    decayParticles.clear();
          }
	  
        } // end for loop
	
        // Find now the new largest width between all current final 
	// state particles.
	Energy newLargestWidth = Energy();
        for(cit = _particles.begin();
	    cit != _particles.end(); ++cit) {
	  if((*cit)->children().size()) {
	    if((*cit)->data().width() > newLargestWidth) { 
	      newLargestWidth = (*cit)->data().width();
	    }
	  }
	}
	if(newLargestWidth < _globalParameters->hadronizationScale()) {  
	  newLargestWidth = Energy();
	  // to avoid to do twice the kinematics reco
	  skipKinReco = true;          
	}                              

	Energy savedVetoAbovePtScale = _showerVariables->vetoAbovePtScale();
	
	_showerVariables->stopShowerAtMassScale(newLargestWidth);
	// ***MAYBE NOT NEEDED***
	_showerVariables->vetoBelowPtScale(newLargestWidth); 
	_showerVariables->vetoAbovePtScale(largestWidth);
	
	_evolver->showerGlobally(ch,_showerVariables,_particles,skipKinReco);
	
	_showerVariables->vetoBelowPtScale( savedVetoAbovePtScale );
	largestWidth = newLargestWidth;

      } // end of while loop (multi-scale loop)

      // In the case Herwig++ Cluster Hadronization model is used for the
      // hadronization, set all final state gluons on the effective mass
      // shell (rather on the physical massless shell).
      if(!_globalParameters->isThePEGStringFragmentationON()) {
	_evolver->setEffectiveGluonMass
	  (_globalParameters->effectiveGluonMass() , _particles);
      } else { 
	_evolver->setEffectiveGluonMass(0 , _particles);
      }

      // Do the final kinematics reconstruction
      bool recons = _evolver->reconstructKinematics(ch);
      
      // set true anyway for test purposes...
      recons = true;


      // Fill the positions information for all the ShowerParticle
      // objects in _particles. It is done at this stage in order to
      // avoid the complications of boosting such positions during the
      // kinematics reshuffling.  ***ACHTUNG!*** fillPostitions()
      // gives segfault when setting labVertex. 
      // fillPositions();

      keepTrying = false;
      if (!recons) {
	if (countFailures++ < maxNumFailures) keepTrying = true; 
	else throw Exception::eventerror; 	
      }
    } // end of try
    catch (Veto &v) {cout << "throwing again!" << endl; throw v;}
    catch ( std::exception & e ) {
      countFailures++;
    }    
  } // end main while loop

  // Debugging
  fillEventRecord(ch);  

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower ) {    
    debuggingInfo();
    if ( generator()->currentEventNumber() < 1000 )
      generator()->log() << "ShowerHandler::debuggingInfo "
			 << " ===> END DEBUGGING <=== "
			 << endl;
  }
//  cerr << *(ch->currentStep()) << endl;
}

void ShowerHandler::
convertToShowerParticles(const tEHPtr ch,
			 ShowerParticleVector & hardProcessParticles) {

  // Incoming (initial state) particles (indeed the partons entering
  // the hard subprocess, not the beam hadrons). 
  //Particle::PrintParticles(cout, ch->currentStep()->particles().begin(),
  //			   ch->currentStep()->particles().end());
  ShowerParticlePtr part;
  if(ch->lastPartons().first ) {
    part = ptr_new<ShowerParticlePtr>(*ch->lastPartons().first);
    if ( part ) {
      part->setFromHardSubprocess(true);
      part->setThePEGBase(ch->lastPartons().first); 
      part->setFinalState(false);
      part->x(ch->lastX1());
      hardProcessParticles.push_back(part);
    }
  } 
  if(ch->lastPartons().second ) {
    part = ptr_new<ShowerParticlePtr>(*ch->lastPartons().second);
    if ( part ) {
      part->setFromHardSubprocess(true);
      part->setThePEGBase(ch->lastPartons().second); 
      part->setFinalState(false);
      part->x(ch->lastX2());
      hardProcessParticles.push_back(part); 
    }
  }
  // Outgoing (final state) particles, excluding the beam remnants. 
  // Notice that we don't need to set true the flag ShowerParticle::isFinalState 
  // because it is set true by default. 
  tParticleSet remnantSet = ch->currentCollision()->getRemnants();
  for ( ParticleSet::const_iterator cit = ch->currentStep()->particles().begin();
	cit != ch->currentStep()->particles().end(); ++cit ) {
    if ( remnantSet.find( (*cit)->original() ) == remnantSet.end() ) {
      part = ptr_new<ShowerParticlePtr>(**cit);
      part->setFromHardSubprocess(true);      
      part->setThePEGBase(*cit); 
      hardProcessParticles.push_back(part);
    }
  }    

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "  Total number of initial Shower particles : " 
		       << hardProcessParticles.size() << endl;
      //<< "  Total number of colour lines : " 
      //	       << mapP7toShowerColine.size() << endl;
    Lorentz5Momentum p_in = Lorentz5Momentum();
    Lorentz5Momentum p_out = Lorentz5Momentum();
    for ( ShowerParticleVector::const_iterator cit = hardProcessParticles.begin();
	  cit != hardProcessParticles.end(); ++cit ) {
      generator()->log() << "\t" << (**cit).data().PDGName() 
			 << "  p=" << (**cit).momentum()
			 << ( (**cit).isFinalState() ? ", out" : ", in" )
			 << endl;
      if ((**cit).isFinalState()) p_out += (**cit).momentum(); 
      else p_in += (**cit).momentum(); 
    }
    generator()->log() << "  p_in  = " << p_in << endl
		       << "  p_out = " << p_out << endl
		       << "  diff  = " << p_in - p_out << endl; 
  }  
}


void ShowerHandler::fillPositions() {

  // Put on the vector  particlesToExamine  the particles that come
  // directly from the hard subprocess. 
  ParticleVector particlesToExamine;
  for(ShowerParticleVector::const_iterator cit = _particles.begin();
      cit != _particles.end(); ++cit) {
    if((*cit)->isFromHardSubprocess()) {
      particlesToExamine.push_back(*cit);
    }
  }

  // Loop over the vector  particlesToExamine , consider its last element
  // and removing it from the vector, then update the position of such
  // element (previous position + calculated displacement), and finally 
  // insert in the vector  particlesToExamine  its children and set their
  // initial positions to the current parent position,
  // until the vector  particlesToExamine  is empty.
  while ( ! particlesToExamine.empty() ) {

    PPtr particle = particlesToExamine.back();
    particlesToExamine.pop_back();

    // See Mark Smith's thesis, page 121, formula (4.106), for the  
    // calculation of the displacement of the particle.
    // Notice that the displacement is already expressed in the Lab frame,
    // therefore no boost is necessary.
    double lambda;
    if ( !particle->data().stable() ) { // unstable, decayed particle
      lambda = particle->momentum().m() / particle->data().width();
    } else {                     // not decayed particle
      lambda = particle->momentum().mag2() / 
	                 _globalParameters->minVirtuality2(); 
    }
    LorentzDistance 
      distance(_globalParameters->conversionFactorGeVtoMillimeter() 
	       * particle->momentum().vect() / GeV, 
	       _globalParameters->conversionFactorGeVtoMillimeter() 
	       * particle->momentum().e() / GeV );
    distance *= log(1.0/rnd()) / 
      (sqrt(sqr(particle->momentum().mag2() - particle->momentum().m2()) +
	    sqr(particle->momentum().mag2() / lambda)) / GeV);

    // update the position.
    LorentzDistance temp = particle->labVertex() + distance; 
    particle->setLabVertex(temp); 
         
    for(ParticleVector::const_iterator cit = particle->children().begin();
	cit != particle->children().end(); ++cit) {
      LorentzDistance temp = particle->labVertex(); 
      (*cit)->setLabVertex(temp); // initialize the position.
      particlesToExamine.push_back(*cit);
    }
  }
}


void ShowerHandler::debuggingInfo() {
  
  if ( generator()->currentEventNumber() < 1000 )
    generator()->log() << "ShowerHandler::debuggingInfo(): ________________________________________________" << endl; 
  if ( generator()->currentEventNumber() == 1 ) {
    // Preliminary information to print only once.
    generator()->log() << "  Info for GlobalParameters:" << endl
		       << "  isOnStringFrag = " 
		       << ( _globalParameters->isThePEGStringFragmentationON() ? 
			    "YES" : "NO" ) 
		       << endl
		       << "  effective gluon mass = " 
		       << _globalParameters->effectiveGluonMass() / GeV 
		       << endl
		       << "  hadronization scale = " 
		       << _globalParameters->hadronizationScale() / GeV 
		       << "  [GeV] " << endl
		       << "  Info for ShowerVariables: " << endl
		       << "  switches:  multi-scale = " 
		       << _showerVariables->isMultiScaleShowerON()
		       << ", decay before = " 
		       << _showerVariables->isDecayBeforeShowerON() 
		       << ", gluino? = " 
		       << ( _showerVariables->hasToDecayBeforeShower(1000021) ? 
			    "Y" : "N" )  
		       << endl
		       << "  cutoff-parameters:  QCD = "
		       << _showerVariables->cutoffMassScale(ShowerIndex::QCD) / GeV 
		       << ", QED = "  
		       << _showerVariables->cutoffMassScale(ShowerIndex::QED) / GeV 
		       << ", EWK = "  
		       << _showerVariables->cutoffMassScale(ShowerIndex::EWK) / GeV 
		       << " [GeV] " 
		       << endl; 
  }

  if ( generator()->currentEventNumber() < 1000 )
    generator()->log() << "  Shower finished - summary:" << endl;
  
  tShowerParticleVector fs;
  for (ShowerParticleVector::const_iterator cit = _particles.begin(); 
       cit != _particles.end(); ++cit ) {
    if ( (*cit)->isFromHardSubprocess() ) {
      if ( generator()->currentEventNumber() < 1000 ) {
	generator()->log() << "  " << (*cit)->data().PDGName() 
			   << ", m = " << (*cit)->momentum().mass()/MeV 
			   << ", Q = " << (*cit)->momentum().m()/MeV 
			   << ", q = " << (*cit)->momentum()/MeV << " (MeV)"
			   <<endl;
      }
      if ( (*cit)->isFinalState() ) {
	tShowerParticleVector thisfs = (*cit)->getFSChildren();
	fs.insert(fs.end(), thisfs.begin(), thisfs.end()); 
      }
    }    
  }
}


void ShowerHandler::fillEventRecord(const tEHPtr ch) {  

  // Transform some of the ShowerParticles object in ThePEG particles,
  // set properly the parent/child relationships and treat carefully 
  // the transformation from ShowerColourLine objects into ThePEG
  // ColourLine ones; and finally then write them in the Event Record     
  StepPtr pstep;
  pstep = ch->newStep(); 
  //  pstep = ch->currentStep(); 
  ShowerParticleVector::iterator it;
  for(it = _particles.begin(); it != _particles.end(); ++it ) {
    PPtr p = dynamic_ptr_cast<PPtr>(*it);
    if((*it)->isFromHardSubprocess() && (*it)->isFinalState()) {
      pstep->setCopy(p, (*it)->getThePEGBase());
      //      pstep->setCopy((*it)->getThePEGBase(), p);
      // (*it)->setMomentum(p->momentum());
      addFinalStateShower((*it),pstep);
      // pstep->addDecayProduct(p);
    } else if((*it)->isFromHardSubprocess()) { 
      pstep->insertIntermediate(p, (*it)->getThePEGBase()->parents()[0],
				(*it)->getThePEGBase()->children()[0]);
      // Otherwise it is a initial state shower particle
      //      pstep->setCopy((*it)->getThePEGBase(), p);
      //      pstep->setCopy(p,(*it)->getThePEGBase());
      //      (*it)->setMomentum(p->momentum());
      //      addInitialStateShower(p,pstep,false);
      addInitialStateShower(p,pstep,false);
    } 
  }  
}


void ShowerHandler::addFinalStateShower(ShowerParticlePtr &p, StepPtr &s) {
  tPPtr dum; 
  tParticleVector yet; 
  tParticleVector addCh;
  ParticleVector::const_iterator cit;
  yet.push_back(p);  
  while(!yet.empty()) { 
    dum = yet.back(); 
    yet.pop_back(); 
    for(cit = dum->children().begin(); cit != dum->children().end(); ++cit) { 
      yet.push_back(*cit); 
      addCh.push_back(*cit); 
    }
    while(!addCh.empty()) {           
      s->addDecayNoCheck(dum, addCh.back());
      addCh.pop_back(); 
    }
  }  
}

// Recursive function
void ShowerHandler::addInitialStateShower(PPtr &p, StepPtr &s, bool doit) {
  // Each parton here should only have one parent
  if(p->parents().size()) { 
    // we have parents, call recursively, should only have one parent however
    PPtr parent = const_ptr_cast<PPtr>(p->parents()[0]);
    addInitialStateShower(parent,s);
  }
  if(doit) {
    for(unsigned int i = 0; i<p->children().size(); i++) {
      if(p->children()[i]->children().size()) {
	// Add it only as an intermediate, if it is a shower intermediate
	// add the child is a final state shower particle,
	// then we add the final state shower from the children
	s->addIntermediate(p->children()[i]);
	ShowerParticlePtr shoChild = 
	  dynamic_ptr_cast<ShowerParticlePtr>(p->children()[i]);
	if(shoChild && shoChild->isFinalState())
	  addFinalStateShower(shoChild,s);
      } else {
	s->addDecayProduct(p->children()[i]);
      }
    }
  }
}


void ShowerHandler::printStep(tStepPtr ptrStep, const string & title) {
  generator()->log() << "[[[ ShowerHandler::printStep " << title << endl;
  for ( ParticleSet::const_iterator it = ptrStep->particles().begin();
  	 it != ptrStep->particles().end(); ++it ) {
    generator()->log() << "  " << (*it)->data().PDGName() 
		       << " " << (*it)->number() 
		       << "<"
		       << (*it)->parents().size()
		       << "> ("
		       << (*it)->children().size()
		       << ") [" << ( (*it)->colourNeighbour() ? 
				     (*it)->colourNeighbour()->number() : 0 ) 
		       << "," << ( (*it)->antiColourNeighbour() ? 
				   (*it)->antiColourNeighbour()->number() : 0 )
                       << "] " << (*it)->momentum().m() 
		       << ", " << (*it)->momentum() 
		       << endl;     
  }
  generator()->log() << "]]] ShowerHandler::printStep " << title << endl;
}



void ShowerHandler::hardMEC(const tEHPtr ch) {
  PVector qq; 
  for (ParticleSet::const_iterator cit=ch->currentStep()->particles().begin();
       cit != ch->currentStep()->particles().end(); ++cit )
    if (abs((*cit)->id()) < 7) qq.push_back(*cit);  
  if (qq.size() == 2) {
    Energy Q = (qq[0]->momentum() + qq[1]->momentum()).m();
    QQbarG *qqg = new QQbarG(Q, qq[0]->momentum().m());
    vector<Lorentz5Momentum> newfs = qqg->applyHard(qq); 
    if(newfs.size() == 3) {
      bool check = true; 
      for (int i=0; i<2; i++)
	if (newfs[i].e() < qq[i]->data().constituentMass()) 
	  check = false; 
      if (newfs[2].e() < _globalParameters->effectiveGluonMass())
	check = false; 
      if (!check) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) { 
	  generator()->log() << "ShowerHandler::hardMEC: " 
			     << "3jet particles too soft to continue!"
			     << endl;
	}
 	newfs.clear();
      } else {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) { 
	  generator()->log() << "ShowerHandler::hardMEC: " 
			     << "replacing hard FS, Q = " 
			     << Q/GeV << "" << endl;
	} 
	for (int i=0; i<3; i++)
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower )  
	    generator()->log() << "  " << newfs[i] 
			       << ", m = " 
			       << newfs[i].m() << endl; 
	// incoming qqbar decay into the final (q qbar g) in current step
	for (int i=0; i<2; i++) 
	  newfs[i].setMass(qq[i]->data().constituentMass());
	newfs[2].setMass(_globalParameters->effectiveGluonMass());
	PPtr newq = getParticleData(
				    abs(qq[0]->id()))->produceParticle(newfs[0]);
	PPtr newa = getParticleData(
				    -abs(qq[0]->id()))->produceParticle(newfs[1]);
	PPtr newg = getParticleData(21)->produceParticle(newfs[2]);
	newq->antiColourNeighbour(newg);
	newa->colourNeighbour(newg);
	for (PVector::iterator it = qq.begin(); it != qq.end(); it++) {
	  (*it)->addChild(newq);
	  (*it)->addChild(newa);
	  (*it)->addChild(newg);
	}
	ch->currentStep()->addDecayProduct(newq);
	ch->currentStep()->addDecayProduct(newa);
	ch->currentStep()->addDecayProduct(newg);

	// safe 'largest pt so far'.  
	Lorentz5Momentum ptot = newfs[0] + newfs[1] + newfs[2];
	double x = 2*newfs[0]*ptot/sqr(qqg->getQ());
	double xb = 2*newfs[1]*ptot/sqr(qqg->getQ());
	//cerr << x << ", " << xb << ", " << newfs[1].m() << ", " 
	//     << qqg->getM() << endl; 
	Energy qt, pt;
	double z;
	qt = qqg->getQ()*sqrt(qqg->getKfromX(x, xb));
	z = qqg->getZfromX(x, xb);
	pt = (1.-z)*sqrt(sqr(z*qt)-sqr(qqg->getM()));
	_showerVariables->setLargestPtQ(pt);
	qt = qqg->getQ()*sqrt(qqg->getKfromX(xb, x));
	z = qqg->getZfromX(xb, x);
	pt = (1.-z)*sqrt(sqr(z*qt)-sqr(qqg->getM()));
	_showerVariables->setLargestPtQbar(pt);
      }
    } else {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) { 
	generator()->log() << "ShowerHandler::hardMEC: " 
			   << "no ME corr" << endl;
	_showerVariables->setLargestPtQ(Energy());
	_showerVariables->setLargestPtQbar(Energy());
      }
    }
    delete qqg;
  } 
}
