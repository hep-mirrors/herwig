// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Pythia7/Handlers/CollisionHandler.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Handlers/Hint.h"
#include "Pythia7/PDF/PDF.h"
#include "ShowerConfig.h"
#include "ShowerIndex.h"
#include "Pythia7/MatrixElement/MEBase.h"
#include "Pythia7/PDT/Decayer.h"
#include "Pythia7/PDT/ParticleData.h"
#include "Pythia7/MatrixElement/PhaseSpaceBase.h"
#include "Pythia7/Handlers/XComb.h"
#include "ShowerParticle.h"
#include "Pythia7/EventRecord/ColourLine.h"
#include "ShowerColourLine.h"
#include "Pythia7/CLHEPWrap/Matrix.h"

using namespace Herwig;


ShowerHandler::~ShowerHandler() {}


void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _pointerGlobalParameters 
     << _pointerMECorrections 
     << _pointerShowerConstrainer
     << _pointerInsideRangeShowerEvolver;
}


void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _pointerGlobalParameters 
     >> _pointerMECorrections 
     >> _pointerShowerConstrainer
     >> _pointerInsideRangeShowerEvolver;
}


ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.


void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.");

  static Reference<ShowerHandler,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
			      "A reference to the GlobalParameters object", 
			      &Herwig::ShowerHandler::_pointerGlobalParameters,
			      false, false, true, false);
  static Reference<ShowerHandler,MECorrections> 
    interfaceMECorrections("MECorrections", 
			   "A reference to the MECorrections object", 
			   &Herwig::ShowerHandler::_pointerMECorrections,
			   false, false, true, false);
  static Reference<ShowerHandler,ShowerConstrainer> 
    interfaceShowerConstrainer("ShowerConstrainer", 
			       "A reference to the ShowerConstrainer object", 
			       &Herwig::ShowerHandler::_pointerShowerConstrainer,
			       false, false, true, false);
  static Reference<ShowerHandler,InsideRangeShowerEvolver> 
    interfaceInsideRangeShowerEvolver("InsideRangeShowerEvolver", 
				      "A reference to the InsideRangeShowerEvolver object", 
				      &Herwig::ShowerHandler::_pointerInsideRangeShowerEvolver,
				      false, false, true, false);

}


void ShowerHandler::cascade() {

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "ShowerHandler::debuggingInfo "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() << endl;
  }

  tPartCollHdlPtr ch = collisionHandler();

  // const Hint & theHint = hint();	    // OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const PDF & theFirstPDF = firstPDF();  // OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const PDF & theSecondPDF = secondPDF();// OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const pair<PDF,PDF> & thePdfs = pdfs();// OK: COMMENTED TO AVOID COMPILATION WARNINGS

  // From the Pythia7 particles entering the hard subprocess, create
  // the corresponding starting ShowerParticle objects and put them
  // in the vector hardProcessParticles (not directly in _particles
  // because we could throw away the all showering and restart it).
  // CollecShoParPtr hardProcessParticles;
  //  createShowerParticlesFromP7Particles( ch, hardProcessParticles );

  const int maxNumFailures = 100;
  int countFailures = 0;
  bool keepTrying = true;
  
  while ( keepTrying  &&  countFailures < maxNumFailures ) {       
    try {
      // Cleaning and initializing
      _particles.clear();
      _pointerShowerConstrainer->reset();     
      _pointerInsideRangeShowerEvolver->clear();

      // even more cleaning... ***ACHTUNG*** find something much, much
      // more sophisticated to save old properties. In particular for
      // the multiscaleshower...
      CollecShoParPtr hardProcessParticles;
      hardProcessParticles.clear();
      createShowerParticlesFromP7Particles( ch, hardProcessParticles );

      // Fill _particles with the particles from the hard subprocess
      for ( CollecShoParPtr::const_iterator cit = hardProcessParticles.begin();
	    cit != hardProcessParticles.end(); ++cit ) {
	_particles.push_back( *cit );      
      }

      // clean up children from previous trial showers:
      //      for ( CollecShoParPtr::iterator it = _particles.begin(); 
      // 	    it != _particles.end(); ++it) {
      //	(*it)->removeChildren(); 
      //}
      
      bool isMECorrectionApplied = false;
      tMECorrectionPtr softMECorrection = tMECorrectionPtr();
      if ( _pointerMECorrections->isMECorrectionsON() &&
	   _pointerMECorrections->getMECorrection( ch->lastME() ) ) { 
        isMECorrectionApplied = true;
	if ( _pointerMECorrections->getMECorrection( ch->lastME() )->hardProcessPlusJetME() ) {
	  _pointerMECorrections->getMECorrection( ch->lastME() )->hardMECorrection();
	} else { 
	  softMECorrection = _pointerMECorrections->getMECorrection( ch->lastME() );
	}
      }

      if ( _pointerShowerConstrainer->isDecayBeforeShowerON() ) {      
        // Loop over _particles and for those that satisfy the condition
        // if ( _pointerShowerConstrainer->hasToDecayBeforeShower( id ) ) 
        // then decay them. Apply the same procedure iteratively for
        // all decay products, until none particles in _particles
        // satisfies the above condition.
        // ***LOOKHERE*** --- code here --- 
      }

      Energy largestWidth = Energy();
      if ( _pointerShowerConstrainer->isMultiScaleShowerON() ) {      
        for ( CollecShoParPtr::const_iterator cit = _particles.begin();
	      cit != _particles.end(); ++cit ) {
	  if ( (*cit)->children().size() ) {
	    if ( (*cit)->data().width() > largestWidth ) { 
	      largestWidth = (*cit)->data().width();
	    }
	  }
	}
	if ( largestWidth < _pointerGlobalParameters->hadronizationScale() ) {  
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
      if ( largestWidth < _pointerGlobalParameters->hadronizationScale() ) {  
	skipKinReco = true;
      }      

      _pointerShowerConstrainer->stopShowerAtMassScale( largestWidth ); 
      _pointerShowerConstrainer->vetoBelowPtScale( largestWidth );  //***MAYBE NOT NEEDED***
      
      _pointerInsideRangeShowerEvolver->showerNormally
	( ch, _pointerShowerConstrainer, softMECorrection, _particles, skipKinReco );

      // The following loop is done only for multi-scale showering.
      while ( largestWidth > _pointerGlobalParameters->hadronizationScale() ) {  

        // Because is forbidden to add elements to a STL container while
	// looping over it, we have to create first a new container which
        // holds the particles that need to be considered, that is the 
        // current initial and final states; then we can loop over the
        // new container, and adding to the initial one the new particles
        // that are produced by the showering.

        CollecShoParPtr currentFinalOrInitialParticles;
        for ( CollecShoParPtr::const_iterator cit = _particles.begin();
	      cit != _particles.end(); ++cit ) {
	  if ( (*cit)->children().size() == 0 ) {
	    currentFinalOrInitialParticles.insert( currentFinalOrInitialParticles.end() , *cit );
	  }
	}

        for ( CollecShoParPtr::iterator it = currentFinalOrInitialParticles.begin();
	      it != currentFinalOrInitialParticles.end(); ++it ) {
	  
	  Energy particleWidth = (*it)->data().width();
	  if ( particleWidth >= largestWidth ) {

            CollecShoParPtr decayParticles;
	    // Decay the particle *it (need pointer to Decayer)
            // and put the decaying ShowerParticle object and all its
            // decay products (for each of them we need to create a
            // ShowerParticle object) in the collection decayParticles.
            // ***LOOKHERE*** --- code here ---

	    tMECorrectionPtr softMEDecayCorrection = tMECorrectionPtr();
	    if ( _pointerMECorrections->isMECorrectionsON()  &&
		 ( ( ! isMECorrectionApplied ) || 
		   ( isMECorrectionApplied  &&  _pointerMECorrections->isComposeMECorrectionsON() )
		   ) ) {
	      // if ( _pointerMECorrections->getMECorrection( "decayer pointer" ) {
	      //    isMECorrectionApplied = true;
	      //    if ( _pointerMECorrections->
              //         getMECorrection( "decayer pointer" )->decayProcessPlusJetME() ) {
	      //      _pointerMECorrections->getMECorrection( ch->lastME() )->hardMECorrection();
	      //    } else { 
	      //      softMEDecayCorrection = _pointerMECorrections->
              //                              getMECorrection( "decayer pointer" );
	      //    }
              // }
	    }

            // In order to proper treat exceptional cases in which 
            // a decay product has larger width than the decaying parent, 
            // we have to keep track of the largest width for each
            // decaying system (including the decaying parent), and
            // also of the global largest width for all possible
            // decaying systems.  
            Energy largestWidthDecayingSystem = Energy();
	    for ( CollecShoParPtr::const_iterator cit = decayParticles.begin();
		  cit != decayParticles.end(); ++cit ) {
	      if ( (*cit)->data().width() > largestWidthDecayingSystem ) { 
		largestWidthDecayingSystem = (*cit)->data().width();
	      }
	    }
	    if ( largestWidthDecayingSystem > largestWidth ) {
	      largestWidth = largestWidthDecayingSystem;
	    }

	    _pointerShowerConstrainer->stopShowerAtMassScale( largestWidthDecayingSystem ); 
	    _pointerShowerConstrainer->vetoBelowPtScale( largestWidthDecayingSystem );  //***MAYBE NOT NEEDED***

	    _pointerInsideRangeShowerEvolver->showerDecay
	      ( ch, _pointerShowerConstrainer, softMEDecayCorrection, decayParticles );
            
	    // Add the particles in  decayParticles  in  _particles .
            // Notice that the first particle in  decayParticles  must be
            // skipped because it corresponds to the decaying particle
            // which was already in  _particles .
	    CollecShoParPtr::const_iterator cit = decayParticles.begin();
            do {
	      ++cit;
	      _particles.insert( _particles.end(), *cit );
	    } while ( cit != decayParticles.end() );
	    decayParticles.clear();
          }

        } // end for loop

        // Find now the new largest width between all current final state particles.
	Energy newLargestWidth = Energy();
        for ( CollecShoParPtr::const_iterator cit = _particles.begin();
	      cit != _particles.end(); ++cit ) {
	  if ( (*cit)->children().size() ) {
	    if ( (*cit)->data().width() > newLargestWidth ) { 
	      newLargestWidth = (*cit)->data().width();
	    }
	  }
	}
	if ( newLargestWidth < _pointerGlobalParameters->hadronizationScale() ) {  
	  newLargestWidth = Energy();
	  skipKinReco = true;          // to avoid to do twice the kinematics reco
	}                              

	Energy savedVetoAbovePtScale = _pointerShowerConstrainer->vetoAbovePtScale();
	
	_pointerShowerConstrainer->stopShowerAtMassScale( newLargestWidth );
	_pointerShowerConstrainer->vetoBelowPtScale( newLargestWidth ); //***MAYBE NOT NEEDED***
	_pointerShowerConstrainer->vetoAbovePtScale( largestWidth );
	
	_pointerInsideRangeShowerEvolver->showerGlobally
	  ( ch, _pointerShowerConstrainer, softMECorrection, _particles, skipKinReco );
	
	_pointerShowerConstrainer->vetoBelowPtScale( savedVetoAbovePtScale );
	largestWidth = newLargestWidth;

      } // end of while loop

      // In the case Herwig++ Cluster Hadronization model is used for the
      // hadronization, set all final state gluons on the effective mass
      // shell (rather on the physical massless shell).
      if ( _pointerGlobalParameters->isPythia7StringFragmentationON() ) {
	_pointerInsideRangeShowerEvolver->setEffectiveGluonMass
	  ( _pointerGlobalParameters->effectiveGluonMass() , _particles );
      }

      // Do the final kinematics reconstruction
      bool recons = _pointerInsideRangeShowerEvolver->reconstructKinematics( ch );

      // Fill the positions information for all the ShowerParticle objects
      // in _particles. It is done at this stage in order to avoid the
      // complications of boosting such positions during the kinematics
      // reshuffling.
      fillPositions();

      keepTrying = false;
      if (!recons) keepTrying = true; 
    } // end of try
    catch ( std::exception & e ) {
      countFailures++;
    }    
  } // end main while loop

  fillEvenRecord( ch );  

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower ) {    
    debuggingInfo();
    if ( generator()->currentEventNumber() < 1000 )
      generator()->log() << "ShowerHandler::debuggingInfo "
			 << " ===> END DEBUGGING <=== "
			 << endl;
  }

}


void ShowerHandler::
createShowerParticlesFromP7Particles( const tPartCollHdlPtr ch, 
				      CollecShoParPtr & hardProcessParticles ) {

  // We need a map which keeps track of the correspondence between
  // the Pythia7 particles which enter the hard subprocess and
  // the ShowerParticles objects we create from these Pythia7 particles.
  // Similarly, in order to properly reproduce the colour connections. 
  // we also need a map of pointers to all Pythia7 ColourLine objects 
  // which are referenced by these Pythia7 particles, and the corresponding
  // ShowerColourLine objects that we will create later on.
  map<tPPtr,tShoParPtr> mapP7toShowerParticle;
  map<tColinePtr,ShoColinePtr> mapP7toShowerColine;

  // Incoming (initial state) particles (indeed the partons entering
  // the hard subprocess, not the beam hadrons). 
  for ( int i=1; i<=2; ++i ) {
    tPPtr p7part = tPPtr();
    if ( i == 1 && ch->lastPartons().first ) {
      p7part = ch->lastPartons().first;
    } else if ( i == 2 && ch->lastPartons().second ) {
      p7part = ch->lastPartons().second;      
    }
    if ( p7part ) {      
      hardProcessParticles.push_back( new_ptr( ShowerParticle( *p7part ) ) );    
      tShoParPtr shopart = hardProcessParticles.back();
      shopart->isFinalState( false );
      mapP7toShowerParticle.insert( pair<tPPtr,tShoParPtr>(p7part,shopart) );
      if ( p7part->colourLine() &&
	   mapP7toShowerColine.find( p7part->colourLine() ) == mapP7toShowerColine.end() ) {
	mapP7toShowerColine.insert( pair<tColinePtr,ShoColinePtr>( p7part->colourLine(), 
								   ShoColinePtr() ) );
      }
      if ( p7part->antiColourLine() &&
	   mapP7toShowerColine.find( p7part->antiColourLine() ) == mapP7toShowerColine.end() ) {
	mapP7toShowerColine.insert( pair<tColinePtr,ShoColinePtr>( p7part->antiColourLine(), 
								   ShoColinePtr() ) );
      }
    }
  }

  // Outgoing (final state) particles, excluding the beam remnants. 
  // Notice that we don't need to set true the flag ShowerParticle::isFinalState 
  // because it is set true by default. 
  tParticleSet remnantSet = ch->currentCollision()->getRemnants();
  for ( ParticleSet::const_iterator cit = ch->currentStep()->particles().begin();
	cit != ch->currentStep()->particles().end(); ++cit ) {
    if ( remnantSet.find( (*cit)->original() ) == remnantSet.end() ) {
      hardProcessParticles.push_back( new_ptr( ShowerParticle( **cit ) ) );    
      mapP7toShowerParticle.insert( pair<tPPtr,tShoParPtr>( *cit, hardProcessParticles.back() ) );
      if ( (*cit)->colourLine() &&
	   mapP7toShowerColine.find( (*cit)->colourLine() ) == mapP7toShowerColine.end() ) {
	mapP7toShowerColine.insert( pair<tColinePtr,ShoColinePtr>( (*cit)->colourLine(), 
								   ShoColinePtr() ) );
      }
      if ( (*cit)->antiColourLine() &&
	   mapP7toShowerColine.find( (*cit)->antiColourLine() ) == mapP7toShowerColine.end() ) {
	mapP7toShowerColine.insert( pair<tColinePtr,ShoColinePtr>( (*cit)->antiColourLine(), 
								    ShoColinePtr() ) );
      }
    }
  }

  // Now that we have the map:  Pythia7 particle  ===>  ShowerParticle object
  // we can complete the other map:
  //        Pythia7 ColourLine object ===>  ShowerColourLine object
  // which, at the moment, only has the first entry (key), whereas
  // the second part (value of the map) is null. We want to create a
  // ShowerColourLine object for each Pythia7 ColourLine object which
  // is actually used for the colour connection between the particles
  // involved in the hard subprocess, but not in the case that such
  // ColourLine object is used to connect with the a beam remnant. 
  for ( map<tColinePtr,ShoColinePtr>::iterator colineIt = mapP7toShowerColine.begin();
	colineIt != mapP7toShowerColine.end(); ++colineIt ) {

    // Loop over  coloured()  Pythia7 particles connected to this coline.
    // If the Pythia7 particle is in the map (therefore it is not a remnant)
    // then create a new ShowerColourLine object if the considered Pythia7
    // ColourLine object has not yet a corresponding ShowerColourLine object,
    // and then add to it, as coloured ShowerParticle object, the one
    // corresponding to the above Pythia7 particle, if this was not already
    // done before.
    for ( tPVector::const_iterator p7partIt = colineIt->first->coloured().begin();
	  p7partIt != colineIt->first->coloured().end(); ++p7partIt ) {
      if ( mapP7toShowerParticle.find( *p7partIt ) != mapP7toShowerParticle.end() ) {
	if ( ! colineIt->second ) {
	  colineIt->second = new_ptr( ShowerColourLine() );
	} 
	if ( ! mapP7toShowerParticle.find( *p7partIt )->second->colourLine() ) { 
	  colineIt->second->addColoured( mapP7toShowerParticle.find( *p7partIt )->second );
	}
      }
    }

    // As above, but considering now anti-coloured particles. 
    for ( tPVector::const_iterator p7partIt = colineIt->first->antiColoured().begin();
	  p7partIt != colineIt->first->antiColoured().end(); ++p7partIt ) {
      if ( mapP7toShowerParticle.find( *p7partIt ) != mapP7toShowerParticle.end() ) {
	if ( ! colineIt->second ) {
	  colineIt->second = new_ptr( ShowerColourLine() );
	} 
	if ( ! mapP7toShowerParticle.find( *p7partIt )->second->antiColourLine() ) { 
	  colineIt->second->addAntiColoured( mapP7toShowerParticle.find( *p7partIt )->second );
	}
      }
    }

  }  

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "  Total number of initial Shower particles : " 
		       << hardProcessParticles.size() << endl
                       << "  Total number of colour lines : " 
		       << mapP7toShowerColine.size() << endl;
    for ( CollecShoParPtr::const_iterator cit = hardProcessParticles.begin();
	  cit != hardProcessParticles.end(); ++cit ) {
      generator()->log() << "\t" << (**cit).data().PDGName() 
			 << "  p=" << (**cit).momentum()
			 << ( (**cit).isFinalState() ? "   OUTGOING" : "   INCOMING" )
			 << endl;
    }
  }  
}


void ShowerHandler::fillPositions() {

  // Put on the vector  particlesToExamine  the particles that come
  // directly from the hard subprocess. 
  vector<tShoParPtr> particlesToExamine;
  for ( CollecShoParPtr::const_iterator cit = _particles.begin();
	cit != _particles.end(); ++cit ) {
    if ( (*cit)->isFromHardSubprocess() ) {
      particlesToExamine.push_back( *cit );
    }
  }

  // Loop over the vector  particlesToExamine , consider its last element
  // and removing it from the vector, then update the position of such
  // element (previous position + calculated displacement), and finally 
  // insert in the vector  particlesToExamine  its children and set their
  // initial positions to the current parent position,
  // until the vector  particlesToExamine  is empty.
  while ( ! particlesToExamine.empty() ) {

    tShoParPtr particle = particlesToExamine.back();
    particlesToExamine.pop_back();

    // See Mark Smith's thesis, page 121, formula (4.106), for the  
    // calculation of the displacement of the particle.
    // Notice that the displacement is already expressed in the Lab frame,
    // therefore no boost is necessary.
    double lambda;
    if ( particle->decayer() ) { // unstable, decayed particle
      lambda = particle->momentum().m() / particle->data().width();
    } else {                     // not decayed particle
      lambda = particle->momentum().mag2() / _pointerGlobalParameters->minVirtuality2(); 
    }
    LorentzDistance distance( _pointerGlobalParameters->conversionFactorGeVtoMillimeter() 
			      * particle->momentum().vect() / GeV, 
			      _pointerGlobalParameters->conversionFactorGeVtoMillimeter() 
			      * particle->momentum().e() / GeV );
    distance *= log ( 1.0 / rnd() ) / 
      ( sqrt( sqr( particle->momentum().mag2() - particle->momentum().m2() ) +
	      sqr( particle->momentum().mag2() / lambda ) ) / GeV );

    particle->position( particle->position() + distance ); // update the position.
         
    for ( CollecShoParPtr::const_iterator cit = particle->children().begin();
	  cit != particle->children().end(); ++cit ) {
      (*cit)->position( particle->position() ); // initialize the position.
      particlesToExamine.push_back( *cit );
    }

  }

}


void ShowerHandler::debuggingInfo() {

  if ( generator()->currentEventNumber() == 1 ) {
    // Preliminary information to print only once.
    generator()->log() << " --- Info for GlobalParameters  --- " << endl
		       << "\t isOnStringFrag = " 
		       << ( _pointerGlobalParameters->isPythia7StringFragmentationON() ? 
			    "YES" : "NO" )
		       << "   effective gluon mass = " 
		       << _pointerGlobalParameters->effectiveGluonMass() / GeV 
		       << "   hadronization scale = " 
		       << _pointerGlobalParameters->hadronizationScale() / GeV 
		       << "  [GeV] " << endl
		       << " --- Info for ShowerConstrainer --- " << endl
		       << "\t SWITCHES :  multi-scale = " 
		       << _pointerShowerConstrainer->isMultiScaleShowerON()
		       << "   decay before = " 
		       << _pointerShowerConstrainer->isDecayBeforeShowerON() 
		       << "   gluino? = " 
		       << ( _pointerShowerConstrainer->hasToDecayBeforeShower(1000021) ? 
			    " YES " : " NO " )  
		       << endl
		       << "\t PARAMETERS : cutoff :  QCD = "
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::QCD) / GeV 
		       << "   QED = "  
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::QED) / GeV 
		       << "   EWK = "  
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::EWK) / GeV 
		       << "   [GeV] " 
		       << endl; 
  }

  if ( generator()->currentEventNumber() < 1000 )
    generator()->log() << "--- Shower finished - summary ---" << endl;

  
  // stuff that goes to STDOUT comes here:
  if ( generator()->currentEventNumber() < 1000 )
    cout << "# event no " << generator()->currentEventNumber() << endl; 

  tCollecShoParPtr fs;
  for (CollecShoParPtr::const_iterator cit = _particles.begin(); 
       cit != _particles.end(); ++cit ) {
    if ( (*cit)->isFromHardSubprocess() ) {
      if ( generator()->currentEventNumber() < 1000 ) {
	generator()->log() << "  " << (*cit)->data().PDGName() 
			   << ", Q = " << (*cit)->momentum().mass()/MeV 
			   << ", q = " << (*cit)->momentum()/MeV << " (MeV)"
			   <<endl;
	(*cit)->deepPrintInfo();
      }
      //cout << (*cit)->sumParentsMomenta() << endl;      
      if ( (*cit)->isFinalState() ) {
	tCollecShoParPtr thisfs = (*cit)->getFSChildren();
	fs.insert(fs.end(), thisfs.begin(), thisfs.end()); 
      }
    }    
  }
  vector<Vector3> n; 
  vector<double> lam; 
  eventShape(fs, lam, n);

  // to play with events in Mathematica, get lists of final state particles...
//   cout << "{";
//   for(tCollecShoParPtr::const_iterator cit=fs.begin(); cit != fs.end(); ++cit) {
//     cout << "{" 
// 	 << (*cit)->momentum().vect().x() << ", "
// 	 << (*cit)->momentum().vect().y() << ", "
// 	 << (*cit)->momentum().vect().z() << "}";
//     if (cit != fs.end()) cout << ", " << endl;
//   }
//   cout << "}"; 

  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  for(tCollecShoParPtr::const_iterator cit=fs.begin(); cit != fs.end(); ++cit) {
    pcm += (*cit)->momentum();     
  }
  Energy root_s = pcm.m();
  if ( generator()->currentEventNumber() < 1000 )
    cout << "# shapes: root(s)*lam[i]*n[i] | lam[i] | n[i]" << endl;
  for (int i=0; i<3; i++) {
    if ( generator()->currentEventNumber() < 1000 )
      cout << "666 " << root_s*lam[i]*n[i][0] << " " 
	   << root_s*lam[i]*n[i][1] << " " 
	   << root_s*lam[i]*n[i][2] << " " 
	   << lam[i] << " " << n[i] << endl;
  }
  double C_parameter = 3.*(lam[0]*lam[1] + lam[1]*lam[2] + lam[2]*lam[0]);
  double D_parameter = 27.*(lam[0]*lam[1]*lam[2]);

  double lam1, lam2, lam3; 
  double dumd; 
  lam1 = lam[0]; lam2 = lam[1]; lam3 = lam[2]; 
  if (lam1 < lam2) {dumd = lam1; lam1 = lam2; lam2 = dumd; }
  if (lam1 < lam3) {dumd = lam1; lam1 = lam3; lam3 = dumd; }
  if (lam2 < lam3) {dumd = lam2; lam2 = lam3; lam3 = dumd; }
  if ( generator()->currentEventNumber() < 1000 )
    cout << "# shapes C = " << C_parameter
	 << ", D = " << D_parameter << endl
	 << "# (lam1, lam2, lam3) = (" << lam1 << ", " << lam2 << ", " 
	 << lam3 << ")" << endl; 
  // this one gives 'blocks':
  if ( generator()->currentEventNumber() < 1000 )
    cout << endl;

  // book some histograms

  HwDebug::lambda1Histo += lam1; 
  HwDebug::lambda2Histo += lam2;
  HwDebug::lambda3Histo += lam3; 
  HwDebug::CparameterHisto += C_parameter; 
  HwDebug::DparameterHisto += D_parameter; 
  HwDebug::multiplicityHisto += fs.size();

  if ( generator()->currentEventNumber() < 1000 || 
       (generator()->currentEventNumber() % 1000) == 0 ) {
    HwDebug::lambda1Histo.printGnuplot("plots/lambda1.dat");
    HwDebug::lambda2Histo.printGnuplot("plots/lambda2.dat");
    HwDebug::lambda3Histo.printGnuplot("plots/lambda3.dat"); 
    HwDebug::CparameterHisto.printGnuplot("plots/Cpara.dat"); 
    HwDebug::DparameterHisto.printGnuplot("plots/Dpara.dat"); 
    HwDebug::multiplicityHisto.printGnuplot("plots/multiplicity.dat");
  }

}


void ShowerHandler::fillEvenRecord( const tPartCollHdlPtr ch ) {  

  // Transform some of the ShowerParticles object in Pythia7 particles,
  // set properly the parent/child relationships and treat carefully 
  // the transformation from ShowerColourLine objects into Pythia7
  // ColourLine ones; and finally then write them in the Event Record     
  //   StepPtr pstep = ch.newStep();
  // ***LOOKHERE*** --- code here ---

}


  // utility method
void ShowerHandler::eventShape(const tCollecShoParPtr & p, 
				vector<double> & lam, vector<Vector3> & n) {

  // get cm-frame
  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  for(tCollecShoParPtr::const_iterator cit=p.begin(); cit != p.end(); ++cit) {
    pcm += (*cit)->momentum();     
  }
  Vector3 beta = pcm.findBoostToCM(); 
  HepSymMatrix Theta = HepSymMatrix(3);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Theta[i][j] = 0.0;
    }
  }
  //  cout << "# beta = " << beta << endl; 
  double sum = 0.; 
  Vector3 sumvec = Vector3();
  
  // get Theta_ij
  for(tCollecShoParPtr::const_iterator cit=p.begin(); cit != p.end(); ++cit) {
    Lorentz5Momentum dum = (*cit)->momentum();
    dum.boost( beta );
    Vector3 pvec = dum.vect();
    sumvec += pvec;
    sum += pvec.mag();
//     cout << "# pvec = " << pvec << endl;
//     cout << "# pvec[i] = (" << pvec[0] << ", " << pvec[1] << ", " << pvec[2] << ")" << endl; 
//     cout << "# pvec.mag() = " << pvec.mag() << endl; 
//     cout << "# sqrt(pvec^2) = " << sqrt(sqr(pvec[0]) + sqr(pvec[1]) + sqr(pvec[2])) << endl; 
    for(int i=0; i<3; i++) {
      for(int j=i; j<3; j++) {
	Theta[i][j] += (pvec[i])*(pvec[j])/(pvec.mag());
      }
    }
    //   cout << "# stepwise tr/sum = " << Theta.trace()/sum << endl; 
  }

  Theta /= sum;  
    
//   cout << "#-- Theta = " << endl; 
//   for(int i=0; i<3; i++) {
//     for(int j=0; j<3; j++) {
//       cout << Theta[i][j] << "\t";
//     }
//     cout << endl; 
//   }    

//  cout << "# pcm = " << pcm << endl
//       << "# sumvec = " << sumvec << endl;
  // diagonalize it
  HepMatrix U = diagonalize(&Theta);

//   cout << "#-- Theta diagonalized = " << endl; 
//   for(int i=0; i<3; i++) {
//     for(int j=0; j<3; j++) {
//       cout << Theta[i][j] << "\t";
//     }
//     cout << endl; 
//   }    
//   cout << "#-- U = " << endl; 
//   for(int i=0; i<3; i++) {
//     for(int j=0; j<3; j++) {
//       cout << U[i][j] << "\t";
//     }
//     cout << endl; 
//   }  

  for(int i=0; i<3; i++) {
    lam.push_back( Theta[i][i] );
    Vector3 ndum;
    for(int j=0; j<3; j++) {
      ndum[j] = U[j][i]; 
    }
    n.push_back( ndum ); 
  }
  // checks
//   double lamsum = 0.;
//   for(int i=0; i<3; i++) {
//     lamsum += lam[i]; 
//     cout << "#(" << i+1 << ") lam = " 
// 	 << lam[i] << ", n = " << n[i] << endl; 
//   }
//  HepVector n1, n2, n3;
//   n1 = n[0]; 
//   n2 = n[1]; 
//   n3 = n[2];
//   cout << "# n1 = " << n1 << endl 
//        << "# n2 = " << n2 << endl 
//        << "# n3 = " << n3 << endl ;
//   cout << "# check othonormality of n's:" << endl 
//        << dot(n1, n1) << "\t" << dot(n1, n2) << "\t" << dot(n1, n3) << endl
//        << dot(n2, n1) << "\t" << dot(n2, n2) << "\t" << dot(n2, n3) << endl
//        << dot(n3, n1) << "\t" << dot(n3, n2) << "\t" << dot(n3, n3) << endl;
//   cout << "#    lamsum = " << lamsum << endl; 
  
}
