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
//#include "ShowerColourLine.h"
#include "Pythia7/CLHEPWrap/Matrix.h"
#include "Pythia7/PDT/DecayMode.h"
#include "Pythia7/EventRecord/Step.h"

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
    generator()->log() << "ShowerHandler::debuggingInfo() Evt #" << generator()->currentEventNumber() 
		       << "\t _______________________________________" << endl; 
		      
  }

  tPartCollHdlPtr ch = collisionHandler();
  //Particle::PrintParticles(cout, ch->currentStep()->particles().begin(),
  //		   ch->currentStep()->particles().end());
  // const Hint & theHint = hint();	    // OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const PDF & theFirstPDF = firstPDF();  // OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const PDF & theSecondPDF = secondPDF();// OK: COMMENTED TO AVOID COMPILATION WARNINGS
  // const pair<PDF,PDF> & thePdfs = pdfs();// OK: COMMENTED TO AVOID COMPILATION WARNINGS

  // From the Pythia7 particles entering the hard subprocess, create
  // the corresponding starting ShowerParticle objects and put them
  // in the vector hardProcessParticles (not directly in _particles
  // because we could throw away the all showering and restart it).
  // ShowerParticleVector hardProcessParticles;
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
      ShowerParticleVector hardProcessParticles;
      hardProcessParticles.clear();
      createShowerParticlesFromP7Particles( ch, hardProcessParticles );

      // Fill _particles with the particles from the hard subprocess
      for ( ShowerParticleVector::const_iterator cit = hardProcessParticles.begin();
	    cit != hardProcessParticles.end(); ++cit ) {
	if(*cit) { 
	  _particles.push_back( *cit );      
	}
      }

      // clean up children from previous trial showers:
      //      for ( ShowerParticleVector::iterator it = _particles.begin(); 
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
        for ( ShowerParticleVector::const_iterator cit = _particles.begin();
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

        ShowerParticleVector currentFinalOrInitialParticles;
        for ( ShowerParticleVector::const_iterator cit = _particles.begin();
	      cit != _particles.end(); ++cit ) {
	  if ( (*cit)->children().size() == 0 ) {
	    currentFinalOrInitialParticles.insert( currentFinalOrInitialParticles.end() , *cit );
	  }
	}

        for ( ShowerParticleVector::iterator it = currentFinalOrInitialParticles.begin();
	      it != currentFinalOrInitialParticles.end(); ++it ) {
	  
	  Energy particleWidth = (*it)->data().width();
	  if ( particleWidth >= largestWidth ) {

            ShowerParticleVector decayParticles;
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
	    for ( ShowerParticleVector::const_iterator cit = decayParticles.begin();
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
	    ShowerParticleVector::const_iterator cit = decayParticles.begin();
            do {
	      ++cit;
	      _particles.insert( _particles.end(), *cit );
	    } while ( cit != decayParticles.end() );
	    decayParticles.clear();
          }

        } // end for loop

        // Find now the new largest width between all current final state particles.
	Energy newLargestWidth = Energy();
        for ( ShowerParticleVector::const_iterator cit = _particles.begin();
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
      if ( !_pointerGlobalParameters->isPythia7StringFragmentationON() ) {
	_pointerInsideRangeShowerEvolver->setEffectiveGluonMass
	  ( _pointerGlobalParameters->effectiveGluonMass() , _particles );
      } else { 
	_pointerInsideRangeShowerEvolver->setEffectiveGluonMass
	  ( 0 , _particles );
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

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower ) {    
    debuggingInfo();
    if ( generator()->currentEventNumber() < 1000 )
      generator()->log() << "ShowerHandler::debuggingInfo "
			 << " ===> END DEBUGGING <=== "
			 << endl;
  }
  fillEventRecord( ch );  
}

void ShowerHandler::
createShowerParticlesFromP7Particles( const tPartCollHdlPtr ch, 
				      ShowerParticleVector & hardProcessParticles ) {

  // Incoming (initial state) particles (indeed the partons entering
  // the hard subprocess, not the beam hadrons). 
  //Particle::PrintParticles(cout, ch->currentStep()->particles().begin(),
  //			   ch->currentStep()->particles().end());
  ShowerParticlePtr part;
  
  if(ch->lastPartons().first ) {
    part = ptr_new<ShowerParticlePtr>(*ch->lastPartons().first);
    if ( part ) {
      part->setFromHardSubprocess(true);
      part->setP7base(ch->lastPartons().first); 
      part->setFinalState(false);
      hardProcessParticles.push_back(part);
    }
  } 
  if(ch->lastPartons().second ) {
    part = ptr_new<ShowerParticlePtr>(*ch->lastPartons().second);
    if ( part ) {
      part->setFromHardSubprocess(true);
      part->setP7base(ch->lastPartons().second); 
      part->setFinalState(false);
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
      part->setP7base(*cit); 
      hardProcessParticles.push_back(part);
    }
  }    

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "  Total number of initial Shower particles : " 
		       << hardProcessParticles.size() << endl;
      //<< "  Total number of colour lines : " 
      //	       << mapP7toShowerColine.size() << endl;
    for ( ShowerParticleVector::const_iterator cit = hardProcessParticles.begin();
	  cit != hardProcessParticles.end(); ++cit ) {
      generator()->log() << "\t" << (**cit).data().PDGName() 
			 << "  p=" << (**cit).momentum()
			 << ( (**cit).isFinalState() ? ", out" : ", in" )
			 << endl;
    }
  }  
}


void ShowerHandler::fillPositions() {

  // Put on the vector  particlesToExamine  the particles that come
  // directly from the hard subprocess. 
  ParticleVector particlesToExamine;
  for ( ShowerParticleVector::const_iterator cit = _particles.begin();
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
      lambda = particle->momentum().mag2() / _pointerGlobalParameters->minVirtuality2(); 
    }
    LorentzDistance distance( _pointerGlobalParameters->conversionFactorGeVtoMillimeter() 
			      * particle->momentum().vect() / GeV, 
			      _pointerGlobalParameters->conversionFactorGeVtoMillimeter() 
			      * particle->momentum().e() / GeV );
    distance *= log ( 1.0 / rnd() ) / 
      ( sqrt( sqr( particle->momentum().mag2() - particle->momentum().m2() ) +
	      sqr( particle->momentum().mag2() / lambda ) ) / GeV );

    particle->setLabVertex( particle->labVertex() + distance ); // update the position.
         
    for ( ParticleVector::const_iterator cit = particle->children().begin();
	  cit != particle->children().end(); ++cit ) {
      (*cit)->setLabVertex( particle->labVertex()); // initialize the position.
      particlesToExamine.push_back( *cit );
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
		       << ( _pointerGlobalParameters->isPythia7StringFragmentationON() ? 
			    "YES" : "NO" ) 
		       << endl
		       << "  effective gluon mass = " 
		       << _pointerGlobalParameters->effectiveGluonMass() / GeV 
		       << endl
		       << "  hadronization scale = " 
		       << _pointerGlobalParameters->hadronizationScale() / GeV 
		       << "  [GeV] " << endl
		       << "  Info for ShowerConstrainer: " << endl
		       << "  switches:  multi-scale = " 
		       << _pointerShowerConstrainer->isMultiScaleShowerON()
		       << ", decay before = " 
		       << _pointerShowerConstrainer->isDecayBeforeShowerON() 
		       << ", gluino? = " 
		       << ( _pointerShowerConstrainer->hasToDecayBeforeShower(1000021) ? 
			    "Y" : "N" )  
		       << endl
		       << "  cutoff-parameters:  QCD = "
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::QCD) / GeV 
		       << ", QED = "  
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::QED) / GeV 
		       << ", EWK = "  
		       << _pointerShowerConstrainer->cutoffMassScale(ShowerIndex::EWK) / GeV 
		       << " [GeV] " 
		       << endl; 
  }

  if ( generator()->currentEventNumber() < 1000 )
    generator()->log() << "  Shower finished - summary:" << endl;
  
  // stuff that goes to STDOUT comes here:
//   if ( generator()->currentEventNumber() < 1000 )
//     cout << "# event no " << generator()->currentEventNumber() << endl; 

  tShowerParticleVector fs;
  for (ShowerParticleVector::const_iterator cit = _particles.begin(); 
       cit != _particles.end(); ++cit ) {
    if ( (*cit)->isFromHardSubprocess() ) {
      if ( generator()->currentEventNumber() < 1000 ) {
	generator()->log() << "  " << (*cit)->data().PDGName() 
			   << ", Q = " << (*cit)->momentum().mass()/MeV 
			   << ", q = " << (*cit)->momentum()/MeV << " (MeV)"
			   <<endl;
	//	(*cit)->deepPrintInfo(); // writes some showerinfo to STDOUT, obsolete (?)
      }
      //cout << (*cit)->sumParentsMomenta() << endl;      
      if ( (*cit)->isFinalState() ) {
	tShowerParticleVector thisfs = (*cit)->getFSChildren();
	fs.insert(fs.end(), thisfs.begin(), thisfs.end()); 
      }
    }    
  }
  vector<Vector3> n; 
  vector<double> lam; 
  eventShape(fs, lam, n);

  // to play with events in Mathematica, get lists of final state particles...
//   cout << "{";
//   for(tShowerParticleVector::const_iterator cit=fs.begin(); cit != fs.end(); ++cit) {
//     cout << "{" 
// 	 << (*cit)->momentum().vect().x() << ", "
// 	 << (*cit)->momentum().vect().y() << ", "
// 	 << (*cit)->momentum().vect().z() << "}";
//     if (cit != fs.end()) cout << ", " << endl;
//   }
//   cout << "}"; 

//   Lorentz5Momentum pcm = Lorentz5Momentum(); 
//   for(tShowerParticleVector::const_iterator cit=fs.begin(); cit != fs.end(); ++cit) {
//     pcm += (*cit)->momentum();     
//   }
//   Energy root_s = pcm.m();
//   if ( generator()->currentEventNumber() < 1000 )
//     cout << "# shapes: root(s)*lam[i]*n[i] | lam[i] | n[i]" << endl;
//   for (int i=0; i<3; i++) {
//     if ( generator()->currentEventNumber() < 1000 )
//       cout << "666 " << root_s*lam[i]*n[i][0] << " " 
// 	   << root_s*lam[i]*n[i][1] << " " 
// 	   << root_s*lam[i]*n[i][2] << " " 
// 	   << lam[i] << " " << n[i] << endl;
//   }
//   double C_parameter = 3.*(lam[0]*lam[1] + lam[1]*lam[2] + lam[2]*lam[0]);
//   double D_parameter = 27.*(lam[0]*lam[1]*lam[2]);

//   double lam1, lam2, lam3; 
//   double dumd; 
//   lam1 = lam[0]; lam2 = lam[1]; lam3 = lam[2]; 
//   if (lam1 < lam2) {dumd = lam1; lam1 = lam2; lam2 = dumd; }
//   if (lam1 < lam3) {dumd = lam1; lam1 = lam3; lam3 = dumd; }
//   if (lam2 < lam3) {dumd = lam2; lam2 = lam3; lam3 = dumd; }
//   if ( generator()->currentEventNumber() < 1000 )
//     cout << "# shapes C = " << C_parameter
// 	 << ", D = " << D_parameter << endl
// 	 << "# (lam1, lam2, lam3) = (" << lam1 << ", " << lam2 << ", " 
// 	 << lam3 << ")" << endl; 
//   // this one gives 'blocks':
//   if ( generator()->currentEventNumber() < 1000 )
//     cout << endl;

  // book some histograms

/*  HwDebug::lambda1Histo += lam1; 
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
*/
}


void ShowerHandler::fillEventRecord( const tPartCollHdlPtr ch ) {  

  // Transform some of the ShowerParticles object in Pythia7 particles,
  // set properly the parent/child relationships and treat carefully 
  // the transformation from ShowerColourLine objects into Pythia7
  // ColourLine ones; and finally then write them in the Event Record     
  StepPtr pstep;
  pstep = ch->newStep();
  for (ShowerParticleVector::const_iterator cit = _particles.begin(); 
       cit != _particles.end(); ++cit ) {
    if ( (*cit)->isFromHardSubprocess() && (*cit)->isFinalState() ) {
      pstep->addDecayNoCol((*cit)->getP7base(), dynamic_ptr_cast<tPPtr>(*cit));
      // pstep->addDecayProduct((*cit)->getP7base(), dynamic_ptr_cast<tPPtr>(*cit));
      (*cit)->addChildrenEvtRec(pstep);   
    }
  }  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) printStep(pstep, "after filling");
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



// utility method
void ShowerHandler::eventShape(const tShowerParticleVector & p, 
				vector<double> & lam, vector<Vector3> & n) {

  // get cm-frame
  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  for(tShowerParticleVector::const_iterator cit=p.begin(); cit != p.end(); ++cit) {
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
  for(tShowerParticleVector::const_iterator cit=p.begin(); cit != p.end(); ++cit) {
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
