// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InsideRangeShowerEvolver class.
//

#include "InsideRangeShowerEvolver.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "MECorrection.h"
#include "PartnerFinder.h"
#include "ShowerParticle.h"
#include "Pythia7/MatrixElement/MEBase.h" 
#include "Pythia7/MatrixElement/PhaseSpaceBase.h"
#include "Pythia7/Handlers/XComb.h"

using namespace Herwig;


InsideRangeShowerEvolver::~InsideRangeShowerEvolver() {}


void InsideRangeShowerEvolver::persistentOutput(PersistentOStream & os) const {
  os << _pointerPartnerFinder
     << _pointerForwardShowerEvolver
     << _pointerBackwardShowerEvolver
     << _pointerKinematicsReconstructor
     << _pointerSplittingGenerator;
}


void InsideRangeShowerEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _pointerPartnerFinder
     >> _pointerForwardShowerEvolver
     >> _pointerBackwardShowerEvolver
     >> _pointerKinematicsReconstructor
     >> _pointerSplittingGenerator;
}


ClassDescription<InsideRangeShowerEvolver> InsideRangeShowerEvolver::initInsideRangeShowerEvolver;
// Definition of the static class description member.


void InsideRangeShowerEvolver::Init() {

  static ClassDocumentation<InsideRangeShowerEvolver> documentation
    ("This class is responsible for carrying out the showering,",
     "including the kinematics reconstruction, in a given scale range.");

  static Reference<InsideRangeShowerEvolver,PartnerFinder> 
    interfacePartnerFinder("PartnerFinder", 
                           "A reference to the PartnerFinder object", 
                           &Herwig::InsideRangeShowerEvolver::_pointerPartnerFinder,
			   false, false, true, false);
  static Reference<InsideRangeShowerEvolver,ForwardShowerEvolver> 
    interfaceForwardShowerEvolver("ForwardShowerEvolver", 
                                  "A reference to the ForwardShowerEvolver object", 
                                  &Herwig::InsideRangeShowerEvolver::_pointerForwardShowerEvolver,
			          false, false, true, false);
  static Reference<InsideRangeShowerEvolver,BackwardShowerEvolver> 
    interfaceBackwardShowerEvolver("BackwardShowerEvolver", 
                                   "A reference to the BackwardShowerEvolver object", 
                                   &Herwig::InsideRangeShowerEvolver::_pointerBackwardShowerEvolver,
			    false, false, true, false);
  static Reference<InsideRangeShowerEvolver,KinematicsReconstructor> 
    interfaceKinematicsReconstructor("KinematicsReconstructor", 
                                     "A reference to the KinematicsReconstructor object", 
                                     &Herwig::InsideRangeShowerEvolver::_pointerKinematicsReconstructor,
			             false, false, true, false);
  static Reference<InsideRangeShowerEvolver,SplittingGenerator> 
    interfaceSplittingGenerator("SplittingGenerator", 
				"A reference to the SplittingGenerator object", 
				&Herwig::InsideRangeShowerEvolver::_pointerSplittingGenerator,
				false, false, true, false);

}

//--------------------------------------------------------------------------------

void InsideRangeShowerEvolver::setDoneMapShower(MapShower & mapShower) {
  for ( MapShower::iterator it = mapShower.begin(); it != mapShower.end(); ++it ) {
    it->second = false;
  }
}


void InsideRangeShowerEvolver::clear() {
  _mapShowerHardJets.clear();
  _collecMapShowerDecayJets.clear();
}


void InsideRangeShowerEvolver::showerNormally( tPartCollHdlPtr ch, 
					       const tShoConstrPtr showerConstrainer, 
					       const tMECorrectionPtr meCorrection,
					       CollecShoParPtr & particles,
					       bool skipKinReco ) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerNormally "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  // Set the initial evolution scales
  if ( _pointerSplittingGenerator->isInteractionON( ShowerIndex::QCD ) ) {
    _pointerPartnerFinder->setQCDInitialEvolutionScales( showerConstrainer, particles );
  }
  if ( _pointerSplittingGenerator->isInteractionON( ShowerIndex::QED ) ) {
    _pointerPartnerFinder->setQEDInitialEvolutionScales( showerConstrainer, particles );
  }
  if ( _pointerSplittingGenerator->isInteractionON( ShowerIndex::EWK ) ) {
    _pointerPartnerFinder->setEWKInitialEvolutionScales( showerConstrainer, particles );
  }

  // Final State Radiation
  if ( _pointerSplittingGenerator->isFSRadiationON() ) {
    // Remember that is not allowed to add element to a STL container
    // while you are iterating over it. Therefore an additional, temporary,
    // container must be used.
    CollecShoParPtr particlesToShower;
    for ( CollecShoParPtr::const_iterator cit = particles.begin();
	  cit != particles.end(); ++cit ) {
      if ( (*cit)->isFinalState() ) {
	particlesToShower.insert( particlesToShower.end(), *cit );
      }
    }
    for ( CollecShoParPtr::const_iterator cit = particlesToShower.begin();
	  cit != particlesToShower.end(); ++cit ) {
      bool hasEmitted = _pointerForwardShowerEvolver->
	timeLikeShower( ch, showerConstrainer, meCorrection, *cit, particles );
      // Fill  _mapShowerHardJets  with  (*cit, hasEmitted);
      if ( _mapShowerHardJets.find( *cit ) != _mapShowerHardJets.end() ) {
	( _mapShowerHardJets.find( *cit ) )->second = hasEmitted;
      } else {
	_mapShowerHardJets.insert( pair<tShoParPtr, bool>( *cit, hasEmitted ) );
      }
    }
  }

  // Initial State Radiation
  if ( _pointerSplittingGenerator->isISRadiationON() ) {
    CollecShoParPtr particlesToShower;
    for ( CollecShoParPtr::const_iterator cit = particlesToShower.begin();
	  cit != particlesToShower.end(); ++cit ) {
      if ( ! (*cit)->isFinalState() ) {
	particlesToShower.insert( particlesToShower.end(), *cit );
      }
    }
    for ( CollecShoParPtr::iterator it = particlesToShower.begin();
	  it != particlesToShower.end(); ++it ) {
      bool hasEmitted = _pointerBackwardShowerEvolver->
	spaceLikeShower( ch, showerConstrainer, meCorrection, *it, particles );
      if ( _mapShowerHardJets.find( *it ) != _mapShowerHardJets.end() ) {
	( _mapShowerHardJets.find( *it ) )->second = hasEmitted;
      } else {
	_mapShowerHardJets.insert( pair<tShoParPtr, bool>( *it, hasEmitted ) );
      }
    }
  }

  //***LOOKHERE**** global update of rhoD matrices? (maybe or maybe not)

  if ( ! skipKinReco ) {
    reconstructKinematics( ch );
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerNormally "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }
  
}


void InsideRangeShowerEvolver::showerDecay( tPartCollHdlPtr ch, 
					    const tShoConstrPtr showerConstrainer, 
					    const tMECorrectionPtr meCorrection,
					    CollecShoParPtr & particles ) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerDecay "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  //***LOOKHERE*** call one or all of the methods below, in order to fix the
  //                 initial qtilda value, according to which interaction is ON.
  //                 for the decaying particle and its decay: 
  //               _pointerPartnerFinder->setQCDInitialEvolutionScales( showerConstrainer,
  //                                                                    particles, true );
  //               _pointerPartnerFinder->setQEDInitialEvolutionScales( showerConstrainer,
  //                                                                    particles, true );
  //               _pointerPartnerFinder->setEWKInitialEvolutionScales( showerConstrainer,
  //                                                                    particles, true );
  //                 products, to set the qtilda initial values;
  //               MapShower mapDecayJets; 
  //               for each decay products {
  //                 bool hasEmitted = _pointerForwardShowerEvolver->
  //                       timeLikeShower( ch, showerConstrainer, meCorrection,
  //					   particle, particles );
  //                 and add to  mapDecayJets  the element (pointer, hasEmitted);
  //               consider the decaying particle {
  //                 bool hasEmitted = _pointerForwardShowerEvolver->
  //                       timeLikeShower( ch, showerConstrainer, meCorrection,
  //					   particle, particles, true );
  //                 and add to  mapDecayJets  the element (pointer, hasEmitted);
  //               }
  //               global update of rhoD matrices (? maybe not);
  //               _pointerKinematicsReconstructor->reconstructDecayJets( mapDecayJets );
  //               setDoneMapShower(mapDecayJets);  
  //               add  mapDecayJets  to the collection  _collecMapShowerDecayJets ;

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerDecay "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

}


void InsideRangeShowerEvolver::showerGlobally( tPartCollHdlPtr & ch, 
					       const tShoConstrPtr showerConstrainer, 
					       const tMECorrectionPtr meCorrection,
					       CollecShoParPtr & particles,
					       bool skipKinReco )
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerGlobally "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  //***LOOKHERE*** call one or all of the methods below, in order to fix the
  //                 initial qtilda value, according to which interaction is ON,
  //                 where in this case the scale specified in showerConstrainer 
  //                 plays an important role;
  //               _pointerPartnerFinder->setQCDInitialEvolutionScales( showerConstrainer,
  //                                                                    particles );
  //               _pointerPartnerFinder->setQEDInitialEvolutionScales( showerConstrainer,
  //                                                                    particles );
  //               _pointerPartnerFinder->setEWKInitialEvolutionScales( showerConstrainer,
  //                                                                    particles );
  //               for each final state particles { 
  //                 bool hasEmitted = _pointerForwardShowerEvolver->
  //                       timeLikeShower( ch, showerConstrainer, meCorrection,
  //					   particle, particles );
  //                 if ( hasEmitted ) then
  //                   find the parent of the particle, parent, which is the closest
  //                     decaying particle or particle entering the hard subprocess;
  //                   set true the flag of  parent  in  _mapShowerHardJets  if the
  //                     parent is from the hard subprocess, or in 
  //                     _collecMapShowerDecayJets  if the parent is a decaying particle;
  //                 }
  //               }
  //               for each initial state particles { 
  //                 bool hasEmitted = _pointerBackwardShowerEvolver->
  //                       spaceLikeShower( ch, showerConstrainer, meCorrection,
  //				            particle, particles );
  //                 if ( hasEmitted ) then
  //                   find the parent of the particle, parent, which is the incoming
  //                     particle entering the hard subprocess;
  //                   set true the flag of  parent  in  _mapShowerHardJets;
  //                 }
  //               }
  //               global update of rhoD matrices (? maybe not);
  //               if ( ! skipKinReco ) {
  //                  reconstructKinematics( particles );
  //               }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::showerGlobally "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

}


void InsideRangeShowerEvolver::setEffectiveGluonMass( const Energy effectiveGluonMass,
						      const CollecShoParPtr & particles ) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::setEffectiveGluonMass "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  //***LOOKHERE*** 
  //               for each final state particles { 
  //                 if ( particle is a gluon ) {
  //                   set the gluon on the effective mass shell;
  //                   find the parent of the gluon, parent, which is the closest
  //                     decaying particle or particle entering the hard subprocess;
  //                   set true the flag of  parent  in  _mapShowerHardJets  if the
  //                     parent is from the hard subprocess, or in 
  //                     _collecMapShowerDecayJets  if the parent is a decaying particle;
  //                 }
  //               }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "InsideRangeShowerEvolver::setEffectiveGluonMass "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

}

    
void InsideRangeShowerEvolver::reconstructKinematics( tPartCollHdlPtr & ch ) 
  throw (Veto, Stop, Exception) {

  for ( CollecMapShower::iterator it = _collecMapShowerDecayJets.begin();
	it != _collecMapShowerDecayJets.end(); ++it ) {    
    bool recoNeeded = false;
    for ( MapShower::const_iterator citer = (*it).begin();
	    citer != (*it).end(); ++citer ) {
	if ( citer->second ) recoNeeded = true;
    }
    if ( recoNeeded ) {
	_pointerKinematicsReconstructor->reconstructDecayJets( *it );
 	setDoneMapShower( *it );
    }
  }
  PPair beamHadrons = ch->currentCollision()->incoming();
  _pointerKinematicsReconstructor->
    reconstructHardJets( _mapShowerHardJets,
			   beamHadrons.first->momentum(),
			   beamHadrons.second->momentum(),
			   ch->lastME() ); 
  setDoneMapShower(_mapShowerHardJets);

}

