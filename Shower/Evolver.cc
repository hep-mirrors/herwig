// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evolver class.
//

#include "Evolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
//#include "MECorrection.h"
#include "PartnerFinder.h"
#include "ShowerParticle.h"
//#include "ThePEG/MatrixElement/MEBase.h" 
// #include "ThePEG/MatrixElement/PhaseSpaceBase.h"
#include "ThePEG/Handlers/XComb.h"

using namespace Herwig;


Evolver::~Evolver() {}


void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _partnerFinder << _forwardEvolver << _backwardEvolver
     << _kinematicsReconstructor << _splittingGenerator;
}


void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _partnerFinder >> _forwardEvolver >> _backwardEvolver
     >> _kinematicsReconstructor >> _splittingGenerator;
}


ClassDescription<Evolver> Evolver::initEvolver;
// Definition of the static class description member.


void Evolver::Init() {

  static ClassDocumentation<Evolver> documentation
    ("This class is responsible for carrying out the showering,",
     "including the kinematics reconstruction, in a given scale range.");

  static Reference<Evolver,PartnerFinder> 
    interfacePartnerFinder("PartnerFinder", 
                           "A reference to the PartnerFinder object", 
                           &Herwig::Evolver::_partnerFinder,
			   false, false, true, false);
  static Reference<Evolver,ForwardEvolver> 
    interfaceForwardEvolver("ForwardEvolver", 
			    "A reference to the ForwardEvolver object", 
			    &Herwig::Evolver::_forwardEvolver,
			    false, false, true, false);
  static Reference<Evolver,BackwardEvolver> 
    interfaceBackwardEvolver("BackwardEvolver", 
			     "A reference to the BackwardEvolver object", 
			     &Herwig::Evolver::_backwardEvolver,
			     false, false, true, false);
  static Reference<Evolver,KinematicsReconstructor> 
    interfaceKinRecon("KinematicsReconstructor", 
		      "A reference to the KinematicsReconstructor object", 
		      &Herwig::Evolver::_kinematicsReconstructor,
		      false, false, true, false);
  static Reference<Evolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::Evolver::_splittingGenerator,
		      false, false, true, false);
}

//--------------------------------------------------------------------------------

void Evolver::setDoneMapShower(MapShower & mapShower) {
  for(MapShower::iterator it=mapShower.begin(); it!=mapShower.end(); ++it) {
    it->second = false;
  }
}


void Evolver::clear() {
  _mapShowerHardJets.clear();
  _mapShowerDecayJets.clear();
}


void Evolver::showerNormally(tPartCollHdlPtr ch, 
			     const tShowerVarsPtr showerVariables, 
			     //const tMECorrectionPtr meCorrection,
			     ShowerParticleVector & particles,
			     bool skipKinReco) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerNormally(): Evt #" 
		       << generator()->currentEventNumber() 
		       << "\t _______________________"
		       << endl;
  }

  // Set the initial evolution scales
  if(_splittingGenerator->isInteractionON(ShowerIndex::QCD)) {
    _partnerFinder->setQCDInitialEvolutionScales(showerVariables, particles);
  }
  if(_splittingGenerator->isInteractionON(ShowerIndex::QED)) {
    _partnerFinder->setQEDInitialEvolutionScales(showerVariables, particles);
  }
  if(_splittingGenerator->isInteractionON(ShowerIndex::EWK)) {
    _partnerFinder->setEWKInitialEvolutionScales(showerVariables, particles);
  }

  bool reconstructed = false; 
  // catch the possibility of an impossible kinematic reconstruction
  ShowerParticleVector::const_iterator cit;
  while(!reconstructed) {

    // Final State Radiation
    if(_splittingGenerator->isFSRadiationON()) {
      // Remember that is not allowed to add element to a STL container
      // while you are iterating over it. Therefore an additional, temporary,
      // container must be used.
      ShowerParticleVector particlesToShower;
      for(cit = particles.begin(); cit != particles.end(); ++cit) {
	if((*cit)->isFinalState()) particlesToShower.push_back(*cit);
      }
      for(cit = particlesToShower.begin(); 
	  cit != particlesToShower.end(); ++cit) {
	_mapShowerHardJets[*cit] = _forwardEvolver->
		 timeLikeShower(ch, showerVariables, *cit, particles);
      }
    }

    // Initial State Radiation
    if(_splittingGenerator->isISRadiationON()) {
      ShowerParticleVector particlesToShower;
      // This may not be right at all!
      for(cit = particles.begin(); cit != particles.end(); ++cit) {
	if(!(*cit)->isFinalState()) {
	  particlesToShower.push_back(*cit);
	}
      }
      for(ShowerParticleVector::iterator it = particlesToShower.begin();
	    it != particlesToShower.end(); ++it ) {
	_mapShowerHardJets[*it] = _backwardEvolver->
		spaceLikeShower(ch, showerVariables, *it, particles);
      }
    }

    if(!skipKinReco) {
      reconstructed = reconstructKinematics( ch );
    } else { 
      reconstructed = true;     
    }
  }
}


void Evolver::showerDecay(tPartCollHdlPtr ch, 
			  const tShowerVarsPtr showerVariables, 
			  //const tMECorrectionPtr meCorrection,
			  ShowerParticleVector & particles ) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerDecay "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerDecay "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

}


void Evolver::showerGlobally(tPartCollHdlPtr & ch, 
			     const tShowerVarsPtr showerVariables, 
			     //const tMECorrectionPtr meCorrection,
			     ShowerParticleVector & particles,
			     bool skipKinReco)
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerGlobally "
		       << " ===> START DEBUGGING <=== "
		       << "    EventNumber=" 
		       << generator()->currentEventNumber() 
		       << endl;
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerGlobally "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

}


void Evolver::setEffectiveGluonMass( const Energy effectiveGluonMass,
						      const ShowerParticleVector & particles ) 
  throw (Veto, Stop, Exception) {

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "Evolver::setEffectiveGluonMass "
// 		       << " ===> START DEBUGGING <=== "
// 		       << "   EventNumber=" << generator()->currentEventNumber() 
// 		       << endl;
//   }

  for(ShowerParticleVector::const_iterator pit = particles.begin(); 
      pit != particles.end(); ++pit) {   
    if ( (*pit)->data().id() == 21 ) {
      Lorentz5Momentum dum = (*pit)->momentum(); 
      dum.setMass( effectiveGluonMass ); 
      (*pit)->set5Momentum( dum );
    }
  }

  // ***ACHTUNG!*** still flags to be set...

  //***LOOKHERE*** 
  //               for each final state particles { 
  //                 if ( particle is a gluon ) {
  //                   set the gluon on the effective mass shell;
  //                   find the parent of the gluon, parent, which is the closest
  //                     decaying particle or particle entering the hard subprocess;
  //                   set true the flag of  parent  in  _mapShowerHardJets  if the
  //                     parent is from the hard subprocess, or in 
  //                     _mapShowerDecayJets  if the parent is a decaying particle;
  //                 }
  //               }

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "Evolver::setEffectiveGluonMass "
// 		       << " ===> END DEBUGGING <=== " 
// 		       << endl;
//   }

}

    
bool Evolver::reconstructKinematics( tPartCollHdlPtr & ch ) 
  throw (Veto, Stop, Exception) {

  for(MapShowerVector::iterator it = _mapShowerDecayJets.begin();
	it != _mapShowerDecayJets.end(); ++it ) {    
    bool recoNeeded = false;
    for ( MapShower::const_iterator citer = (*it).begin();
	    citer != (*it).end(); ++citer ) {
	if ( citer->second ) recoNeeded = true;
    }
    if ( recoNeeded ) {
	_kinematicsReconstructor->reconstructDecayJets( *it );
 	setDoneMapShower( *it );
    }
  }
  PPair beamHadrons = ch->currentCollision()->incoming();
  bool ok = _kinematicsReconstructor->
    reconstructHardJets( _mapShowerHardJets,
			 beamHadrons.first->momentum(),
			 beamHadrons.second->momentum(),
			 ch->lastME() ); 
  setDoneMapShower(_mapShowerHardJets);
  return ok; 
  
}

