// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evolver class.
//

#include "Evolver.h"
#include "Herwig++/Hadronization/Remnant.h"
#include "Herwig++/Utilities/EnumParticles.h"
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

//-----------------------------------------------------------------------------

void Evolver::setDoneMapShower(MapShower & mapShower) {
  for(MapShower::iterator it=mapShower.begin(); it!=mapShower.end(); ++it) {
    it->second = false;
  }
}


void Evolver::clear() {
  _mapShowerHardJets.clear();
  _mapShowerDecayJets.clear();
}


bool Evolver::showerNormally(tEHPtr ch, 
			     const tShowerVarsPtr showerVariables, 
			     //const tMECorrectionPtr meCorrection,
			     ShowerParticleVector & particles,
			     bool skipKinReco) 
  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Evolver::showerNormally(): Evt #" 
      		       << generator()->currentEventNumber() 
		       << " ________________________________________" 
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
  bool showerOK = true;
  // catch the possibility of an impossible kinematic reconstruction
  ShowerParticleVector::const_iterator cit;
  while(!reconstructed && showerOK) {
    // Initial State Radiation
    if(_splittingGenerator->isISRadiationON()) {
      ShowerParticleVector particlesToShower;
      // This may not be right at all!
      for(cit = particles.begin(); cit != particles.end(); ++cit) {
	if(!(*cit)->isFinalState()) particlesToShower.push_back(*cit);
      }
      for(ShowerParticleVector::iterator it = particlesToShower.begin();
	    it != particlesToShower.end(); ++it ) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log()  << "  calling backwardsEvolver with " 
			      << *it << endl;
	}
	int status =
	  _backwardEvolver->spaceLikeShower(ch, showerVariables, 
					    *it, particles);
	if (status == 1) _mapShowerHardJets[*it] = true;
	else if (status == 0) _mapShowerHardJets[*it] = false;
	else if (status == -1) showerOK = false;
      }
      reconstructed = reconstructISKinematics(ch);
    } else reconstructed = true;
  }

  if (!showerOK) {
    //    cerr << "Evolver::showerNormally: !showerOK from spaceLikeShower, will return false." << endl;
    return false;
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Successful reconstruction of IS shower!" << endl;
    for(cit = particles.begin(); cit != particles.end(); ++cit) {
      if ((*cit)->isFromHardSubprocess()) 
	generator()->log() << (*cit)->id() << " " 
			   << (*cit)->momentum() << endl;
    }
  }

  makeRemnants(particles);
  //  cerr << "before setQCDInitialEvolutionScales." << endl;
  //   if(_splittingGenerator->isInteractionON(ShowerIndex::QCD)) {
  //     _partnerFinder->setQCDInitialEvolutionScales(showerVariables, particles);
  //   }

  reconstructed = false;
  while(!reconstructed) {
    // Final State Radiation
    if(_splittingGenerator->isFSRadiationON()) {
      // Remember that is not allowed to add element to a STL Vector
      // while you are iterating over it. Therefore an additional, temporary,
      // container must be used.
      ShowerParticleVector particlesToShower;
      for(cit = particles.begin(); cit != particles.end(); ++cit) {
	if((*cit)->isFinalState()) particlesToShower.push_back(*cit);

      }
// has to be switched on for e+e-!
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log()  << "  Resetting partners..." << endl;
      }
      if(_splittingGenerator->isInteractionON(ShowerIndex::QCD)) {
	_partnerFinder->setQCDInitialEvolutionScales(showerVariables, 
						     particlesToShower);
      }
      
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "FS shower loop... " << flush << endl;
      }
      for(cit = particlesToShower.begin(); 
	  cit != particlesToShower.end(); ++cit) {
	//cout << (*cit)->id() << ", " << *cit << " " 
	//     << (*cit)->momentum()/GeV << flush << endl;
	_mapShowerHardJets[*cit] = _forwardEvolver->
	  timeLikeShower(ch, showerVariables, *cit, particles);
      }
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "...FS shower complete." << endl << flush;
      }
      if(!skipKinReco) reconstructed = reconstructKinematics( ch );
      else reconstructed = true;     
    } else reconstructed = true;
  }
  return true;
}



void Evolver::showerDecay(tEHPtr ch, 
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


void Evolver::showerGlobally(tEHPtr & ch, 
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
    if ( (*pit)->data().id() == ParticleID::g && (*pit)->isFinalState()) {
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

    
bool Evolver::reconstructISKinematics( tEHPtr & ch ) 
  throw (Veto, Stop, Exception) {
  bool ok = _kinematicsReconstructor->
    reconstructHardISJets(_mapShowerHardJets);  
  setDoneMapShower(_mapShowerHardJets);
  return ok;   
}


bool Evolver::reconstructKinematics( tEHPtr & ch ) 
  throw (Veto, Stop, Exception) 
{
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
  //  cerr << "mapshowerhardjets.size() " << _mapShowerHardJets.size() << '\n';

  bool ok = _kinematicsReconstructor->
    reconstructHardJets( _mapShowerHardJets,
			 beamHadrons.first->momentum(),
			 beamHadrons.second->momentum());
  setDoneMapShower(_mapShowerHardJets);
  return ok;   
}

void Evolver::makeRemnants(ShowerParticleVector & particles)
{
  //  cerr << "================ makeRemnants() ============\n";
  //cerr << *generator()->currentEvent() << '\n';

  // loop over the jets
  ShowerParticleVector::iterator it;
  for(it = particles.begin(); it != particles.end(); ++it )
    {
      // if initial-state, from hard process and
      if(!(*it)->isFinalState()&&(*it)->isFromHardSubprocess()&&
	 splittingGenerator()->isISRadiationON())
	{
	  tShowerParticlePtr current(*it),next;
	  Lorentz5Momentum ptotal(current->momentum());
	  do
	    {
	      // find the parent of the particle
	      if(!current->parents().empty())
		{
		  next=dynamic_ptr_cast<tShowerParticlePtr>(current->parents()[0]);
		  if(next)
		    {
		      tParticleSet siblings=current->siblings();
		      tParticleSet::iterator sib;
		      for(sib=siblings.begin();sib!=siblings.end();++sib)
			{ptotal+=(*sib)->momentum();}
		      current=next;
		    }
		}
	      else
		next=tShowerParticlePtr();
	    }
	  while(next);
	  // now we have the total momentum of the shower and are at the top
	  // of the tree so find the remnant and update it
	  tParticleSet siblings=current->siblings();
	  tParticleSet::iterator sib;
	  for(sib=siblings.begin();sib!=siblings.end();++sib)
	    {
	      // if the remnant perform cast and update
	      if((*sib)->id()==ExtraParticleID::Remnant)
		{
		  tRemnantPtr rem=dynamic_ptr_cast<tRemnantPtr>(*sib);
		  if(rem)
		    {rem->regenerate(current,current->parents()[0]->momentum()-ptotal);}
		}
	    }
	}
    }
  tParticleSet remnantSet
    = generator()->currentEventHandler()->currentCollision()->getRemnants();
  tParticleSet::iterator rem;
}
