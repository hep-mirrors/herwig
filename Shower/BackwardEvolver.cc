// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BackwardEvolver class.
//

#include "BackwardEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerParticle.h"
#include "ShowerKinematics.h"

using namespace Herwig;


BackwardEvolver::~BackwardEvolver() {}


void BackwardEvolver::persistentOutput(PersistentOStream & os) const {
  os << _splittingGenerator << _forwardEvolver;
}


void BackwardEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _splittingGenerator >> _forwardEvolver;
}


ClassDescription<BackwardEvolver> BackwardEvolver::initBackwardEvolver;
// Definition of the static class description member.

void BackwardEvolver::Init() {

  static ClassDocumentation<BackwardEvolver> documentation
    ("This class is responsible for the backward showering of space-like particles");

  static Reference<BackwardEvolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::BackwardEvolver::_splittingGenerator,
		      false, false, true, false);
  static Reference<BackwardEvolver,ForwardEvolver> 
    interfaceForwardEvolver("ForwardEvolver", 
			    "A reference to the ForwardEvolver object", 
			    &Herwig::BackwardEvolver::_forwardEvolver,
			    false, false, true, false);

}

//------------------------------------------------------------------------------

bool BackwardEvolver::spaceLikeShower(tPartCollHdlPtr ch, 
				      const tShowerVarsPtr showerVars, 
				      //const tMECorrectionPtr meCorrectionPtr,
				      tShowerParticlePtr particle, 
				      ShowerParticleVector &allShowerParticles)
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  bool hasEmitted = false;
  tShowerParticlePtr spaceLikePart = particle;
  tShowerParticleVector particlesYetToShower;   // only time-like particles

  do {

    Branching bb = _splittingGenerator->chooseBackwardBranching(ch, *spaceLikePart);
    //                accept it according to the  showerVariables  and soft correction;
    //                if ( does not branch ) {
    //                  rhoD propagation;
    //                  spaceLikePart = tShowerParticlePtr();
    //                } else {
    //                  hasEmitted = true;
    _splittingGenerator->generateBranchingKinematics(ch, *spaceLikePart, bb);
    //                  create the new ShowerParticles and then store the
    //                    unique space-like one in spaceLikePart, whereas 
    //                    the others are stored into  particleYetToShower;
    //                  store also the shoKin;
    //                }
    // 
    //                NB) To access the PDF:
    //                      PDF myPDF = ch.pdf( parton );
    //                    where parton is a pointer to a ThePEG particle
    //                    and the method returns a PDF object for the given
    //                    particle. Such method is defined in LastXCombInfo
    //                    from which PartialCollisionHandler inherits from.
    //***endLOOKHERE***

  } while(!spaceLikePart);

  while(!particlesYetToShower.empty()) {
    tShowerParticlePtr part = particlesYetToShower.back();
    particlesYetToShower.pop_back();
    hasEmitted = hasEmitted || 
     _forwardEvolver->timeLikeShower(ch, showerVars, part, allShowerParticles);
  } 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> END DEBUGGING <=== "
		       << endl;
  }
  return hasEmitted;
}


