// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranReconstructor class.
//

#include "FortranReconstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Timer.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"

using namespace Herwig;

NoPIOClassDescription<FortranReconstructor> FortranReconstructor::initFortranReconstructor;
// Definition of the static class description member.

void FortranReconstructor::Init() {

  static ClassDocumentation<FortranReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

}

bool FortranReconstructor::reconstructHardJets(ShowerTreePtr hard) const {
  throw Exception() << "FortranReconstructor::reconstructHardJets()"
		    << " not implemented yet" << Exception::runerror;
}

bool FortranReconstructor::reconstructDecayJets(ShowerTreePtr decay) const {
  throw Exception() << "FortranReconstructor::reconstructDecayJets()"
		    << " not implemented yet" << Exception::runerror;
}
