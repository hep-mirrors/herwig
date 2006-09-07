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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FortranReconstructor.tcc"
#endif


using namespace Herwig;

NoPIOClassDescription<FortranReconstructor> FortranReconstructor::initFortranReconstructor;
// Definition of the static class description member.


void FortranReconstructor::Init() {

  static ClassDocumentation<FortranReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

}
