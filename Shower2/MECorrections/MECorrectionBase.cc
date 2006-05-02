// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECorrectionBase class.
//

#include "MECorrectionBase.h"
#include "Herwig++/Shower2/Kinematics/ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MECorrectionBase.tcc"
#endif


using namespace Herwig;

MECorrectionBase::~MECorrectionBase() {}

AbstractClassDescription<MECorrectionBase> MECorrectionBase::initMECorrectionBase;
// Definition of the static class description member.

void MECorrectionBase::Init() {

  static ClassDocumentation<MECorrectionBase> documentation
    ("The MECorrectionBase class is base class for implementation of the classic"
     " matrix element correction in Herwig++");

}

void MECorrectionBase::persistentOutput(PersistentOStream & os) const {
  os << _showerVariables;
}

void MECorrectionBase::persistentInput(PersistentIStream & is, int) {
  is >> _showerVariables;
}
