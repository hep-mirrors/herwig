// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KinematicsReconstructor class.
//

#include "KinematicsReconstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

AbstractClassDescription<KinematicsReconstructor> KinematicsReconstructor::initKinematicsReconstructor;
// Definition of the static class description member.

void KinematicsReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _showerVariables;
}

void KinematicsReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _showerVariables;
}

void KinematicsReconstructor::Init() {

  static ClassDocumentation<KinematicsReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

}
