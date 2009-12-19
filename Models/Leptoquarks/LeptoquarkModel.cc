// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModel class.
//

#include "LeptoquarkModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LeptoquarkModel::LeptoquarkModel() {}

LeptoquarkModel::~LeptoquarkModel() {}

IBPtr LeptoquarkModel::clone() const {
  return new_ptr(*this);
}

IBPtr LeptoquarkModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LeptoquarkModel::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void LeptoquarkModel::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<LeptoquarkModel> LeptoquarkModel::initLeptoquarkModel;
// Definition of the static class description member.

void LeptoquarkModel::Init() {

  static ClassDocumentation<LeptoquarkModel> documentation
    ("There is no documentation for the LeptoquarkModel class");

}

