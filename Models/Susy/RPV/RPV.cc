// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPV class.
//

#include "RPV.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RPV::persistentOutput(PersistentOStream & ) const {
}

void RPV::persistentInput(PersistentIStream & , int) {
}

ClassDescription<RPV> RPV::initRPV;
// Definition of the static class description member.

void RPV::Init() {

  static ClassDocumentation<RPV> documentation
    ("There is no documentation for the RPV class");

}

