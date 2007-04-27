// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NBodyDecayConstructorBase class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig; 
using namespace ThePEG;

void NBodyDecayConstructorBase::persistentOutput(PersistentOStream & ) const {
}

void NBodyDecayConstructorBase::persistentInput(PersistentIStream & , int) {
}

AbstractClassDescription<NBodyDecayConstructorBase> NBodyDecayConstructorBase::initNBodyDecayConstructorBase;
// Definition of the static class description member.

void NBodyDecayConstructorBase::Init() {

  static ClassDocumentation<NBodyDecayConstructorBase> documentation
    ("There is no documentation for the NBodyDecayConstructorBase class");

}

