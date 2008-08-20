// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEBase class.
//

#include "UEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


AbstractClassDescription<UEBase> UEBase::initUEBase;
// Definition of the static class description member.

void UEBase::Init() {

  static ClassDocumentation<UEBase> documentation
    ("There is no documentation for the UEBase class");

}

