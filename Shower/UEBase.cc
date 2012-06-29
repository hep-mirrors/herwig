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

AbstractClassDescription<UEBase> UEBase::initUEBase;
// Definition of the static class description member.

void UEBase::Init() {

  static ClassDocumentation<UEBase> documentation
    ("The UEBase class is an abstract base class used to minimize the"
     " dependence between the MPIHandler and all Shower classes");

}

