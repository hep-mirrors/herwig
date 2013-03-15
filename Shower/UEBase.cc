// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEBase class.
//

#include "UEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractClass<UEBase,Interfaced>
describeUEBase ("Herwig::UEBase","HwShower.so");

void UEBase::Init() {

  static ClassDocumentation<UEBase> documentation
    ("The UEBase class is an abstract base class used to minimize the"
     " dependence between the MPIHandler and all Shower classes");

}

