// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGVertex class.
//

#include "SMGGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;

void SMGGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM;
}

void SMGGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM;
}

ClassDescription<SMGGGVertex> SMGGGVertex::initSMGGGVertex;
// Definition of the static class description member.

void SMGGGVertex::Init() {
 static ClassDocumentation<SMGGGVertex> documentation
    ("The SMGGGVertex class is the implementation"
     " of the Standard Model triple gluon vertex.");
  
}

// couplings for the GGG vertex
void SMGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr, tcPDPtr)
{
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaS(q2);
      _couplast = sqrt(4.0*Constants::pi*alpha);
      _q2last=q2;
    }
  setNorm(_couplast);
}

}
}

