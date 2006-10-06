// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelVVGRVertex class.
//

#include "RSModelVVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

void RSModelVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << _theKappa;
}

void RSModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _theKappa;
}

ClassDescription<RSModelVVGRVertex> RSModelVVGRVertex::initRSModelVVGRVertex;
// Definition of the static class description member.

void RSModelVVGRVertex::Init() {
 static ClassDocumentation<RSModelVVGRVertex> documentation
    ("The RSModelVVGRVertex class is the implementation"
     " of the RSModel vector-vector-graviton vertex");
  
}
  
void RSModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr)
{setNorm(_theKappa);}

}
}

