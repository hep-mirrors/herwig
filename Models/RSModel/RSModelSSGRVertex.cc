// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelSSGRVertex class.
//

#include "RSModelSSGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void RSModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << ounit(_theKappa,InvGeV);
}

void RSModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelSSGRVertex> RSModelSSGRVertex::initRSModelSSGRVertex;
// Definition of the static class description member.

void RSModelSSGRVertex::Init() {
  static ClassDocumentation<RSModelSSGRVertex> documentation
    ("The RSModelSSGRVertex class is the implementation of"
     " the RSModel scalar-scalar-graviton vertex");
  
}
void RSModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr)
{setNorm(Complex(_theKappa * UnitRemoval::E));}
}

