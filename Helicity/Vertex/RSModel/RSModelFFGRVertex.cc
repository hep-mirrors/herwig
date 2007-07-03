// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGRVertex class.
//

#include "RSModelFFGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

void RSModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << ounit(_theKappa,InvGeV);
}

void RSModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelFFGRVertex> RSModelFFGRVertex::initRSModelFFGRVertex;
// Definition of the static class description member.

void RSModelFFGRVertex::Init() {
  static ClassDocumentation<RSModelFFGRVertex> documentation
    ("The RSModelFFGRVertex class is the RSModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
// couplings
void RSModelFFGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr)
{setNorm(Complex(_theKappa * UnitRemoval::E));}

}
}

