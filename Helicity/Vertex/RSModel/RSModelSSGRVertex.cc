// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelSSGRVertex class.
//

#include "RSModelSSGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

RSModelSSGRVertex::~RSModelSSGRVertex() {}

void RSModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << _theKappa;
}

void RSModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _theKappa;
}

ClassDescription<RSModelSSGRVertex> RSModelSSGRVertex::initRSModelSSGRVertex;
// Definition of the static class description member.

void RSModelSSGRVertex::Init() {
  
  static Reference<RSModelSSGRVertex,StandardModelBase> interfaceModel
    ("Model",
     "Reference to the Model object",
     &RSModelSSGRVertex::_theModel, false, false, true, false, false);
  
  static ClassDocumentation<RSModelSSGRVertex> documentation
    ("The RSModelSSGRVertex class is the implementation of"
     " the RSModel scalar-scalar-graviton vertex");
  
}
void RSModelSSGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c)
{setNorm(_theKappa);}

}
}

