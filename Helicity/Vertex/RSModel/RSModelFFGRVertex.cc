// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGRVertex class.
//

#include "RSModelFFGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;

RSModelFFGRVertex::~RSModelFFGRVertex() {}

void RSModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << _theKappa;
}

void RSModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _theKappa;
}

ClassDescription<RSModelFFGRVertex> RSModelFFGRVertex::initRSModelFFGRVertex;
// Definition of the static class description member.

void RSModelFFGRVertex::Init() {
  
  
  static Reference<RSModelFFGRVertex,StandardModelBase> interfaceModel
    ("Model",
     "Reference to the Model object",
     &RSModelFFGRVertex::_theModel, false, false, true, false, false);
  
  
  static ClassDocumentation<RSModelFFGRVertex> documentation
    ("The \\classname{RSModelFFGRVertex} class is the RSModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
// couplings
void RSModelFFGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c)
{setNorm(_theKappa);}

}
}

