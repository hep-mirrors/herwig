// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGGVertex class.
//

#include "SMGGGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

SMGGGGVertex::~SMGGGGVertex() {}

void SMGGGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM;
}

void SMGGGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM;
}

ClassDescription<SMGGGGVertex>
SMGGGGVertex::initSMGGGGVertex;
// Definition of the static class description member.

void SMGGGGVertex::Init() {

  
  static Reference<SMGGGGVertex,StandardModelBase> interfaceSM
    ("StandardModel",
     "Reference to the Standard Model",
     &SMGGGGVertex::_theSM, false, false, true, false);
  
  static ClassDocumentation<SMGGGGVertex> documentation
    ("The \\classname{SMGGGGVertex} class is the implementation of the"
     " Standard Model quartic gluon coupling");
  
}


// couplings for the GGGG vertex
void SMGGGGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				      tcPDPtr c,tcPDPtr d)
{
  // set the order and type
  setType(1);setOrder(0,1,2,3);
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaS(q2);
      _couplast = 4.0*3.14159265*alpha;
      _q2last=q2;
    }
  setNorm(_couplast);
}

}
}

