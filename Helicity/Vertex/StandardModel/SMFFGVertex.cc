// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFGVertex class.
//

#include "SMFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

void SMFFGVertex::persistentOutput(PersistentOStream & os) const {
  os <<  _theSM;
}

void SMFFGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM;
  _couplast=0.;_q2last=0.;
}

ClassDescription<SMFFGVertex> 
SMFFGVertex::initSMFFGVertex;
// Definition of the static class description member.

void SMFFGVertex::Init() {
  static ClassDocumentation<SMFFGVertex> documentation
    ("The SMFFGVertex class is the implementation of"
     "the coupling of the gluon to the Standard Model fermions");
  
}

// coupling for FFG vertex
void SMFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c)
{
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alphas = _theSM->alphaS(q2);
      _couplast = -sqrt(4.0*3.14159265*alphas);
      _q2last=q2;
    }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if(iferm>=1 && iferm<=6)
    {
      setLeft(1.);
      setRight(1.);
  }
  else
    {
      throw HelicityConsistencyError() << "SMFFGVertex::setCoupling" 
				       << "Unknown particle in gluon vertex" 
				       << Exception::warning;
      setLeft(0.);setRight(0.);
    }
}

}
}
