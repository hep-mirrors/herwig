// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFZVertex class.
//

#include "SMFFZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

SMFFZVertex::~SMFFZVertex() {}

void SMFFZVertex::persistentOutput(PersistentOStream & os) const {
  for(int ix=0;ix<17;++ix){os << _gl[ix] << _gr[ix];}
  os <<  _theSM;
}

void SMFFZVertex::persistentInput(PersistentIStream & is, int) {
  for(int ix=0;ix<17;++ix){is >> _gl[ix] >> _gr[ix];}
  is >> _theSM;
  _couplast=0.;_q2last=0.;
}

ClassDescription<SMFFZVertex> 
SMFFZVertex::initSMFFZVertex;
// Definition of the static class description member.

void SMFFZVertex::Init() {
  
  static Reference<SMFFZVertex,StandardModelBase> interfaceSM
    ("StandardModel",
     "Reference to the SM object",
     &SMFFZVertex::_theSM, false, false, true, false);
    
  static ClassDocumentation<SMFFZVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the Z boson to the Standard Model fermions");
}

// coupling for FFZ vertex
void SMFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c)
{
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = -sqrt(4.0*3.14159265*alpha);
      _q2last=q2;
    }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
    {
      setLeft(_gl[iferm]);
      setRight(_gr[iferm]);
    }
  else
    {throw HelicityConsistencyError() << "SMFFZVertex::setCoupling "
				      << "Unknown particle in Z vertex" 
				      << Exception::warning;
      setLeft(0.);setRight(0.);
    }
}
  
}
}

