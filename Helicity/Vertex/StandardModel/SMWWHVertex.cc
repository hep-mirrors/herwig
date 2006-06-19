// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WWHVertex class.
//
#include "SMWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
void SMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _mw << _zfact << _sw;
}

void SMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _mw >> _zfact >> _sw;
  _couplast=0.; _q2last=0.;
}
    
ClassDescription<SMWWHVertex>SMWWHVertex::initSMWWHVertex;
// Definition of the static class description member.

void SMWWHVertex::Init() {
            
  static Reference<SMWWHVertex,StandardModelBase> interfaceSM
    ("StandardModel",
     "Reference to the SM object",
     &SMWWHVertex::_theSM, false, false, true, false);
  
  
  static ClassDocumentation<SMWWHVertex> documentation
    ("The SMWWHVertex class is the implementation"
     " of the helicity amplitude calculation for the coupling of the Standard"
     " Model electroweak gauge bosons to the Higgs.");
}


void SMWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c)
{
  int ibos=abs(a->id());
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = sqrt(4.0*3.14159265*alpha)*_mw/_sw;
      _q2last=q2;
    }
  if(ibos==24){setNorm(_couplast);}
  else if(ibos==23){setNorm(_couplast*_zfact);}
  else
    {
      throw HelicityConsistencyError() << "SMWWHVertex::setCoupling "
				       << "Invalid particles in WWH Vertex" 
				       << Exception::warning;
      setNorm(0.);
    }
}
  
}
}

