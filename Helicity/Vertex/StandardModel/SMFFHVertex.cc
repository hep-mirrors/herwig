// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFHVertex class.
//

#include "SMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Config/Constants.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
using ThePEG::Constants::pi;
    
SMFFHVertex::~SMFFHVertex() {}

void SMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _mw << _sw;
}

void SMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _mw >> _sw;
  _couplast=0.;_idlast=0;_q2last=0.;_masslast=0.;
}

ClassDescription<SMFFHVertex> 
SMFFHVertex::initSMFFHVertex;
// Definition of the static class description member.

void SMFFHVertex::Init() {
  
  
  static Reference<SMFFHVertex,Herwig::StandardModel>
    interfaceSM
    ("StandardModel",
     "Reference to the SM object",
     &SMFFHVertex::_theSM, false, false, true, false);

  static ClassDocumentation<SMFFHVertex> documentation
    ("The \\classname{SMFFHVertex} class is the implementation"
     " of the helicity amplitude calculation of the Standard Model Higgs"
     " fermion-antiferiom vertex.");
  
}
void SMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c)
{
  int iferm=abs(a->id());
  // left and right couplings set to one
  setLeft(1.); setRight(1.);
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = -0.5*sqrt(4.0*pi*alpha)/_sw/_mw;
      _q2last=q2;
      _idlast=iferm;
      if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
        {_masslast=_theSM->mass(q2,a);}
      else
        {
	  throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
					   << "Unknown particle in Higgs vertex" 
					   << Exception::warning;
	  _masslast = 0;
        }
    }
  else if(iferm!=_idlast)
    {
      _idlast=iferm;
      if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
        {_masslast=_theSM->mass(q2,a);}
      else
        {
	  throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
					   << "Unknown particle in Higgs vertex" 
					   << Exception::warning;
	  _masslast = 0;
        }
    }
  setNorm(_couplast*_masslast);
} 
 
}
}
