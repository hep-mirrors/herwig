// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFWVertex class.
//

#include "SMFFWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
void SMFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _theCKM << _theSM << _ckm;
}
  
void SMFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theCKM >> _theSM >> _ckm;
  _couplast=0.;_q2last=0.;
}
  
ClassDescription<SMFFWVertex>SMFFWVertex::initSMFFWVertex;
  
// Definition of the static class description member.
  
void SMFFWVertex::Init() {
  static ClassDocumentation<SMFFWVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the W boson to the Standard Model fermions");
  
}
  
// coupling for FFW vertex
void SMFFWVertex::setCoupling(Energy2 q2, tcPDPtr a, tcPDPtr b, tcPDPtr c)
{
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      double sw    = sqrt(2.*(_theSM->sin2ThetaW()));
      _couplast = -sqrt(4.0*3.14159265*alpha)/sw;
      _q2last=q2;
    }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  // quarks
  if(iferm>=1 && iferm <=6)
    {
      int iu,id;
      // up type first
      if(iferm%2==0)
        {
          iu = iferm/2;
          id = (ianti+1)/2;
        }
      // down type first
      else
        {
          iu = ianti/2;
          id = (iferm+1)/2;
        }
      if( iu<1 || iu>3 || id<1 || id>3)
        {
          throw HelicityConsistencyError() << "SMFFWVertex::setCoupling "
      				     << "Unknown particle in W vertex" 
      				     << Exception::warning;
          setLeft(0.);setRight(0.);return;
        }
      setLeft(_ckm[iu-1][id-1]);
      setRight(0.);
    }
  // leptons
  else if(iferm>=11 && iferm <=16)
    {
      setLeft(1.);
      setRight(0.);
    }
  else
    {
      throw HelicityConsistencyError() << "SMFFWVertex::setCoupling "
      				 << "Unknown particle in W vertex" 
      				 << Exception::warning;
      setLeft(0.);setRight(0.);
    }
}
  
}
}





