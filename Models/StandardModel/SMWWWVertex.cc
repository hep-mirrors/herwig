// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWVertex class.
//

#include "SMWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

void SMWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _zfact; 
}

void SMWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _zfact;
  _couplast=0.;_q2last=0.*GeV2;
}

ClassDescription<SMWWWVertex>
SMWWWVertex::initSMWWWVertex;
// Definition of the static class description member.

void SMWWWVertex::Init() {
  static ClassDocumentation<SMWWWVertex> documentation
    ("The SMWWWVertex class is the implementation of the "
     "Standard Model triple electroweak boson coupling.");
  
}
    
// couplings for the WWW vertex
void SMWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c)
{
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = sqrt(4.0*Constants::pi*alpha);
      _q2last=q2;
    }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )
    {setNorm(_couplast);}
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) )
    {setNorm(-_couplast);}
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) || 
          (ida== 23 && idb==-24 && idc== 24) || 
          (ida== 24 && idb== 23 && idc==-24) )
    {setNorm(_couplast*_zfact);}
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) || 
          (ida== 23 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 23 && idc== 24) )
    {setNorm(-_couplast*_zfact);}
  else
    {
      throw HelicityConsistencyError() << "SMWWWVertex::setCoupling "
				       << "Invalid particles in WWW Vertex" 
				       << Exception::warning;
    }
}

}
}
