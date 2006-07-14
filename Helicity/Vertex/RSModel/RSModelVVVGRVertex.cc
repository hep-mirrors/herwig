// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelVVVGRVertex class.
//

#include "RSModelVVVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

void RSModelVVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << _theKappa << _zfact;
}
void RSModelVVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _theKappa >> _zfact;
  for(int ix=0;ix<2;++ix){_couplast[ix]=0.;_q2last[ix]=0.;}
}

ClassDescription<RSModelVVVGRVertex> RSModelVVVGRVertex::initRSModelVVVGRVertex;
// Definition of the static class description member.

void RSModelVVVGRVertex::Init() {
  
  
  static Reference<RSModelVVVGRVertex,StandardModelBase> interfaceModel
    ("Model",
     "Reference to the model object",
     &RSModelVVVGRVertex::_theModel, false, false, true, false, false);
  
  static ClassDocumentation<RSModelVVVGRVertex> documentation
    ("The RSModelVVVGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the VVVGR vertex
void RSModelVVVGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				 tcPDPtr c, tcPDPtr d)
{
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(ida==21 && idb==21 && idc==21)
    {
  if(q2!=_q2last[1])
    {
      double alphaS = _theModel->alphaS(q2);
      _couplast[1] = sqrt(4.0*3.14159265*alphaS);
      _q2last[1]=q2;
    }
  setNorm(_couplast[1]*_theKappa);
    }
  else
    {
      if(q2!=_q2last[0])
        {
          double alpha = _theModel->alphaEM(q2);
          _couplast[0] = sqrt(4.0*3.14159265*alpha);
          _q2last[0]=q2;
        }
      // W- W+ photon and cylic perms
      if((ida==-24 && idb== 24 && idc== 22) ||
         (ida== 22 && idb==-24 && idc== 24) || 
         (ida== 24 && idb== 22 && idc==-24) )
        {setNorm(_couplast[0]*_theKappa);}
      // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) ||
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) )
    {setNorm(-_couplast[0]*_theKappa);}
      // W- W+ Z and cylic perms
      else if((ida==-24 && idb== 24 && idc== 23) ||
    	  (ida== 23 && idb==-24 && idc== 24) || 
    	  (ida== 24 && idb== 23 && idc==-24) )
        {setNorm(_couplast[0]*_zfact*_theKappa);}
      // W+ W- Z (anticylic perms of above)
      else if((ida== 24 && idb==-24 && idc== 23) ||
    	  (ida== 23 && idb== 24 && idc==-24) || 
    	  (ida==-24 && idb== 23 && idc== 24) )
        {setNorm(-_couplast[0]*_zfact*_theKappa);}
      else
        {
	  throw HelicityConsistencyError() << "RSModelVVVGRVertex::setCoupling " 
					   << "Invalid particles in VVVGR Vertex" 
					   << Exception::warning;
	  setNorm(0.);
        }
    }
}
    
}
}
