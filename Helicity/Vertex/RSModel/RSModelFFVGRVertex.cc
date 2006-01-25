// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFVGRVertex class.
//

#include "RSModelFFVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
RSModelFFVGRVertex::~RSModelFFVGRVertex() {}

void RSModelFFVGRVertex::persistentOutput(PersistentOStream & os) const {
  for(int ix=0;ix<17;++ix){os << _charge[ix];}
  os <<  _theModel << _theKappa;
}

void RSModelFFVGRVertex::persistentInput(PersistentIStream & is, int) {
  for(int ix=0;ix<17;++ix){is >> _charge[ix];}
  is >> _theModel >> _theKappa;
  for(int ix=0;ix<2;++ix){_couplast[ix]=0.;_q2last[ix]=0.;}
}

ClassDescription<RSModelFFVGRVertex> RSModelFFVGRVertex::initRSModelFFVGRVertex;
// Definition of the static class description member.

void RSModelFFVGRVertex::Init() {
  
  static Reference<RSModelFFVGRVertex,StandardModelBase> interfaceModel
    ("Model",
     "Reference to the Model object",
     &RSModelFFVGRVertex::_theModel, false, false, true, false);

  static ClassDocumentation<RSModelFFVGRVertex> documentation
    ("The RSModelFFVGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}
// FFVGR coupling
void RSModelFFVGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				 tcPDPtr c, tcPDPtr d)
{
  // work out the particles
  int iferm=abs(a->id());
  int ibos =c->id();
  Complex norm;
  // overall factor
  // photon
  if(ibos==22)
    {
      // alpha
      if(_q2last[0]!=q2)
        {
          double alpha = _theModel->alphaEM(q2);
          _couplast[0] = sqrt(4.0*3.14159265*alpha);
          _q2last[0]=q2;
        }
      norm = _theKappa*_couplast[0];
      // _charge of particle
      if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
        {norm = norm*_charge[iferm];}
      else
        {
	  throw HelicityConsistencyError() << "RSModelFFVGRVertex::setCoupling " 
					   << "Unknown particle in FFVGR vertex" 
					   << Exception::warning;
	  setNorm(0.);
	  return;
	}
    }
  // gluon
  else if (ibos==21||ibos==9)
    {
      if(_q2last[1]!=q2)
        {
          double alphas= _theModel->alphaS(q2);
          _couplast[1] = sqrt(4.0*3.14159265*alphas);
          _q2last[1]=q2;
        }
      norm = _theKappa*_couplast[1];
    }
  else
    {
      throw HelicityConsistencyError() << "RSModelFFVGRVertex::setCoupling " 
				       << "Unknown boson in FFVGR vertex" 
				       << Exception::warning;
      setNorm(0.);
      return;
    }
  // set the coupling
  setNorm(norm);
}

}
}

