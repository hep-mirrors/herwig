// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayVertex class.
//
//  Author: Peter Richardson
//

#include "DecayVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/SpinInfo.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayVertex.tcc"
#endif

namespace Herwig {

using ThePEG::Helicity::SpinInfo;
using ThePEG::Helicity::tcSpinfoPtr;

namespace Helicity {
using namespace ThePEG;

DecayVertex::~DecayVertex() {}
    
NoPIOClassDescription<DecayVertex> DecayVertex::initDecayVertex;
  // Definition of the static class description member.
    
void DecayVertex::Init() {
  
  static ClassDocumentation<DecayVertex> documentation
    ("The DecayVertex class is the implementation of a "
     "vertex for a decay for the Herwig++ spin correlation algorithm");
  
}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix DecayVertex::getRhoMatrix(int i)
{
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    {
      if(int(ix)!=i)
	{rhoout.push_back(dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix());}
    }
  // calculate the spin density matrix
  RhoDMatrix input=dynamic_ptr_cast<tcSpinfoPtr>(incoming()[0])->rhoMatrix();
  RhoDMatrix temp=_matrixelement.calculateRhoMatrix(i,input,rhoout);
  return temp;
}

// method to get the D matrix for an incoming particle
RhoDMatrix DecayVertex::getDMatrix(int i)
{
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    {
      Dout.push_back(dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix());
    }
  // calculate the spin density matrix and return the answer
  RhoDMatrix temp = _matrixelement.calculateDMatrix(Dout);
  return temp;
    }
}

}
