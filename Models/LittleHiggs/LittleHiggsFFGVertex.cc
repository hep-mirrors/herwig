// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsFFGVertex class.
//

#include "LittleHiggsFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

void LittleHiggsFFGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM;
}

void LittleHiggsFFGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM;
}

ClassDescription<LittleHiggsFFGVertex> LittleHiggsFFGVertex::initLittleHiggsFFGVertex;
// Definition of the static class description member.

void LittleHiggsFFGVertex::Init() {

  static ClassDocumentation<LittleHiggsFFGVertex> documentation
    ("There is no documentation for the LittleHiggsFFGVertex class");

}

// coupling for FFG vertex
void LittleHiggsFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alphas = _theSM->alphaS(q2);
    _couplast = -sqrt(4.0*Constants::pi*alphas);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6) || iferm==8) {
    setLeft(1.);
    setRight(1.);
  }
  else
    throw HelicityConsistencyError() << "LittleHiggsFFGVertex::setCoupling" 
				     << "Unknown particle in gluon vertex" 
				     << Exception::runerror;
}

LittleHiggsFFGVertex::LittleHiggsFFGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  // PDG codes for the particles
  vector<int> first,second,third;
  for(unsigned int ix=1;ix<9;++ix) {
    if(ix==7) ++ix;
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(21);
  }
  setList(first,second,third);
}
  
void LittleHiggsFFGVertex::doinit() throw(InitException) {
  _theSM = generator()->standardModel();
  orderInGs(1);
  orderInGem(0);
  FFVVertex::doinit();
}
