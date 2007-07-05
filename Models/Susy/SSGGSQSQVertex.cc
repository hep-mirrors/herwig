// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGGSQSQVertex class.
//

#include "SSGGSQSQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

SSGGSQSQVertex::SSGGSQSQVertex() : _q2last(),_couplast(0.) {
  vector<int> first,second,third,fourth;
  //L-L squarks
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(21);
    second.push_back(21);
    third.push_back(ix);
    fourth.push_back(-ix);
  }
  //R-R squarks
  for(unsigned int ix=2000001;ix<2000007;++ix) {
    first.push_back(21);
    second.push_back(21);
    third.push_back(ix);
    fourth.push_back(-ix);
  }
  setList(first,second,third,fourth);
}

void SSGGSQSQVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS;
}

void SSGGSQSQVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS;
  _couplast = 0.;
  _q2last = 0.*GeV2;
}

ClassDescription<SSGGSQSQVertex> SSGGSQSQVertex::initSSGGSQSQVertex;
// Definition of the static class description member.

void SSGGSQSQVertex::Init() {

  static ClassDocumentation<SSGGSQSQVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");

}

void SSGGSQSQVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				 tcPDPtr) { 
  if(q2 != _q2last) {
    double alphaStr = _theSS->alphaS(q2);
    _couplast = 4.*Constants::pi*alphaStr;
  }
  setNorm(_couplast);
}
