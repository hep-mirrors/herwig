// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHGGVertex class.
//

#include "SMHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM 
     << ounit(_mw,GeV) 
     << _sw << _minloop << _maxloop;
}

void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM 
     >> iunit(_mw,GeV) 
     >> _sw >> _minloop >> _maxloop;
  _couplast =Complex(0.);
  _q2last = 0.*GeV2;
}

ClassDescription<SMHGGVertex> SMHGGVertex::initSMHGGVertex;
// Definition of the static class description member.

void SMHGGVertex::Init() {
  
  static ClassDocumentation<SMHGGVertex> documentation
    ("This class implements the h->g,g vertex");

  static Parameter<SMHGGVertex,unsigned int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHGGVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHGGVertex,unsigned int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHGGVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);
}
  
void SMHGGVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, tcPDPtr part3) {
  int delta = _maxloop - _minloop;
  if (0>delta) {
    _maxloop=_maxloop-delta;
    _minloop=_minloop+delta;
    delta*=-1;
  }
  ++delta;

  type.resize(delta,PDT::SpinUnknown);
  masses.resize(delta,0.*GeV);
  left.resize(delta,0.);
  right.resize(delta,0.);
  for (int i = 0; i < delta; ++i) {
    tcPDPtr q = getParticleData(_minloop+i);
    type[i] = PDT::Spin1Half;
    masses[i] = q->mass();
    left[i] = q->mass()/_mw;
  }
  right=left;

  if(q2 != _q2last) {
    double g = sqrt(4.*Constants::pi*_theSM->alphaEM(q2)/_theSM->sin2ThetaW());
    double gs2 = 4.*Constants::pi*_theSM->alphaS(q2);
    _couplast = Complex(UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi));
    _q2last = q2;
  }
  setNorm(_couplast);

  //calculate tensor coefficients
//  SVVLoopVertex::setCoupling(q2,part1,part2,part3);
  SimpleSVVLoopVertex::setCoupling(q2,part1,part2,part3);
}

}
 
