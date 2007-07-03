// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHPPVertex class.
//

#include "SMHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
void SMHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _sw;
}

void SMHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _sw;
  _couplast = 0.;
  _q2last = 0.*GeV2;
}

ClassDescription<SMHPPVertex> SMHPPVertex::initSMHPPVertex;
// Definition of the static class description member.

void SMHPPVertex::Init() {

  static ClassDocumentation<SMHPPVertex> documentation
    ("This class implements the h0->gamma,gamma vertex.");
  
}

void SMHPPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
                              tcPDPtr part2, tcPDPtr part3) {
  type.resize(2,PDT::SpinUnknown);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1;
  masses.resize(2,0.*GeV);
  left.resize(2,0.);right.resize(2,0.);
  masses[0] = _theSM->mass(q2,getParticleData(ParticleID::t));
  masses[1] = _mw;
  left[0] = -_theSM->mass(q2,getParticleData(6))*(4./3.)/_mw/2.;
  left[1] =  _mw * UnitRemoval::InvE;
  right = left;
  if(q2 != _q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = 4.*Constants::pi*alpha*sqrt(4*Constants::pi*alpha)/_sw;
    _q2last = q2;
  }
  setNorm(_couplast);
  SVVLoopVertex::setCoupling(q2,part1,part2,part3);
}
  
}
}
