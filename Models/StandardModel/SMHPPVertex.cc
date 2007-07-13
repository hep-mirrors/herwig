// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHPPVertex class.
//

#include "SMHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void SMHPPVertex::doinit() throw(InitException) {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _top = getParticleData(ThePEG::ParticleID::t);
  _sw = sqrt(_theSM->sin2ThetaW());
  //resize vectors to correct size so that it is not done for each pass 
  //through setCoupling
  masses.resize(2);
  couplings.resize(2);
  type.resize(2);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1;
  orderInGs(0);
  orderInGem(3);
  SVVLoopVertex::doinit();
}
    
void SMHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _sw << _top << _haveCoeff;
}

void SMHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw, GeV) >> _sw >> _top >> _haveCoeff;
}

ClassDescription<SMHPPVertex> SMHPPVertex::initSMHPPVertex;
// Definition of the static class description member.

void SMHPPVertex::Init() {

  static ClassDocumentation<SMHPPVertex> documentation
    ("This class implements the h0->gamma,gamma vertex.");
  
}

void SMHPPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
                              tcPDPtr part2, tcPDPtr part3) {
  if( part1->id() != ParticleID::h0 && 
      part2->id() != ParticleID::gamma &&
      part3->id() != ParticleID::gamma ) {
    throw HelicityConsistencyError() 
      << "SMHPPVertex::setCoupling() - The particle content of this vertex "
      << "is incorrect: " << part1->id() << " " << part2->id()
      << part3->id() << Exception::warning;
    setNorm(0.);
    return;
  }
  if(q2 != _q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = 4.*Constants::pi*alpha*sqrt(4*Constants::pi*alpha)/_sw;
    Energy mt = _theSM->mass(q2,_top);
    masses[0] = mt;
    masses[1] = _mw;
    double topc = -2.*mt/3./_mw;
    couplings[0] = make_pair(topc, topc);
    couplings[1] = make_pair(_mw*UnitRemoval::InvE, _mw*UnitRemoval::InvE);
    _haveCoeff = false;
    _q2last = q2;
  }
  setNorm(_couplast);
  if( !_haveCoeff ) {
    SVVLoopVertex::setCoupling(q2,part1,part2,part3);
    _haveCoeff = true;
  }

}

