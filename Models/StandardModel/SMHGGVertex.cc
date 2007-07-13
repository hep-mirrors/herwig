// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHGGVertex class.
//

#include "SMHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void SMHGGVertex::doinit() throw(InitException) {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _sw = sqrt(_theSM->sin2ThetaW());
  
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _top = getParticleData(ParticleID::t);
  _bottom = getParticleData(ParticleID::b);
  if( _qopt == 0 ) {
    setNParticles(1);
    type.resize(1, PDT::Spin1Half);
  }
  else {
    setNParticles(2);
    type.resize(2, PDT::Spin1Half);
  }
  orderInGs(2);
  orderInGem(1);
  SVVLoopVertex::doinit();
}


void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw, GeV) << _sw << _top << _qopt << _haveCoeff;
}

void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw, GeV) >> _sw >> _top >> _qopt >> _haveCoeff;
}

ClassDescription<SMHGGVertex> SMHGGVertex::initSMHGGVertex;
// Definition of the static class description member.

void SMHGGVertex::Init() {
  
  static ClassDocumentation<SMHGGVertex> documentation
    ("This class implements the h->g,g vertex");
  
  static Switch<SMHGGVertex,int> interfaceQuarks
    ("Quarks",
     "Which quark loops to include in the coupling.",
     &SMHGGVertex::_qopt, 0, false, false);
  static SwitchOption interfaceQuarksTop
    (interfaceQuarks,
     "Top",
     "Only include top loop",
     0);
  static SwitchOption interfaceQuarksBottom
    (interfaceQuarks,
     "Bottom",
     "Include top and bottom",
     1);
}
  
void SMHGGVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr part3) {
  if( part1->id() != ParticleID::h0 && 
      part2->id() != ParticleID::g &&
      part3->id() != ParticleID::g ) {
    throw HelicityConsistencyError() 
      << "SMHGGVertex::setCoupling() - The particle content of this vertex "
      << "is incorrect: " << part1->id() << " " << part2->id()
      << part3->id() << Exception::warning;
    setNorm(0.);
    return;
  }
  if(q2 != _q2last) {
    double alphaStr = _theSM->alphaS(q2);
    double alpha = _theSM->alphaEM(q2);
    _couplast = 2.*Constants::pi*alphaStr*sqrt(4*Constants::pi*alpha)/_sw;
    Energy mt = _theSM->mass(q2, _top);
    if(_qopt == 0 ) {
      masses.push_back(mt);
      couplings.push_back(make_pair(mt/_mw, mt/_mw));
    }
    else if( _qopt == 1 ) {
      Energy mb = _theSM->mass(q2, _bottom);
      masses.push_back(mb);
      couplings.push_back(make_pair(mb/_mw, mb/_mw));
    }
    else {
      throw InitException() 
	<< "SMHGGVertex::setCoupling() - Unknown option for particles "
	<< "in the loop " << _qopt
	<< Exception::warning;
      setNorm(0.);
      return;
    }
    _haveCoeff = false;
    _q2last = q2;
  }
  setNorm(_couplast);
  //calculate tensor coefficients
  if( !_haveCoeff ) {
    SVVLoopVertex::setCoupling(q2,part1,part2,part3);
    _haveCoeff = true;
  }

}
 
