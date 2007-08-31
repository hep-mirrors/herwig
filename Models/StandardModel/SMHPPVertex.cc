// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHPPVertex class.
//

#include "SMHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;

void SMHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << massopt << _minloop << _maxloop << _CoefRepresentation;
}

void SMHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw, GeV) >> massopt >> _minloop >> _maxloop >> _CoefRepresentation;
  setCoefScheme(_CoefRepresentation);
}

ClassDescription<SMHPPVertex> SMHPPVertex::initSMHPPVertex;
// Definition of the static class description member.

void SMHPPVertex::Init() {
  static ClassDocumentation<SMHPPVertex> documentation
    ("This class implements the h0->gamma,gamma vertex.");

  static Parameter<SMHPPVertex,unsigned int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHPPVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHPPVertex,unsigned int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHPPVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHPPVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHPPVertex::massopt, 2, false, false);
  static SwitchOption interfaceHeavyMass
    (interfaceMassOption,
     "PoleMasses",
     "The loop is calculcated with the pole quark masses",
     1);
  static SwitchOption interfaceNormalMass
    (interfaceMassOption,
     "RunningMasses",
     "running quark masses are taken in the loop",
     2);

  static Switch<SMHPPVertex,unsigned int> interfaceScheme
    ("CoefficientScheme",
     "Which scheme for the tensor coefficients is applied",
     &SMHPPVertex::_CoefRepresentation, 1, false, false);
  static SwitchOption interfaceSchemeSimplified
    (interfaceScheme,
     "Simplified",
     "Represection suitable for the simplified the H-g-g and H-gamma-gamma vertices",
     1);
  static SwitchOption interfaceSchemeGeneral
    (interfaceScheme,
     "General",
     "Represection suitable for the Passarino-Veltman tensor reduction scheme",
     2);
}


void SMHPPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
                              tcPDPtr part2, tcPDPtr part3) {
  if( part1->id() != ParticleID::h0 && 
      part2->id() != ParticleID::gamma &&
      part3->id() != ParticleID::gamma ) {
    throw HelicityConsistencyError() 
      << "SMHPPVertex::setCoupling() - The particle content of this vertex "
      << "is incorrect: " << part1->id() << " " << part2->id() << " "
      << part3->id() << Exception::warning;
    setNorm(0.);
    return;
  }

  orderInGs(0);
  orderInGem(3);

  unsigned int Qminloop = _minloop;
  unsigned int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }

  switch (_CoefRepresentation) {
    case 1: {
      if(q2 != _q2last) {
        double g = sqrt(4.*Constants::pi*_theSM->alphaEM(q2)/_theSM->sin2ThetaW());
        _couplast = UnitRemoval::E * pow(g,3)/_mw/sqr(Constants::pi)/sqrt(2.)/16.;
        _q2last = q2;
      }
      setNorm(_couplast);

      Complex loop(0.);
      // quark loops
      for (unsigned int i = Qminloop; i <= Qmaxloop; ++i) {
        tcPDPtr qrk = getParticleData(i);
        Energy mass = qrk->mass();
        Charge charge = qrk->charge();
        if (2 == massopt) {
          mass = _theSM->mass(q2,qrk);
        }
        loop += sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/q2);
      }
      // lepton loops
      unsigned int Lminloop = 3; // still fixed value
      unsigned int Lmaxloop = 3; // still fixed value
      for (unsigned int i = Lminloop; i <= Lmaxloop; ++i) {
        tcPDPtr lpt = getParticleData(9 + 2*i);
        Energy mass = lpt->mass();
        Charge charge = lpt->charge();
        if (2 == massopt) {
          mass = _theSM->mass(q2,lpt);
        }
        loop += sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/q2)/3.;  // 3. -> no color!
      }

      // W loop
      loop += Aw(sqr(_mw)/q2);

      a00(loop);
      a11(0.0);
      a12(0.0);
      a21(-loop);
      a22(0.0);
      aEp(0.0);
      break;}
    case 2: {
      if(q2 != _q2last) {
        double e = sqrt(4.*Constants::pi*_theSM->alphaEM(q2));
        _couplast = pow(e,3)/_theSM->sin2ThetaW();
        _q2last = q2;
      }
      setNorm(_couplast);

      type.resize(3,PDT::SpinUnknown);
      type[0] = PDT::Spin1Half;
      type[1] = PDT::Spin1Half;
      type[2] = PDT::Spin1;
      masses.resize(3,0.*GeV);
      masses[0] = _theSM->mass(q2,getParticleData(ParticleID::t));
      masses[1] = _theSM->mass(q2,getParticleData(ParticleID::b));
      masses[2] = _mw;
      double copl = -_theSM->mass(q2,getParticleData(6))*(4./3.)/_mw/2.;
      couplings.push_back(make_pair(copl, copl));
      copl = -_theSM->mass(q2,getParticleData(5))*(4./3.)/_mw/2.;
      couplings.push_back(make_pair(copl, copl));
      couplings.push_back(make_pair(UnitRemoval::InvE*_mw, UnitRemoval::InvE*_mw));

      SVVLoopVertex::setCoupling(q2, part1, part2, part3);
      break;}
  }

  return;
}

Complex SMHPPVertex::Af(const double tau) const {
  return tau*(4. - W2(tau)*(1. - 4.*tau));
}

Complex SMHPPVertex::Aw(const double tau) const {
  return 12.*W2(tau)*tau*(4.*tau - 2.) - 12.*tau - 2.;
}

Complex SMHPPVertex::W2(double lambda) const {
  double pi = Constants::pi;

  if (0.0 == lambda) 
    return 0.0;

  if (lambda < 0.0) 
    return 4.*sqr(asinh(0.5*sqrt(-1./lambda)));

  double root(0.5*sqrt(1./lambda));
  Complex ac(0.);
  if(root < 1.) {
    ac = -sqr(asin(root));
  } else {
// formulae from NPB297,221
    double ex = acosh(root);
    ac = sqr(ex) - 0.25*sqr(pi) - pi*ex*Complex(0.,1.);
/*
// formulae from Higgs hunter's guide (gives the same result).
    double pl = .5 + .5*sqrt(1. - 4.*lambda);
    double ms = .5 - .5*sqrt(1. - 4.*lambda);
    double lg = 0.5*log(pl/ms);
    ac = sqr(lg) - 0.25*sqr(pi) - pi*lg*Complex(0.,1.);
*/
  }
  return 4.*ac;
}

}
