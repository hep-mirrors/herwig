// -*- C++ -*-
//
// SMHGGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHGGVertex class.
//

#include "SMHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;  

SMHGGVertex::SMHGGVertex()
  :_couplast(0.),_q2last(ZERO),_mw(),massopt(1),_minloop(6),
   _maxloop(6),_CoefRepresentation(1) {
  //PDG codes for particles at vertices
  vector<long> first(1,21),second(1,21),third(1,25);
  setList(first,second,third);
}

void SMHGGVertex::doinit() throw(InitException) {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!_theSM) 
    throw InitException();
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  orderInGs(2);
  orderInGem(1);
  VVSLoopVertex::doinit();
}

void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << massopt 
     << _minloop << _maxloop << _CoefRepresentation;
}

void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> massopt 
     >> _minloop >> _maxloop >> _CoefRepresentation;
  setCoefScheme(_CoefRepresentation);
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

  static Switch<SMHGGVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHGGVertex::massopt, 1, false, false);
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

  static Switch<SMHGGVertex,unsigned int> interfaceScheme
    ("CoefficientScheme",
     "Which scheme for the tensor coefficients is applied",
     &SMHGGVertex::_CoefRepresentation, 1, false, false);
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

void SMHGGVertex::setCoupling(Energy2 q2, tcPDPtr part2, tcPDPtr part3, tcPDPtr part1) {
  assert(part1 && part2 && part3);
  if (part1->id() != ParticleID::h0 && part2->id() != ParticleID::g &&
      part3->id() != ParticleID::g ) 
    throw HelicityConsistencyError() 
      << "SMHGGVertex::setCoupling() - The particle content of this vertex "
      << "is incorrect: " << part1->id() << " " << part2->id() << " " << part3->id() 
      << Exception::runerror;
  unsigned int Qminloop = _minloop;
  unsigned int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }
  switch (_CoefRepresentation) {
  case 1: {
    if(q2 != _q2last) {
      double g   = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    setNorm(_couplast);
    Complex loop(0.);
    for (unsigned int i = Qminloop; i <= Qmaxloop; ++i) {
      tcPDPtr qrk = getParticleData(i);
      Energy mass = (2 == massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
      loop += Af(sqr(mass)/q2);
    }
    a00(loop);
    a11(0.0);
    a12(0.0);
    a21(-loop);
    a22(0.0);
    aEp(0.0);
    break;
  }
  case 2: {
    if (q2 != _q2last) {
      double g = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = gs2*g/Constants::pi/sqrt(0.5*Constants::pi);
    }
    setNorm(_couplast);
    unsigned int delta = Qmaxloop - Qminloop + 1;
    type.resize(delta,PDT::SpinUnknown);
    masses.resize(delta,ZERO);
    for (unsigned int i = 0; i < delta; ++i) {
      tcPDPtr q = getParticleData(_minloop+i);
      type[i] = PDT::Spin1Half;
      masses[i] = (2 == massopt) ? _theSM->mass(q2,q) : q->mass();
      couplings.push_back(make_pair(masses[i]/_mw, masses[i]/_mw));
    }
    setNParticles(delta);
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
    break;
  }
  }
}

Complex SMHGGVertex::Af(double tau) const {
  return tau*(4.- W2(tau)*(1.-4.*tau));
}

Complex SMHGGVertex::W2(double lambda) const {
  double pi = Constants::pi;

  if (0.0 == lambda)     return 0.0;
  else if (lambda < 0.0) return 4.*sqr(asinh(0.5*sqrt(-1./lambda)));

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
