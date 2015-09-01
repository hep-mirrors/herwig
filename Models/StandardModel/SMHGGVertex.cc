// -*- C++ -*-
//
// SMHGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;  

SMHGGVertex::SMHGGVertex()
  :_couplast(0.),_q2last(ZERO),_mw(),massopt(1),_minloop(6),
   _maxloop(6),_CoefRepresentation(1) {
  orderInGs(2);
  orderInGem(1);
}

void SMHGGVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(21,21,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!_theSM) 
    throw InitException();
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  VVSLoopVertex::doinit();
  // code to test the partial width
//   Energy mh = getParticleData(25)->mass();
//   Complex I(0.);
//   for(long ix=int(_minloop);ix<=int(_maxloop);++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double lambda = sqr(mt/mh);
//     Complex fl;
//     if(lambda>=0.25) {
//       fl = -2.*sqr(asin(0.5/sqrt(lambda)));
//     }
//     else {
//       double etap = 0.5+sqrt(0.25-lambda);
//       double etam = 0.5-sqrt(0.25-lambda);
//       fl = 0.5*sqr(log(etap/etam))-0.5*sqr(Constants::pi)
// 	-Complex(0.,1.)*Constants::pi*log(etap/etam);
//     }
//     I += 3.*(2.*lambda+lambda*(4.*lambda-1)*fl);
//   }
//   Energy width = sqr(weakCoupling(sqr(mh))*sqr(strongCoupling(sqr(mh))))/36./8.*sqr(mh/_mw)*mh
//     /sqr(4.*sqr(Constants::pi))*std::norm(I)/Constants::pi;
//   cerr << "testing anal " << width/GeV << "\n";
  if(loopToolsInitialized()) Looptools::ltexi();
}

void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << massopt 
     << _minloop << _maxloop << _CoefRepresentation;
}

void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> massopt 
     >> _minloop >> _maxloop >> _CoefRepresentation;
}

ClassDescription<SMHGGVertex> SMHGGVertex::initSMHGGVertex;
// Definition of the static class description member.

void SMHGGVertex::Init() {
  
  static ClassDocumentation<SMHGGVertex> documentation
    ("This class implements the h->g,g vertex");

  static Parameter<SMHGGVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHGGVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHGGVertex,int> interfaceMaxQuarkInLoop
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
  static SwitchOption interfaceInfiniteTopMass
    (interfaceMassOption,
     "InfiniteTopMass",
     "the loop consists of an infinitely massive top quark",
     3);

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
  assert(part1->id() == ParticleID::h0 &&
	 part2->id() == ParticleID::g  && part3->id() == ParticleID::g );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }
  if(massopt==3) {
    if(q2 != _q2last) {
      double g   = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(2./3.);
    a00( loop);    a11(0.0);   a12(0.0);
    a21(-loop);    a22(0.0);   aEp(0.0);
    return;
  }
  switch (_CoefRepresentation) {
  case 1: {
    if(q2 != _q2last||_couplast==0.) {
      double g   = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(0.);
    for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
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
      Looptools::clearcache();
      _couplast = 0.25*sqr(strongCoupling(q2))*weakCoupling(q2);
      _q2last = q2;
    }
    norm(_couplast);
    int delta = Qmaxloop - Qminloop + 1;
    type.resize(delta,PDT::SpinUnknown);
    masses.resize(delta,ZERO);
    couplings.clear();
    for (int i = 0; i < delta; ++i) {
      tcPDPtr q = getParticleData(_minloop+i);
      type[i] = PDT::Spin1Half;
      masses[i] = (2 == massopt) ? _theSM->mass(q2,q) : q->mass();
      const double ratio = masses[i]/_mw;
      couplings.push_back(make_pair(ratio, ratio));
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
    // formulae from NPB297,221
  if(root < 1.) {
    ac = -sqr(asin(root));
  } 
  else {
    double ex = acosh(root);
    ac = sqr(ex) - 0.25*sqr(pi) - pi*ex*Complex(0.,1.);
  }
  return 4.*ac;
}
