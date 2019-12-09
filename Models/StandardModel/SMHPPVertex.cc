// -*- C++ -*-
//
// SMHPPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHPPVertex class.
//

#include "SMHPPVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;

void SMHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _massopt << _minloop << _maxloop 
     << _CoefRepresentation;
}

void SMHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw, GeV) >> _massopt >> _minloop >> _maxloop 
     >> _CoefRepresentation;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHPPVertex,VVSLoopVertex>
describeHerwigSMHPPVertex("Herwig::SMHPPVertex", "Herwig.so");

void SMHPPVertex::Init() {
  
  static ClassDocumentation<SMHPPVertex> documentation
    ("This class implements the h0->gamma,gamma vertex.");

  static Parameter<SMHPPVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHPPVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHPPVertex,int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHPPVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHPPVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHPPVertex::_massopt, 2, false, false);
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


void SMHPPVertex::setCoupling(Energy2 q2, tcPDPtr part2,
                              tcPDPtr part3, tcPDPtr part1) {
  assert( part1->id() == ParticleID::h0 &&
	  part2->id() == ParticleID::gamma && part3->id() == ParticleID::gamma );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }
  switch (_CoefRepresentation) {
  case 1: {
    if(q2 != _q2last||_couplast==0.) {
      double g   = weakCoupling(q2);
      double e2 = sqr(electroMagneticCoupling(q2));
      _couplast = UnitRemoval::E * e2 * g / 8. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(0.);
    // quark loops
    for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
      tcPDPtr qrk = getParticleData(i);
      Energy mass = (2 == _massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
      Charge charge = qrk->charge();
      loop += Complex(3.*sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/invariant(0,0)));
    }
    // lepton loops
    int Lminloop = 3; // still fixed value
    int Lmaxloop = 3; // still fixed value
    for (int i = Lminloop; i <= Lmaxloop; ++i) {
      tcPDPtr lpt = getParticleData(9 + 2*i);
      Energy mass = (2 == _massopt) ? _theSM->mass(q2,lpt) : lpt->mass();
      Charge charge = lpt->charge();
      loop += Complex(sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/invariant(0,0)));
    }
    // W loop
    loop += Aw(sqr(_mw)/invariant(0,0));
    a00(loop);
    a11(0.0);
    a12(0.0);
    a21(-loop);
    a22(0.0);
    aEp(0.0);
    break;
  }
  case 2: {
    if(q2 != _q2last||_couplast==0.) {
      Looptools::clearcache();
      double e = electroMagneticCoupling(q2);
      _couplast = pow(e,3)/sqrt(sin2ThetaW());
      _q2last = q2;
    }
    norm(_couplast);
    // quarks
    int delta = Qmaxloop - Qminloop + 1;
    type.resize(delta,PDT::SpinUnknown);
    masses.resize(delta,ZERO);
    for (int i = 0; i < delta; ++i) {
      tcPDPtr q = getParticleData(_minloop+i);
      type[i] = PDT::Spin1Half;
      masses[i] = (2 == _massopt) ? _theSM->mass(q2,q) : q->mass();
      double copl = -masses[i]*3.*sqr(q->iCharge()/3.)/_mw/2.;
      couplings.push_back(make_pair(copl, copl));
    }
    // tau
    type.push_back(PDT::Spin1Half);
    tcPDPtr tau = getParticleData(ParticleID::tauminus);
    masses.push_back(_theSM->mass(q2,tau));
    double copl = -masses.back()*sqr(tau->iCharge()/3.)/_mw/2.;
    couplings.push_back(make_pair(copl, copl));
    // W
    type.push_back(PDT::Spin1);
    masses.push_back(_mw);
    const double mw = UnitRemoval::InvE*_mw;
    couplings.push_back(make_pair(mw,mw));
    setNParticles(delta+2);
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
    break;
  }
  }
}

Complex SMHPPVertex::Af(const double tau) const {
  return tau*(4. - W2(tau)*(1. - 4.*tau));
}

Complex SMHPPVertex::Aw(const double tau) const {
  return 0.5*(-3.*W2(tau)*tau*(4.*tau - 2.) - 12.*tau - 2.);
}

Complex SMHPPVertex::W2(double lambda) const {
  double pi = Constants::pi;

  if (0.0 == lambda) 
    return 0.0;

  if (lambda < 0.0) 
    return 4.*sqr(asinh(0.5*sqrt(-1./lambda)));

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

SMHPPVertex::SMHPPVertex() 
  :_couplast(0.),_q2last(),_mw(),_massopt(1),
   _minloop(6),_maxloop(6),_CoefRepresentation(1) {
  orderInGs(0);
  orderInGem(3);
  colourStructure(ColourStructure::SINGLET);
}


// functions for loops for testing
// namespace {

// Complex F0(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return tau*(1.-tau*ft);
// }

// Complex FHalf(double tau,double eta) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return -2.*tau*(eta+(1.-tau*eta)*ft);
// }

// Complex F1(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return 2.+3.*tau+3.*tau*(2.-tau)*ft;
// }
// }


void SMHPPVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(22,22,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  VVSLoopVertex::doinit();
//   // code to test the partial width
//   Energy mh = getParticleData(25)->mass();
//   Complex I(0.);
//   for(long ix=int(_minloop);ix<=int(_maxloop);++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == _massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double tau = sqr(2.*mt/mh);
//     I += 3.*sqr(double(qrk->iCharge())/3.)*FHalf(tau,1.);
//     cerr << "testing half " << FHalf(tau,1) << " " << Af(0.25*tau) << "\n";
//   }
//   for(long ix=15;ix<=15;++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == _massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double tau = sqr(2.*mt/mh);
//     I += sqr(double(qrk->iCharge())/3.)*FHalf(tau,1.);
//   }
//   I += F1(sqr(2.*_mw/mh));
//   Energy width = sqr(weakCoupling(sqr(mh))*sqr(electroMagneticCoupling(sqr(mh))))
//     /1024./pow(Constants::pi,5)/16.*sqr(mh/_mw)*mh*std::norm(I);
//   cerr << "testing anal " << width/GeV << "\n";
  if(loopToolsInitialized()) Looptools::ltexi();
}

