// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KornerKramerCharmDecayer class.
//

#include "KornerKramerCharmDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KornerKramerCharmDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

KornerKramerCharmDecayer::~KornerKramerCharmDecayer() {}

KornerKramerCharmDecayer::KornerKramerCharmDecayer() 
{
  // default value of the fermi constant taken from PDG 2002
  _GF = 1.16639E-5/GeV2;
  // one over the number of colours 
  _oneNC=0.;
  // pseudoscalar meson decay constants
  _fpi = 131.7*MeV;
  _FK  = 160.6*MeV;
  // vector decay constants
  _frho   = 0.272;
  _fKstar = 0.238;
  // masses for the form-factors for the factorizing diagrams
  _mdcplus  = 2.42*GeV;
  _mdcminus = 2.01*GeV;
  _mscplus  = 2.54*GeV;
  _mscminus = 2.11*GeV;
  // perturbative factors
  _cplus  = 0.73;
  _cminus = 1.90;
  // factors for the non-factorizing diagrams
  _H2 = 0.119*GeV;
  _H3 =-0.011*GeV;
  // lambda_c to lambda pi+
  _incoming.push_back(4122);_outgoingB.push_back(3122);_outgoingM.push_back(211);
  _maxweight.push_back(0.0153611);
  _I1.push_back(-1./3.);_I2.push_back(-1./3.);
  _I3.push_back(-1./3.);_I4.push_back( 2./3.);_I5.push_back( 1./6.);
  _Ihat3.push_back(-1./3.);_Ihat4.push_back(2./3.);
  // lambda_c to sigma_0 pi+
  _incoming.push_back(4122);_outgoingB.push_back(3212);_outgoingM.push_back(211);
  _maxweight.push_back(0.00684877);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(-1./sqrt(3.));_I4.push_back(0.);_I5.push_back(-1./sqrt(12.));
  _Ihat3.push_back(1./sqrt(3.));_Ihat4.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ pi0
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(111);
  _maxweight.push_back(0.00686766);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(1./sqrt(3.));_I4.push_back(0.);_I5.push_back(1./sqrt(12.));
  _Ihat3.push_back(-1./sqrt(3.));_Ihat4.push_back(2./sqrt(3.));
  // lambda_c+ to sigma+ eta
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(221);
  _maxweight.push_back(0.00311807);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(-0.169101981);_I4.push_back(0.577350262);_I5.push_back(0.204124141);
  _Ihat3.push_back(0.408248281);_Ihat4.push_back(-0.816496562);
  // lambda_c+ to sigma+ eta'
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(331);
  _maxweight.push_back(0.0267591);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.985598555);_I4.push_back(-0.577350261);_I5.push_back(0.204124147);
  _Ihat3.push_back( 0.408248294);_Ihat4.push_back(-0.816496587);
  // lambda_c to p k0
  _incoming.push_back(4122);_outgoingB.push_back(2212);_outgoingM.push_back(-311);
  _maxweight.push_back(0.0444906);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(-2./sqrt(6.));_I4.push_back(2./sqrt(6.));_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c+ to xi K+
  _incoming.push_back(4122);_outgoingB.push_back(3322);_outgoingM.push_back(321);
  _maxweight.push_back(0.00533608);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(2./sqrt(6.));_I5.push_back(1./sqrt(6.));
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c to sigma+ Kbar0
  _incoming.push_back(4232);_outgoingB.push_back(3222);_outgoingM.push_back(-311);
  _maxweight.push_back(0.135814);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(-2./sqrt(6.));_Ihat4.push_back(4./sqrt(6.));
  // xi_c to xi_0 pi+
  _incoming.push_back(4232);_outgoingB.push_back(3322);_outgoingM.push_back(211);
  _maxweight.push_back(0.0745085);
  _I1.push_back(1./sqrt(6.));_I2.push_back(1./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(2./sqrt(6.));_Ihat4.push_back(-4./sqrt(6.));
  // xi_c to lambda kbar0
  _incoming.push_back(4132);_outgoingB.push_back(3122);_outgoingM.push_back(-311);
  _maxweight.push_back(0.0025821);
  _I1.push_back(1./6.);_I2.push_back(1./6.);
  _I3.push_back(2./3.);_I4.push_back(-1./3.);_I5.push_back(1./6.);
  _Ihat3.push_back(-1./3.);_Ihat4.push_back(2./3.);
  // xi_c to sigma k0
  _incoming.push_back(4132);_outgoingB.push_back(3212);_outgoingM.push_back(-311);
  _maxweight.push_back(0.024961);
  _I1.push_back(1./sqrt(12.));_I2.push_back(1./sqrt(12.));
  _I3.push_back(0.);_I4.push_back(-2./sqrt(12.));_I5.push_back(-1./sqrt(12.));
  _Ihat3.push_back(2./sqrt(12.));_Ihat4.push_back(-4./sqrt(12.));
  // xi_c to sigma k-
  _incoming.push_back(4132);_outgoingB.push_back(3222);_outgoingM.push_back(-321);
  _maxweight.push_back(0.00250939);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(2./sqrt(6.));_I5.push_back(1./sqrt(6.));
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi0 pi0
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(111);
  _maxweight.push_back(0.000739771);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(2./sqrt(12.));_I4.push_back(-2./sqrt(12.));_I5.push_back(0.);
  _Ihat3.push_back(-2./sqrt(12.));_Ihat4.push_back(4./sqrt(12.));
  // xi_c0 to xi0 eta
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(221);
  _maxweight.push_back(0.00490303);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(-0.169101981);_I4.push_back(-0.408248281);_I5.push_back(-0.288675131);
  _Ihat3.push_back(0.408248281);_Ihat4.push_back(-0.816496562);
  // xi_c0 to xi0 eta'
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(331);
  _maxweight.push_back(0.0178693);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.985598555);_I4.push_back(-0.408248294);_I5.push_back(0.288675131);
  _Ihat3.push_back(0.408248294);_Ihat4.push_back(-0.816496587);
  // xi_c0 to xi- pi+
  _incoming.push_back(4132);_outgoingB.push_back(3312);_outgoingM.push_back(211);
  _maxweight.push_back(0.0220038);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(-2./sqrt(6.));_I4.push_back(2./sqrt(6.));_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c to xi0 K0
  _incoming.push_back(4332);_outgoingB.push_back(3322);_outgoingM.push_back(-311);
  _maxweight.push_back(0.0250268);
  _I1.push_back(1.);_I2.push_back(-1.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(2.);_Ihat4.push_back(0.);
  // lambda_c to lambda rho+
  _incoming.push_back(4122);_outgoingB.push_back(3122);_outgoingM.push_back(213);
  _maxweight.push_back(0.467008);
  _I1.push_back(-1./3.);_I2.push_back(-1./3.);
  _I3.push_back(-1./3.);_I4.push_back( 2./3.);_I5.push_back( 1./6.);
  _Ihat3.push_back(-1./3.);_Ihat4.push_back(2./3.);
  // lambda_c to sigma_0 rho+
  _incoming.push_back(4122);_outgoingB.push_back(3212);_outgoingM.push_back(213);
  _maxweight.push_back(0.0702498);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(-1./sqrt(3.));_I4.push_back(0.);_I5.push_back(-1./sqrt(12.));
  _Ihat3.push_back(1./sqrt(3.));_Ihat4.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ rho
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(113);
  _maxweight.push_back(0.0697802);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(1./sqrt(3.));_I4.push_back(0.);
  _I5.push_back(1./sqrt(12.));
  _Ihat3.push_back(-1./sqrt(3.));_Ihat4.push_back(2./sqrt(3.));
  // lambda_c to sigma+ omega
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(223);
  _maxweight.push_back(0.252355);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.57735033);_I4.push_back(0.);_I5.push_back( 0.288675151);
  _Ihat3.push_back(0.577350303);_Ihat4.push_back( -1.15470061);
  // lambda_c to simga+ phi
  _incoming.push_back(4122);_outgoingB.push_back(3222);_outgoingM.push_back(333);
  _maxweight.push_back(0.0063064);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.816496614);_I4.push_back(-0.816496624);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to p k*0bar
  _incoming.push_back(4122);_outgoingB.push_back(2212);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0996461);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(-2./sqrt(6.));_I4.push_back(2./sqrt(6.));_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c+ to xi K*+
  _incoming.push_back(4122);_outgoingB.push_back(3322);_outgoingM.push_back(323);
  _maxweight.push_back(0.00413946);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(2./sqrt(6.));_I5.push_back(1./sqrt(6.));
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c to sigma+ K*bar0
  _incoming.push_back(4232);_outgoingB.push_back(3222);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0583248);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(-2./sqrt(6.));_Ihat4.push_back(4./sqrt(6.));
  // xi_c to xi_0 rho+
  _incoming.push_back(4232);_outgoingB.push_back(3322);_outgoingM.push_back(213);
  _maxweight.push_back(2.18754);
  _I1.push_back(1./sqrt(6.));_I2.push_back(1./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(2./sqrt(6.));_Ihat4.push_back(-4./sqrt(6.));
  // xi_c to lambda k*bar0
  _incoming.push_back(4132);_outgoingB.push_back(3122);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0373243);
  _I1.push_back(1./6.);_I2.push_back(1./6.);
  _I3.push_back(2./3.);_I4.push_back(-1./3.);_I5.push_back(1./6.);
  _Ihat3.push_back(-1./3.);_Ihat4.push_back(2./3.);
  // xi_c to sigma k*0
  _incoming.push_back(4132);_outgoingB.push_back(3212);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0135937);
  _I1.push_back(1./sqrt(12.));_I2.push_back(1./sqrt(12.));
  _I3.push_back(0.);_I4.push_back(-2./sqrt(12.));_I5.push_back(-1./sqrt(12.));
  _Ihat3.push_back(2./sqrt(12.));_Ihat4.push_back(-4./sqrt(12.));
  // xi_c to sigma k*-
  _incoming.push_back(4132);_outgoingB.push_back(3222);_outgoingM.push_back(-323);
  _maxweight.push_back(0.00877507);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(2./sqrt(6.));_I5.push_back(1./sqrt(6.));
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi0 rho0
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(113);
  _maxweight.push_back(0.0364247);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(2./sqrt(12.));_I4.push_back(-2./sqrt(12.));_I5.push_back(0.);
  _Ihat3.push_back(-2./sqrt(12.));_Ihat4.push_back(4./sqrt(12.));
  // xi_c0 to xi0 omega
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(223);
  _maxweight.push_back(0.134395);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back( 0.57735033);_I4.push_back(-0.577350303);_I5.push_back(0.);
  _Ihat3.push_back(0.577350303);_Ihat4.push_back(-1.15470061);
  // xi_c0 to xi0 phi
  _incoming.push_back(4132);_outgoingB.push_back(3322);_outgoingM.push_back(333);
  _maxweight.push_back(0.00676745);
  _I1.push_back(0.);_I2.push_back(0.);
  _I3.push_back(0.816496614);_I4.push_back(0.);_I5.push_back(0.408248312);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi- rho+
  _incoming.push_back(4132);_outgoingB.push_back(3312);_outgoingM.push_back(213);
  _maxweight.push_back(0.309733);
  _I1.push_back(-1./sqrt(6.));_I2.push_back(-1./sqrt(6.));
  _I3.push_back(-2./sqrt(6.));_I4.push_back(2./sqrt(6.));_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c to xi0 K*0
  _incoming.push_back(4332);_outgoingB.push_back(3322);_outgoingM.push_back(-313);
  _maxweight.push_back(0.019967);
  _I1.push_back(1.);_I2.push_back(-1.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(2.);_Ihat4.push_back(0.);
  // lambda_c to sigma*0 pi+
  _incoming.push_back(4122);_outgoingB.push_back(3214);_outgoingM.push_back(211);
  _maxweight.push_back(0.0254397);
  _I1.push_back(0.);_I2.push_back(1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ pi0
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(111);
  _maxweight.push_back(0.0256531);
  _I1.push_back(0.);_I2.push_back(1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ eta
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(221);
  _maxweight.push_back(0.0299659);
  _I1.push_back(0.);_I2.push_back( 0.569036);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ eta'
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(331);
  _maxweight.push_back(4.4186e-07);
  _I1.push_back(0.);_I2.push_back(-0.097631);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to delta+ kbar0
  _incoming.push_back(4122);_outgoingB.push_back(2214);_outgoingM.push_back(-311);
  _maxweight.push_back(0.0215801);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to delta++ k-
  _incoming.push_back(4122);_outgoingB.push_back(2224);_outgoingM.push_back(-321);
  _maxweight.push_back(0.0649794);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to xi*0 k+
  _incoming.push_back(4122);_outgoingB.push_back(3324);_outgoingM.push_back(321);
  _maxweight.push_back(0.0198121);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  _incoming.push_back(4132);_outgoingB.push_back(3214);_outgoingM.push_back(-311);
  _maxweight.push_back(0.0118739);
  _I1.push_back(0.);_I2.push_back(-1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to sigma*+ Kbar-
  _incoming.push_back(4132);_outgoingB.push_back(3224);_outgoingM.push_back(-321);
  _maxweight.push_back(0.0240231);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 pi0
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(111);
  _maxweight.push_back(0.0173433);
  _I1.push_back(0.);_I2.push_back(-1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 eta
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(221);
  _maxweight.push_back(0.000939099);
  _I1.push_back(0.);_I2.push_back(0.097631);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 eta'
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(331);
  _maxweight.push_back(1.31513e-05);
  _I1.push_back(0.);_I2.push_back(-0.569036);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*- pi+
  _incoming.push_back(4132);_outgoingB.push_back(3314);_outgoingM.push_back(211);
  _maxweight.push_back(0.0338039);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to omega- K+
  _incoming.push_back(4132);_outgoingB.push_back(3334);_outgoingM.push_back(321);
  _maxweight.push_back(0.00715116);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c0 to xi*0 kbar0
  _incoming.push_back(4332);_outgoingB.push_back(3324);_outgoingM.push_back(-311);
  _maxweight.push_back(0.00850994);
  _I1.push_back(-1./sqrt(3.));_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c0 to omega pi+
  _incoming.push_back(4332);_outgoingB.push_back(3334);_outgoingM.push_back(211);
  _maxweight.push_back(0.012698);
  _I1.push_back(-1./sqrt(3.));_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*0 rho+
  _incoming.push_back(4122);_outgoingB.push_back(3214);_outgoingM.push_back(213);
  _maxweight.push_back(0.0478329);
  _I1.push_back(0.);_I2.push_back(1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ rho0
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(113);
  _maxweight.push_back(0.0476973);
  _I1.push_back(0.);_I2.push_back(1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ omega
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(223);
  _maxweight.push_back(0.0260017);
  _I1.push_back(0.);_I2.push_back(1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to sigma*+ phi
  _incoming.push_back(4122);_outgoingB.push_back(3224);_outgoingM.push_back(333);
  _maxweight.push_back(3.55986e-07);
  _I1.push_back(0.);_I2.push_back(-sqrt(2.)/3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to delta+ kbar0
  _incoming.push_back(4122);_outgoingB.push_back(2214);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0988993);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to delta++ k-
  _incoming.push_back(4122);_outgoingB.push_back(2224);_outgoingM.push_back(-323);
  _maxweight.push_back(0.303435);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(6.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // lambda_c to xi*0 k+
  _incoming.push_back(4122);_outgoingB.push_back(3324);_outgoingM.push_back(323);
  _maxweight.push_back(9.70866e-05);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  _incoming.push_back(4132);_outgoingB.push_back(3214);_outgoingM.push_back(-313);
  _maxweight.push_back(0.0281627);
  _I1.push_back(0.);_I2.push_back(-1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to sigma*- K-
  _incoming.push_back(4132);_outgoingB.push_back(3224);_outgoingM.push_back(-323);
  _maxweight.push_back(0.0575937);
  _I1.push_back(0.);_I2.push_back(-2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 pi0
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(113);
  _maxweight.push_back(0.0480608);
  _I1.push_back(0.);_I2.push_back(-1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 omega
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(223);
  _maxweight.push_back(0.0223044);
  _I1.push_back(0.);_I2.push_back(-1./3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*0 phi
  _incoming.push_back(4132);_outgoingB.push_back(3324);_outgoingM.push_back(333);
  _maxweight.push_back(3.187e-09);
  _I1.push_back(0.);_I2.push_back(-sqrt(2.)/3.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to xi*- rho+
  _incoming.push_back(4132);_outgoingB.push_back(3314);_outgoingM.push_back(213);
  _maxweight.push_back(0.0958052);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // xi_c0 to omega- K*+
  _incoming.push_back(4132);_outgoingB.push_back(3334);_outgoingM.push_back(323);
  _maxweight.push_back(0.000205851);
  _I1.push_back(0.);_I2.push_back(2./sqrt(18.));
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c0 to xi*0 k*bar0
  _incoming.push_back(4332);_outgoingB.push_back(3324);_outgoingM.push_back(-313);
  _maxweight.push_back(0.050904);
  _I1.push_back(-1./sqrt(3.));_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // omega_c0 to omega rho+
  _incoming.push_back(4332);_outgoingB.push_back(3334);_outgoingM.push_back(213);
  _maxweight.push_back(0.0814021);
  _I1.push_back(-1./sqrt(3.));_I2.push_back(0.);
  _I3.push_back(0.);_I4.push_back(0.);_I5.push_back(0.);
  _Ihat3.push_back(0.);_Ihat4.push_back(0.);
  // initial size of the vectors
  _initsize=_incoming.size();
}

void KornerKramerCharmDecayer::doinit() throw(InitException) {
  Baryon1MesonDecayerBase::doinit();
  // check the vectors have the same size
  unsigned int isize=_incoming.size();
  if(isize!=_I1.size()||isize!=_I2.size()||isize!=_I3.size()||isize!=_I4.size()||
     isize!=_I5.size()||isize!=_Ihat3.size()||isize!=_Ihat4.size()||
     isize!=_incoming.size()||isize!=_outgoingB.size()||isize!=_outgoingM.size()||
     isize!=_maxweight.size())
    {throw InitException() << "Inconsistent parameters in"
			   << " KornerKramerCharmDecayer::doinit()" 
			   << Exception::abortnow;}
  // compute the various coefficients
  Energy m1,m2,m3,fmes; 
  Energy2 P1P2,Qplus,Qminus,gmes;
  double Fnonfact,H2,H3,A,B,A2,B2,A3,B3,Ffact[2];
  int mspin,bspin;
  double chi,
    chiplus( 0.5*(_cplus*(1.+_oneNC)+_cminus*(1.-_oneNC))),
    chiminus(0.5*(_cplus*(1.+_oneNC)-_cminus*(1.-_oneNC))); 
  InvEnergy2 pre(_GF*sqrt(SM().CKM(1,1)*SM().CKM(0,0)/2.)),mform2[2];

  // testing only
  pre = _GF*0.974/sqrt(2.);
  vector<double> wgt(1,1.);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  unsigned int iy;
  for(unsigned int ix=0;ix<isize;++ix)
    {
      // get the mass of the particles
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoingB[ix]);
      extpart[2]=getParticleData(_outgoingM[ix]);
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,_maxweight[ix],wgt);
      m1=extpart[0]->mass();
      m2=extpart[1]->mass();
      m3=extpart[2]->mass();
      // dot product P1P2
      P1P2 = 0.5*(m1*m1+m2*m2-m3*m3);
      // formfactor for the non-factorizating diagrams and coefficients
      Fnonfact = 1./(1.-m3*m3/_mscminus/_mscminus)*
	(1.-(m1-m2)*(m1-m2)/_mscminus/_mscminus);
      Fnonfact *= Fnonfact*Fnonfact;
      H2 = _H2*Fnonfact;
      H3 = _H3*Fnonfact;
      if(abs(_outgoingM[ix])==211)
	{
	  fmes=_fpi;
	  chi=chiplus;
	  mform2[0]=1./(_mscminus*_mscminus);
	  mform2[1]=1./(_mscplus *_mscplus);
	}
      else if(abs(_outgoingM[ix])==311)
	{
	  fmes=_FK;
	  chi=chiminus;
	  mform2[0]=1./(_mdcminus*_mdcminus);
	  mform2[1]=1./(_mdcplus *_mdcplus);
	}
      else if(abs(_outgoingM[ix])==213)
	{
	  gmes = _frho*m3*m3;
	  chi=chiplus;
	  mform2[0]=1./(_mscminus*_mscminus);
	  mform2[1]=1./(_mscplus *_mscplus);
	}
      else if(abs(_outgoingM[ix])==313)
	{
	  gmes = _fKstar*m3*m3;
	  chi=chiminus;
	  mform2[0]=1./(_mdcminus*_mdcminus);
	  mform2[1]=1./(_mdcplus *_mdcplus);
	}
      else
	{
	  fmes=0.;chi=0.;
	  mform2[0]=1./(_mscminus*_mscminus);
	  mform2[1]=1./(_mscplus *_mscplus);
	  gmes=0.;
	}
      // form factor for the factorising diagrams
      for(iy=0;iy<2;++iy)
	{
	  Ffact[iy]  = 1./(1.-m3*m3*mform2[iy])*(1.-(m1-m2)*(m1-m2)*mform2[iy]);
	  Ffact[iy] *= Ffact[iy]*Ffact[iy];
	}
      // invariants
      Qplus  = (m1+m2)*(m1+m2)-m3*m3;
      Qminus = (m1-m2)*(m1-m2)-m3*m3;
      // decide which type of decay
      mspin=getParticleData(_outgoingM[ix])->iSpin();
      bspin=getParticleData(_outgoingB[ix])->iSpin();
      if(bspin==2)
	{
	  if(mspin==1)
	    {
	      // non-factorizing piece of the diagrams
	      A = -0.75*H2/m1/m2*_cminus*( (m1*P1P2-m1*m2*(m2+m3))*_I3[ix]
					  -(m2*P1P2-m1*m2*(m1+m3))*_Ihat3[ix]);
	      B = 0.25/m1/m2*_cminus*( 0.5*H2*Qplus*( m1*(   _I3[ix]+2.*    _I4[ix])
						     +m2*(_Ihat3[ix]+2.*_Ihat4[ix]))
				       +H3*m1*m2*(m1+m2+m3)*12.*_I5[ix]);
	      // the factorizing piece
	      A +=0.25/m1/m2*chi*fmes*Qplus*(m1-m2)*(2.*_I1[ix]+_I2[ix])*Ffact[0];
	      B +=0.25/m1/m2*chi*fmes*Qplus*(m1+m2)/3.*(4.*_I1[ix]+5.*_I2[ix])*Ffact[1];
	      // add to vectors
	      _A1.push_back(A*pre);_B1.push_back(B*pre);
	      _A2.push_back(0.);_B2.push_back(0.);
	      _A3.push_back(0.);_B3.push_back(0.);
	    }
	  else if(mspin==3)
	    {
	      // non-factorizing terms
	      A  = -0.25*H2/m1/m2*_cminus*
		( (m1*P1P2-m1*m2*(m2+m3))*(   _I3[ix]+2.*   _I4[ix])
		 +(m2*P1P2-m1*m2*(m1+m3))*(_Ihat3[ix]+2.*_Ihat4[ix]));
	      A2 = -0.25*H2/m1/m2*_cminus*(m1+m2+m3)*(+m1*(   _I3[ix]+2.*   _I4[ix])
						      -m2*(_Ihat3[ix]+2.*_Ihat4[ix]));
	      B  = +0.25/m1/m2*_cminus*(0.5*H2*Qplus*(+m1*(   _I3[ix]+2.*   _I4[ix])
						      +m2*(_Ihat3[ix]+2.*_Ihat4[ix]))
					+H3*m1*m2*(m1+m2+m3)*12.*_I5[ix]);
	      B2 = -0.25/m1/m2*_cminus*
		(+H2*(+m1*(m1+m2)*(   _I3[ix]+2.*   _I4[ix])-3.*m1*m3*   _I3[ix]
		      +m2*(m1+m2)*(_Ihat3[ix]+2.*_Ihat4[ix])-3.*m2*m3*_Ihat3[ix])
		 +H3*m1*m2*24.*_I5[ix]);
	      // the factorizing piece
	      A  += -0.25/m1/m2*chi*gmes*Qplus/3.*(4.*_I1[ix]+5.*_I2[ix])*Ffact[1];
	      B  +=  0.25/m1/m2*chi*gmes*Qplus/3.*(4.*_I1[ix]+5.*_I2[ix])*Ffact[0];
	      B2 +=  0.25/m1/m2*chi*gmes*(m1+m2)*4./3.*(_I1[ix]-_I2[ix])*Ffact[0];
	      // add to vectors
	      _A1.push_back(A*pre);_B1.push_back(B*pre);
	      _A2.push_back(A2*pre);_B2.push_back(B2*pre);
	      _A3.push_back(0.);_B3.push_back(0.);
	    }
	  else
	    {throw InitException() << "Invalid outgoing meson spin in"
				   << " KornerKramerCharmDecayer::doinit()" 
				   << Exception::abortnow;}
	}
      else if(bspin==4)
	{
	  if(mspin==1)
	    {
	      // first the non-factorizing piece
	      B  = -1.5*_cminus*H2*_I2[ix];
	      // then the factorizing piece
	      B  += chi*fmes*(1.+m2/m1)*_I1[ix]*Ffact[1];
	      // add to vectors
	      // make the coupling dimensionless
	      _A1.push_back(0.);_B1.push_back(0.);
	      _A2.push_back(0.);_B2.push_back(B*pre);
	      _A3.push_back(0.);_B3.push_back(0.);
	    }
	  else if(mspin==3)
	    {
	      double norm(0.75/m1/m2*_cminus*H2*_I2[ix]);
	      // first the non-factorizing piece
	      A  = -norm*m1*(P1P2-m2*(m2+m3))*2.;
	      A2 = 0.;
	      A3 =  norm*m1*2.;
	      B  = -norm*m1*Qplus;
	      B2 = -norm*m1*m2*2.;
	      B3 =  norm*m1*2.;
	      // then the factorizing piece
	      norm = 0.5/m1/m2*chi*gmes*_I1[ix];
	      A  += norm*Qplus*Ffact[1];
	      A2 += 0.;
	      A3 +=-norm*2.*Ffact[1];
	      B  += norm*Qplus*Ffact[0];
	      B2 += norm*m2*2.*Ffact[0];
	      B3 +=-norm*2.*Ffact[0];
	      // add to vectors
	      _A1.push_back( A*pre);_B1.push_back( B*pre);
	      _A2.push_back(2.*A2*pre);_B2.push_back(2.*B2*pre);
	      _A3.push_back(A3*pre);_B3.push_back(B3*pre);
	    }
	  else
	    {throw InitException() << "Invalid outgoing meson spin in"
				   << " KornerKramerCharmDecayer::doinit()" 
				   << Exception::abortnow;}
	}
      else
	{throw InitException() << "Invalid outgoing baryon spin in"
			       << " KornerKramerCharmDecayer::doinit()" 
			       << Exception::abortnow;}
    }
}


bool KornerKramerCharmDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do
    {
      if(id0==_incoming[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incoming[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){allowed=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector KornerKramerCharmDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1),id(parent.id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);bool cc;
  do 
    {
      if(id==_incoming[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incoming[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;cc=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}

void KornerKramerCharmDecayer::persistentOutput(PersistentOStream & os) const {
  os << _GF << _oneNC << _fpi << _FK << _frho << _fKstar << _mdcplus << _mdcminus 
     << _mscplus << _mscminus << _cplus << _cminus << _H2 << _H3 << _I1 << _I2 
     << _I3 << _I4 << _I5 << _Ihat3 << _Ihat4 << _incoming << _outgoingB 
     << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3 << _initsize;
}

void KornerKramerCharmDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _GF >> _oneNC >> _fpi >> _FK >> _frho >> _fKstar >> _mdcplus >> _mdcminus 
     >> _mscplus >> _mscminus >> _cplus >> _cminus >> _H2 >> _H3 >> _I1 >> _I2 
     >> _I3 >> _I4 >> _I5 >> _Ihat3 >> _Ihat4 >> _incoming >> _outgoingB 
     >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3 >> _initsize;
}

ClassDescription<KornerKramerCharmDecayer> 
KornerKramerCharmDecayer::initKornerKramerCharmDecayer;
// Definition of the static class description member.

void KornerKramerCharmDecayer::Init() {

  static ClassDocumentation<KornerKramerCharmDecayer> documentation
    ("The \\classname{KornerKramerCharmDecayer} class implements the"
     " non-leptonic weak decay of charm baryons using the results of Z.Phys.C55,659.");

  static Parameter<KornerKramerCharmDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &KornerKramerCharmDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);

  static Parameter<KornerKramerCharmDecayer,double> interfaceOneOverNc
    ("OneOverNc",
     "One divided by the number of colours, the default is to take N_c to infinity.",
     &KornerKramerCharmDecayer::_oneNC, 0.0, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFpi
    ("Fpi",
     "The decay constant for the pi meson.",
     &KornerKramerCharmDecayer::_fpi, MeV, 131.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFK
    ("FK",
     "The decay constant for the K meson",
     &KornerKramerCharmDecayer::_FK, MeV, 160.6*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceFrho
    ("Frho",
     "The decay constant for the rho meson",
     &KornerKramerCharmDecayer::_frho, 0.272, 0.0, 1.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfacefKstar
    ("fKstar",
     "The decay constant for the K star meson",
     &KornerKramerCharmDecayer::_fKstar, 0.238, 0.0, 1.0,
     false, false, false);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcplus
    ("Mdcplus",
     "The mass of the 1+ dc meson for the form-factors",
     &KornerKramerCharmDecayer::_mdcplus, GeV, 2.42*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscplus
    ("Mscplus",
     "The mass of the 1+ sc meson for the form-factors",
     &KornerKramerCharmDecayer::_mscplus, GeV, 2.54*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcminus
    ("Mdcminus",
     "The mass of the 1- dc meson for the form-factors",
     &KornerKramerCharmDecayer::_mdcminus, GeV, 2.01*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscminus
    ("Mscminus",
     "The mass of the 1- sc meson for the form-factors",
     &KornerKramerCharmDecayer::_mscminus, GeV, 2.11*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCplus
    ("Cplus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::_cplus, 0.73, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCminus
    ("Cminus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::_cminus, 1.90, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH2
    ("H2",
     "The H2 parameter",
     &KornerKramerCharmDecayer::_H2, GeV, 0.119*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH3
    ("H3",
     "The H3 parameter",
     &KornerKramerCharmDecayer::_H3, GeV,-0.011*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI1
    ("I1",
     "The I_1 invariant for the decay modes",
     &KornerKramerCharmDecayer::_I1, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI2
    ("I2",
     "The I_2 invariant for the decay modes",
     &KornerKramerCharmDecayer::_I2, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI3
    ("I3",
     "The I_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::_I3, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI4
    ("I4",
     "The I_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::_I4, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI5
    ("I5",
     "The I_5 invariant for the decay modes",
     &KornerKramerCharmDecayer::_I5, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat3
    ("Ihat3",
     "The Ihat_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::_Ihat3, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat4
    ("Ihat4",
     "The Ihat_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::_Ihat4, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &KornerKramerCharmDecayer::_incoming, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingB
    ("OutgoingB",
     "The PDG code of the outgoing baryon",
     &KornerKramerCharmDecayer::_outgoingB, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingM
    ("OutgoingM",
     "The PDG code of the outgoing meson",
     &KornerKramerCharmDecayer::_outgoingM, -1, 0, -1000000, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &KornerKramerCharmDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void KornerKramerCharmDecayer::
 halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
			Complex& A,Complex&B) const
{A =_A1[imode];B=_B1[imode];}

// couplings for spin-1/2 to spin-1/2 spin-1
void KornerKramerCharmDecayer::halfHalfVectorCoupling(int imode,
						      Energy m0,Energy m1,Energy m2,
						      Complex& A1,Complex& A2,
						      Complex& B1,Complex& B2) const
{
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A1 = _B1[imode]; B1 = _A1[imode];
  A2 = _B2[imode]*(m0+m1); B2 = _A2[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-0
void KornerKramerCharmDecayer::halfThreeHalfScalarCoupling(int imode,
							   Energy m0,Energy m1,Energy m2,
							   Complex&  A, Complex& B) const
{
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A = _B2[imode]*(m0+m1);B=_A2[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-1
void KornerKramerCharmDecayer::
halfThreeHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
			    Complex& A1, Complex& A2, Complex& A3,
			    Complex& B1, Complex& B2, Complex& B3) const
{
  A1=_A1[imode];A2=_A2[imode]*(m0+m1);A3=_A3[imode]*(m0+m1)*(m0+m1);
  B1=_B1[imode];B2=_B2[imode]*(m0+m1);B3=_B3[imode]*(m0+m1)*(m0+m1);
}
 
void KornerKramerCharmDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  output << "set " << fullName() << ":GFermi " << _GF*GeV2 << "\n";
  output << "set " << fullName() << ":OneOverNc " <<  _oneNC<< "\n";
  output << "set " << fullName() << ":Fpi " << _fpi/MeV << "\n";
  output << "set " << fullName() << ":FK " <<_FK/MeV  << "\n";
  output << "set " << fullName() << ":Frho " << _frho << "\n";
  output << "set " << fullName() << ":fKstar " << _fKstar << "\n";
  output << "set " << fullName() << ":Mdcplus " << _mdcplus/GeV << "\n";
  output << "set " << fullName() << ":Mscplus " << _mscplus/GeV << "\n";
  output << "set " << fullName() << ":Mdcminus " << _mdcminus/GeV << "\n";
  output << "set " << fullName() << ":Mscminus " << _mscminus/GeV << "\n";
  output << "set " << fullName() << ":Cplus " << _cplus << "\n";
  output << "set " << fullName() << ":Cminus " << _cminus << "\n";
  output << "set " << fullName() << ":H2 " << _H2/GeV << "\n";
  output << "set " << fullName() << ":H3 " << _H3/GeV << "\n";
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":I1 " 
		 << ix << " " << _I1[ix] << "\n";
	  output << "set " << fullName() << ":I2 " 
		 << ix << " " << _I2[ix] << "\n";
	  output << "set " << fullName() << ":I3 " 
		 << ix << " " << _I3[ix] << "\n";
	  output << "set " << fullName() << ":I4 " 
		 << ix << " " << _I4[ix] << "\n";
	  output << "set " << fullName() << ":I5 " 
		 << ix << " " << _I5[ix] << "\n";
	  output << "set " << fullName() << ":Ihat3 " 
		 << ix << " " << _Ihat3[ix] << "\n";
	  output << "set " << fullName() << ":Ihat4 " 
		 << ix << " " << _Ihat4[ix] << "\n";
	  output << "set " << fullName() << ":Incoming " 
		 << ix << " " << _incoming[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingB " 
		 << ix << " " << _outgoingB[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingM " 
		 << ix << " " << _outgoingM[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " 
		 << ix << " " << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":I1 " 
		 << ix << " " << _I1[ix] << "\n";
	  output << "insert " << fullName() << ":I2 " 
		 << ix << " " << _I2[ix] << "\n";
	  output << "insert " << fullName() << ":I3 " 
		 << ix << " " << _I3[ix] << "\n";
	  output << "insert " << fullName() << ":I4 " 
		 << ix << " " << _I4[ix] << "\n";
	  output << "insert " << fullName() << ":I5 " 
		 << ix << " " << _I5[ix] << "\n";
	  output << "insert " << fullName() << ":Ihat3 " 
		 << ix << " " << _Ihat3[ix] << "\n";
	  output << "insert " << fullName() << ":Ihat4 " 
		 << ix << " " << _Ihat4[ix] << "\n";
	  output << "insert " << fullName() << ":Incoming " 
		 << ix << " " << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingB " 
		 << ix << " " << _outgoingB[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingM " 
		 << ix << " " << _outgoingM[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " 
		 << ix << " " << _maxweight[ix] << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
