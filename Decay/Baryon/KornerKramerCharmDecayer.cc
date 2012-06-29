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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

KornerKramerCharmDecayer::KornerKramerCharmDecayer() {
  // one over the number of colours 
  oneNC_=0.;
  // pseudoscalar meson decay constants
  fpi_ = 131.7*MeV;
  fk_  = 160.6*MeV;
  // vector decay constants
  frho_   = 0.272;
  fKstar_ = 0.238;
  // masses for the form-factors for the factorizing diagrams
  mdcplus_  = 2.42*GeV;
  mdcminus_ = 2.01*GeV;
  mscplus_  = 2.54*GeV;
  mscminus_ = 2.11*GeV;
  // perturbative factors
  cplus_  = 0.73;
  cminus_ = 1.90;
  // factors for the non-factorizing diagrams
  H2_ = 0.119*GeV;
  H3_ =-0.011*GeV;
  // lambda_c to lambda pi+
  incoming_.push_back(4122);outgoingB_.push_back(3122);outgoingM_.push_back(211);
  maxweight_.push_back(0.0153611);
  I1_.push_back(-1./3.);I2_.push_back(-1./3.);
  I3_.push_back(-1./3.);I4_.push_back( 2./3.);I5_.push_back( 1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // lambda_c to sigma_0 pi+
  incoming_.push_back(4122);outgoingB_.push_back(3212);outgoingM_.push_back(211);
  maxweight_.push_back(0.00684877);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-1./sqrt(3.));I4_.push_back(0.);I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(1./sqrt(3.));Ihat4_.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ pi0
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(111);
  maxweight_.push_back(0.00686766);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(1./sqrt(3.));I4_.push_back(0.);I5_.push_back(1./sqrt(12.));
  Ihat3_.push_back(-1./sqrt(3.));Ihat4_.push_back(2./sqrt(3.));
  // lambda_c+ to sigma+ eta
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(221);
  maxweight_.push_back(0.00311807);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-0.169101981);I4_.push_back(0.577350262);I5_.push_back(0.204124141);
  Ihat3_.push_back(0.408248281);Ihat4_.push_back(-0.816496562);
  // lambda_c+ to sigma+ eta'
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(331);
  maxweight_.push_back(0.0267591);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.985598555);I4_.push_back(-0.577350261);I5_.push_back(0.204124147);
  Ihat3_.push_back( 0.408248294);Ihat4_.push_back(-0.816496587);
  // lambda_c to p k0
  incoming_.push_back(4122);outgoingB_.push_back(2212);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0444906);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c+ to xi K+
  incoming_.push_back(4122);outgoingB_.push_back(3322);outgoingM_.push_back(321);
  maxweight_.push_back(0.00533608);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c to sigma+ Kbar0
  incoming_.push_back(4232);outgoingB_.push_back(3222);outgoingM_.push_back(-311);
  maxweight_.push_back(0.135814);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(6.));Ihat4_.push_back(4./sqrt(6.));
  // xi_c to xi_0 pi+
  incoming_.push_back(4232);outgoingB_.push_back(3322);outgoingM_.push_back(211);
  maxweight_.push_back(0.0745085);
  I1_.push_back(1./sqrt(6.));I2_.push_back(1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2./sqrt(6.));Ihat4_.push_back(-4./sqrt(6.));
  // xi_c to lambda kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3122);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0025821);
  I1_.push_back(1./6.);I2_.push_back(1./6.);
  I3_.push_back(2./3.);I4_.push_back(-1./3.);I5_.push_back(1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // xi_c to sigma k0
  incoming_.push_back(4132);outgoingB_.push_back(3212);outgoingM_.push_back(-311);
  maxweight_.push_back(0.024961);
  I1_.push_back(1./sqrt(12.));I2_.push_back(1./sqrt(12.));
  I3_.push_back(0.);I4_.push_back(-2./sqrt(12.));I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(2./sqrt(12.));Ihat4_.push_back(-4./sqrt(12.));
  // xi_c to sigma k-
  incoming_.push_back(4132);outgoingB_.push_back(3222);outgoingM_.push_back(-321);
  maxweight_.push_back(0.00250939);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(111);
  maxweight_.push_back(0.000739771);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(2./sqrt(12.));I4_.push_back(-2./sqrt(12.));I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(12.));Ihat4_.push_back(4./sqrt(12.));
  // xi_c0 to xi0 eta
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(221);
  maxweight_.push_back(0.00490303);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-0.169101981);I4_.push_back(-0.408248281);I5_.push_back(-0.288675131);
  Ihat3_.push_back(0.408248281);Ihat4_.push_back(-0.816496562);
  // xi_c0 to xi0 eta'
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(331);
  maxweight_.push_back(0.0178693);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.985598555);I4_.push_back(-0.408248294);I5_.push_back(0.288675131);
  Ihat3_.push_back(0.408248294);Ihat4_.push_back(-0.816496587);
  // xi_c0 to xi- pi+
  incoming_.push_back(4132);outgoingB_.push_back(3312);outgoingM_.push_back(211);
  maxweight_.push_back(0.0220038);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c to xi0 K0
  incoming_.push_back(4332);outgoingB_.push_back(3322);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0250268);
  I1_.push_back(1.);I2_.push_back(-1.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2.);Ihat4_.push_back(0.);
  // lambda_c to lambda rho+
  incoming_.push_back(4122);outgoingB_.push_back(3122);outgoingM_.push_back(213);
  maxweight_.push_back(0.467008);
  I1_.push_back(-1./3.);I2_.push_back(-1./3.);
  I3_.push_back(-1./3.);I4_.push_back( 2./3.);I5_.push_back( 1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // lambda_c to sigma_0 rho+
  incoming_.push_back(4122);outgoingB_.push_back(3212);outgoingM_.push_back(213);
  maxweight_.push_back(0.0702498);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-1./sqrt(3.));I4_.push_back(0.);I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(1./sqrt(3.));Ihat4_.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ rho
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(113);
  maxweight_.push_back(0.0697802);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(1./sqrt(3.));I4_.push_back(0.);
  I5_.push_back(1./sqrt(12.));
  Ihat3_.push_back(-1./sqrt(3.));Ihat4_.push_back(2./sqrt(3.));
  // lambda_c to sigma+ omega
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(223);
  maxweight_.push_back(0.252355);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.57735033);I4_.push_back(0.);I5_.push_back( 0.288675151);
  Ihat3_.push_back(0.577350303);Ihat4_.push_back( -1.15470061);
  // lambda_c to simga+ phi
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(333);
  maxweight_.push_back(0.0063064);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.816496614);I4_.push_back(-0.816496624);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to p k*0bar
  incoming_.push_back(4122);outgoingB_.push_back(2212);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0996461);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c+ to xi K*+
  incoming_.push_back(4122);outgoingB_.push_back(3322);outgoingM_.push_back(323);
  maxweight_.push_back(0.00413946);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c to sigma+ K*bar0
  incoming_.push_back(4232);outgoingB_.push_back(3222);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0583248);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(6.));Ihat4_.push_back(4./sqrt(6.));
  // xi_c to xi_0 rho+
  incoming_.push_back(4232);outgoingB_.push_back(3322);outgoingM_.push_back(213);
  maxweight_.push_back(2.18754);
  I1_.push_back(1./sqrt(6.));I2_.push_back(1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2./sqrt(6.));Ihat4_.push_back(-4./sqrt(6.));
  // xi_c to lambda k*bar0
  incoming_.push_back(4132);outgoingB_.push_back(3122);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0373243);
  I1_.push_back(1./6.);I2_.push_back(1./6.);
  I3_.push_back(2./3.);I4_.push_back(-1./3.);I5_.push_back(1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // xi_c to sigma k*0
  incoming_.push_back(4132);outgoingB_.push_back(3212);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0135937);
  I1_.push_back(1./sqrt(12.));I2_.push_back(1./sqrt(12.));
  I3_.push_back(0.);I4_.push_back(-2./sqrt(12.));I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(2./sqrt(12.));Ihat4_.push_back(-4./sqrt(12.));
  // xi_c to sigma k*-
  incoming_.push_back(4132);outgoingB_.push_back(3222);outgoingM_.push_back(-323);
  maxweight_.push_back(0.00877507);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi0 rho0
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(113);
  maxweight_.push_back(0.0364247);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(2./sqrt(12.));I4_.push_back(-2./sqrt(12.));I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(12.));Ihat4_.push_back(4./sqrt(12.));
  // xi_c0 to xi0 omega
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(223);
  maxweight_.push_back(0.134395);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back( 0.57735033);I4_.push_back(-0.577350303);I5_.push_back(0.);
  Ihat3_.push_back(0.577350303);Ihat4_.push_back(-1.15470061);
  // xi_c0 to xi0 phi
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(333);
  maxweight_.push_back(0.00676745);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.816496614);I4_.push_back(0.);I5_.push_back(0.408248312);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi- rho+
  incoming_.push_back(4132);outgoingB_.push_back(3312);outgoingM_.push_back(213);
  maxweight_.push_back(0.309733);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c to xi0 K*0
  incoming_.push_back(4332);outgoingB_.push_back(3322);outgoingM_.push_back(-313);
  maxweight_.push_back(0.019967);
  I1_.push_back(1.);I2_.push_back(-1.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2.);Ihat4_.push_back(0.);
  // lambda_c to sigma*0 pi+
  incoming_.push_back(4122);outgoingB_.push_back(3214);outgoingM_.push_back(211);
  maxweight_.push_back(0.0254397);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ pi0
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(111);
  maxweight_.push_back(0.0256531);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ eta
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(221);
  maxweight_.push_back(0.0299659);
  I1_.push_back(0.);I2_.push_back( 0.569036);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ eta'
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(331);
  maxweight_.push_back(4.4186e-07);
  I1_.push_back(0.);I2_.push_back(-0.097631);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta+ kbar0
  incoming_.push_back(4122);outgoingB_.push_back(2214);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0215801);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta++ k-
  incoming_.push_back(4122);outgoingB_.push_back(2224);outgoingM_.push_back(-321);
  maxweight_.push_back(0.0649794);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to xi*0 k+
  incoming_.push_back(4122);outgoingB_.push_back(3324);outgoingM_.push_back(321);
  maxweight_.push_back(0.0198121);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3214);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0118739);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*+ Kbar-
  incoming_.push_back(4132);outgoingB_.push_back(3224);outgoingM_.push_back(-321);
  maxweight_.push_back(0.0240231);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(111);
  maxweight_.push_back(0.0173433);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 eta
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(221);
  maxweight_.push_back(0.000939099);
  I1_.push_back(0.);I2_.push_back(0.097631);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 eta'
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(331);
  maxweight_.push_back(1.31513e-05);
  I1_.push_back(0.);I2_.push_back(-0.569036);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*- pi+
  incoming_.push_back(4132);outgoingB_.push_back(3314);outgoingM_.push_back(211);
  maxweight_.push_back(0.0338039);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to omega- K+
  incoming_.push_back(4132);outgoingB_.push_back(3334);outgoingM_.push_back(321);
  maxweight_.push_back(0.00715116);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to xi*0 kbar0
  incoming_.push_back(4332);outgoingB_.push_back(3324);outgoingM_.push_back(-311);
  maxweight_.push_back(0.00850994);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to omega pi+
  incoming_.push_back(4332);outgoingB_.push_back(3334);outgoingM_.push_back(211);
  maxweight_.push_back(0.012698);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*0 rho+
  incoming_.push_back(4122);outgoingB_.push_back(3214);outgoingM_.push_back(213);
  maxweight_.push_back(0.0478329);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ rho0
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(113);
  maxweight_.push_back(0.0476973);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ omega
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(223);
  maxweight_.push_back(0.0260017);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ phi
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(333);
  maxweight_.push_back(3.55986e-07);
  I1_.push_back(0.);I2_.push_back(-sqrt(2.)/3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta+ kbar0
  incoming_.push_back(4122);outgoingB_.push_back(2214);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0988993);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta++ k-
  incoming_.push_back(4122);outgoingB_.push_back(2224);outgoingM_.push_back(-323);
  maxweight_.push_back(0.303435);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to xi*0 k+
  incoming_.push_back(4122);outgoingB_.push_back(3324);outgoingM_.push_back(323);
  maxweight_.push_back(9.70866e-05);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3214);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0281627);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*- K-
  incoming_.push_back(4132);outgoingB_.push_back(3224);outgoingM_.push_back(-323);
  maxweight_.push_back(0.0575937);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(113);
  maxweight_.push_back(0.0480608);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 omega
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(223);
  maxweight_.push_back(0.0223044);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 phi
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(333);
  maxweight_.push_back(3.187e-09);
  I1_.push_back(0.);I2_.push_back(-sqrt(2.)/3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*- rho+
  incoming_.push_back(4132);outgoingB_.push_back(3314);outgoingM_.push_back(213);
  maxweight_.push_back(0.0958052);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to omega- K*+
  incoming_.push_back(4132);outgoingB_.push_back(3334);outgoingM_.push_back(323);
  maxweight_.push_back(0.000205851);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to xi*0 k*bar0
  incoming_.push_back(4332);outgoingB_.push_back(3324);outgoingM_.push_back(-313);
  maxweight_.push_back(0.050904);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to omega rho+
  incoming_.push_back(4332);outgoingB_.push_back(3334);outgoingM_.push_back(213);
  maxweight_.push_back(0.0814021);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // initial size of the vectors
  initsize_=incoming_.size();
  // intermediates
  generateIntermediates(false);
}

void KornerKramerCharmDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void KornerKramerCharmDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the vectors have the same size
  unsigned int isize=incoming_.size();
  if(isize!=I1_.size()||isize!=I2_.size()||isize!=I3_.size()||isize!=I4_.size()||
     isize!=I5_.size()||isize!=Ihat3_.size()||isize!=Ihat4_.size()||
     isize!=incoming_.size()||isize!=outgoingB_.size()||isize!=outgoingM_.size()||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in"
			  << " KornerKramerCharmDecayer::doinit()" 
			  << Exception::abortnow;
  // compute the various coefficients
  Energy m1,m2,m3,fmes(ZERO); 
  Energy2 P1P2,Qplus,gmes(ZERO);
  double Fnonfact,A3,B3,Ffact[2];
  Energy H2,H3,A2,B2;
  Energy2 A,B;
  int mspin,bspin;
  double chi,
    chiplus( 0.5*(cplus_*(1.+oneNC_)+cminus_*(1.-oneNC_))),
    chiminus(0.5*(cplus_*(1.+oneNC_)-cminus_*(1.-oneNC_))); 
  InvEnergy2 pre(SM().fermiConstant()*sqrt(SM().CKM(1,1)*SM().CKM(0,0)/2.)),mform2[2];
  // testing only
//   pre = SM().fermiConstant()*0.974/sqrt(2.);
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  unsigned int iy;
  for(unsigned int ix=0;ix<isize;++ix) {
    // get the mass of the particles
    extpart[0]=getParticleData(incoming_[ix]);
    extpart[1]=getParticleData(outgoingB_[ix]);
    extpart[2]=getParticleData(outgoingM_[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    addMode(mode,maxweight_[ix],wgt);
    m1=extpart[0]->mass();
    m2=extpart[1]->mass();
    m3=extpart[2]->mass();
    // dot product P1P2
    P1P2 = 0.5*(m1*m1+m2*m2-m3*m3);
    // formfactor for the non-factorizating diagrams and coefficients
    Fnonfact = 1./(1.-m3*m3/mscminus_/mscminus_)*
      (1.-(m1-m2)*(m1-m2)/mscminus_/mscminus_);
    Fnonfact *= Fnonfact*Fnonfact;
    H2 = H2_*Fnonfact;
    H3 = H3_*Fnonfact;
    if(abs(outgoingM_[ix])==211) {
      fmes=fpi_;
      chi=chiplus;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
    }
    else if(abs(outgoingM_[ix])==311) {
      fmes=fk_;
      chi=chiminus;
      mform2[0]=1./(mdcminus_*mdcminus_);
      mform2[1]=1./(mdcplus_ *mdcplus_);
    }
    else if(abs(outgoingM_[ix])==213) {
      gmes = frho_*m3*m3;
      chi=chiplus;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
    }
    else if(abs(outgoingM_[ix])==313) {
      gmes = fKstar_*m3*m3;
      chi=chiminus;
      mform2[0]=1./(mdcminus_*mdcminus_);
      mform2[1]=1./(mdcplus_ *mdcplus_);
    }
    else {
      fmes=ZERO;
      chi=0.;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
      gmes=ZERO;
    }
    // form factor for the factorising diagrams
    for(iy=0;iy<2;++iy) {
      Ffact[iy]  = 1./(1.-m3*m3*mform2[iy])*(1.-(m1-m2)*(m1-m2)*mform2[iy]);
      Ffact[iy] *= Ffact[iy]*Ffact[iy];
    }
    // invariants
    Qplus  = (m1+m2)*(m1+m2)-m3*m3;
    // decide which type of decay
    mspin=getParticleData(outgoingM_[ix])->iSpin();
    bspin=getParticleData(outgoingB_[ix])->iSpin();
    if(bspin==2) {
      if(mspin==1) {
	// non-factorizing piece of the diagrams
	A = -0.75*H2/m1/m2*cminus_*( (m1*P1P2-m1*m2*(m2+m3))*I3_[ix]
				     -(m2*P1P2-m1*m2*(m1+m3))*Ihat3_[ix]);
	B = 0.25/m1/m2*cminus_*( 0.5*H2*Qplus*( m1*(   I3_[ix]+2.*    I4_[ix])
						+m2*(Ihat3_[ix]+2.*Ihat4_[ix]))
				 +H3*m1*m2*(m1+m2+m3)*12.*I5_[ix]);
	// the factorizing piece
	A +=0.25/m1/m2*chi*fmes*Qplus*(m1-m2)*(2.*I1_[ix]+I2_[ix])*Ffact[0];
	B +=0.25/m1/m2*chi*fmes*Qplus*(m1+m2)/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[1];
	// add to vectors
	A1_.push_back(A*pre);B1_.push_back(B*pre);
	A2_.push_back(ZERO);B2_.push_back(ZERO);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else if(mspin==3) {
	// non-factorizing terms
	A  = -0.25*H2/m1/m2*cminus_*
	  ( (m1*P1P2-m1*m2*(m2+m3))*(   I3_[ix]+2.*   I4_[ix])
	    +(m2*P1P2-m1*m2*(m1+m3))*(Ihat3_[ix]+2.*Ihat4_[ix]));
	A2 = -0.25*H2/m1/m2*cminus_*(m1+m2+m3)*( m1*(   I3_[ix]+2.*   I4_[ix])
						 -m2*(Ihat3_[ix]+2.*Ihat4_[ix]));
	B  = +0.25/m1/m2*cminus_*(0.5*H2*Qplus*( m1*(   I3_[ix]+2.*   I4_[ix])
						 +m2*(Ihat3_[ix]+2.*Ihat4_[ix]))
				  +H3*m1*m2*(m1+m2+m3)*12.*I5_[ix]);
	B2 = -0.25/m1/m2*cminus_*
	  ( H2*( m1*(m1+m2)*(   I3_[ix]+2.*   I4_[ix])-3.*m1*m3*   I3_[ix]
		 +m2*(m1+m2)*(Ihat3_[ix]+2.*Ihat4_[ix])-3.*m2*m3*Ihat3_[ix])
	    +H3*m1*m2*24.*I5_[ix]);
	// the factorizing piece
	A  += -0.25/m1/m2*chi*gmes*Qplus/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[1];
	B  +=  0.25/m1/m2*chi*gmes*Qplus/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[0];
	B2 +=  0.25/m1/m2*chi*gmes*(m1+m2)*4./3.*(I1_[ix]-I2_[ix])*Ffact[0];
	// add to vectors
	A1_.push_back(A*pre);B1_.push_back(B*pre);
	A2_.push_back(A2*pre);B2_.push_back(B2*pre);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else
	throw InitException() << "Invalid outgoing meson spin in"
			      << " KornerKramerCharmDecayer::doinit()" 
			      << Exception::abortnow;
    }
    else if(bspin==4) {
      if(mspin==1) {
	// first the non-factorizing piece
	B2  = -1.5*cminus_*H2*I2_[ix];
	// then the factorizing piece
	B2  += chi*fmes*(1.+m2/m1)*I1_[ix]*Ffact[1];
	// add to vectors
	// make the coupling dimensionless
	A1_.push_back(0.);B1_.push_back(0.);
	A2_.push_back(ZERO);B2_.push_back(B2*pre);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else if(mspin==3) {
	InvEnergy norm(0.75/m1/m2*cminus_*H2*I2_[ix]);
	// first the non-factorizing piece
	A  = -norm*m1*(P1P2-m2*(m2+m3))*2.;
	A2 = ZERO;
	A3 =  norm*m1*2.;
	B  = -norm*m1*Qplus;
	B2 = -norm*m1*m2*2.;
	B3 =  norm*m1*2.;
	// then the factorizing piece
	double norm2 = 0.5/m1/m2*chi*gmes*I1_[ix];
	A  += norm2*Qplus*Ffact[1];
      //A2 += ZERO;
	A3 +=-norm2*2.*Ffact[1];
	B  += norm2*Qplus*Ffact[0];
	B2 += norm2*m2*2.*Ffact[0];
	B3 +=-norm2*2.*Ffact[0];
	// add to vectors
	A1_.push_back( A*pre);B1_.push_back( B*pre);
	A2_.push_back(2.*A2*pre);B2_.push_back(2.*B2*pre);
	A3_.push_back(A3*pre);B3_.push_back(B3*pre);
      }
      else
	throw InitException() << "Invalid outgoing meson spin in"
			      << " KornerKramerCharmDecayer::doinit()" 
			      << Exception::abortnow;
    }
    else
      throw InitException() << "Invalid outgoing baryon spin in"
			    << " KornerKramerCharmDecayer::doinit()" 
			    << Exception::abortnow;
  }
}

int KornerKramerCharmDecayer::modeNumber(bool & cc,tcPDPtr parent,
					 const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoingB_[ix]&&id2==outgoingM_[ix])||
	 (id2==outgoingB_[ix]&&id1==outgoingM_[ix])) imode=ix;
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoingB_[ix]&&id2==-outgoingM_[ix])||
	 (id2==-outgoingB_[ix]&&id1==-outgoingM_[ix])) imode=ix;
      if(((id1==-outgoingB_[ix]&&id2==outgoingM_[ix])||
	  (id2==-outgoingB_[ix]&&id1==outgoingM_[ix]))&&
	 (outgoingM_[ix]==111||outgoingM_[ix]==221||outgoingM_[ix]==331||
	  outgoingM_[ix]==113||outgoingM_[ix]==223||outgoingM_[ix]==333))
	imode=ix;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  // charge conjugation if needed
  cc = id0<0;
  // return mode number
  return imode;
}

void KornerKramerCharmDecayer::persistentOutput(PersistentOStream & os) const {
  os << oneNC_ << ounit(fpi_,GeV) << ounit(fk_,GeV) << frho_ 
     << fKstar_ << ounit(mdcplus_,GeV) << ounit(mdcminus_,GeV) << ounit(mscplus_,GeV) 
     << ounit(mscminus_,GeV) << cplus_ << cminus_ << ounit(H2_,GeV) << ounit(H3_,GeV) 
     << I1_ << I2_ << I3_ << I4_ << I5_ << Ihat3_ << Ihat4_ << incoming_ << outgoingB_ 
     << outgoingM_ << maxweight_ << A1_ << ounit(A2_,1./GeV) << ounit(A3_,1./GeV2) 
     << B1_ << ounit(B2_,1./GeV) << ounit(B3_,1./GeV2) << initsize_;
}

void KornerKramerCharmDecayer::persistentInput(PersistentIStream & is, int) {
  is >> oneNC_ >> iunit(fpi_,GeV) >> iunit(fk_,GeV) >> frho_ 
     >> fKstar_ >> iunit(mdcplus_,GeV) >> iunit(mdcminus_,GeV) >> iunit(mscplus_,GeV) 
     >> iunit(mscminus_,GeV) >> cplus_ >> cminus_ >> iunit(H2_,GeV) >> iunit(H3_,GeV) 
     >> I1_ >> I2_ >> I3_ >> I4_ >> I5_ >> Ihat3_ >> Ihat4_ >> incoming_ >> outgoingB_ 
     >> outgoingM_ >> maxweight_ >> A1_ >> iunit(A2_,1./GeV) >> iunit(A3_,1./GeV2) 
     >> B1_ >> iunit(B2_,1./GeV) >> iunit(B3_,1./GeV2) >> initsize_;
}

ClassDescription<KornerKramerCharmDecayer> 
KornerKramerCharmDecayer::initKornerKramerCharmDecayer;
// Definition of the static class description member.

void KornerKramerCharmDecayer::Init() {

  static ClassDocumentation<KornerKramerCharmDecayer> documentation
    ("The KornerKramerCharmDecayer class implements the"
     " non-leptonic weak decay of charm baryons using the results of Z.Phys.C55,659.",
     "The non-leptonic charm decays were simulated using the KornerKramerCharmDecayer"
     "class which implements the model of \\cite{Korner:1992wi}.",
     "\\bibitem{Korner:1992wi}\n"
     "J.~G.~Korner and M.~Kramer,\n"
     "Z.\\ Phys.\\  C {\\bf 55} (1992) 659.\n"
     "%%CITATION = ZEPYA,C55,659;%%\n");

  static Parameter<KornerKramerCharmDecayer,double> interfaceOneOverNc
    ("OneOverNc",
     "One divided by the number of colours, the default is to take N_c to infinity.",
     &KornerKramerCharmDecayer::oneNC_, 0.0, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFpi
    ("Fpi",
     "The decay constant for the pi meson.",
     &KornerKramerCharmDecayer::fpi_, MeV, 131.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFK
    ("FK",
     "The decay constant for the K meson",
     &KornerKramerCharmDecayer::fk_, MeV, 160.6*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceFrho
    ("Frho",
     "The decay constant for the rho meson",
     &KornerKramerCharmDecayer::frho_, 0.272, 0.0, 1.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfacefKstar
    ("fKstar",
     "The decay constant for the K star meson",
     &KornerKramerCharmDecayer::fKstar_, 0.238, 0.0, 1.0,
     false, false, false);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcplus
    ("Mdcplus",
     "The mass of the 1+ dc meson for the form-factors",
     &KornerKramerCharmDecayer::mdcplus_, GeV, 2.42*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscplus
    ("Mscplus",
     "The mass of the 1+ sc meson for the form-factors",
     &KornerKramerCharmDecayer::mscplus_, GeV, 2.54*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcminus
    ("Mdcminus",
     "The mass of the 1- dc meson for the form-factors",
     &KornerKramerCharmDecayer::mdcminus_, GeV, 2.01*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscminus
    ("Mscminus",
     "The mass of the 1- sc meson for the form-factors",
     &KornerKramerCharmDecayer::mscminus_, GeV, 2.11*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCplus
    ("Cplus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::cplus_, 0.73, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCminus
    ("Cminus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::cminus_, 1.90, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH2
    ("H2",
     "The H2 parameter",
     &KornerKramerCharmDecayer::H2_, GeV, 0.119*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH3
    ("H3",
     "The H3 parameter",
     &KornerKramerCharmDecayer::H3_, GeV,-0.011*GeV,-10.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI1
    ("I1",
     "The I_1 invariant for the decay modes",
     &KornerKramerCharmDecayer::I1_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI2
    ("I2",
     "The I_2 invariant for the decay modes",
     &KornerKramerCharmDecayer::I2_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI3
    ("I3",
     "The I_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::I3_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI4
    ("I4",
     "The I_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::I4_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI5
    ("I5",
     "The I_5 invariant for the decay modes",
     &KornerKramerCharmDecayer::I5_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat3
    ("Ihat3",
     "The Ihat_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::Ihat3_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat4
    ("Ihat4",
     "The Ihat_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::Ihat4_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &KornerKramerCharmDecayer::incoming_, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingB
    ("OutgoingB",
     "The PDG code of the outgoing baryon",
     &KornerKramerCharmDecayer::outgoingB_, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingM
    ("OutgoingM",
     "The PDG code of the outgoing meson",
     &KornerKramerCharmDecayer::outgoingM_, -1, 0, -1000000, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &KornerKramerCharmDecayer::maxweight_,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void KornerKramerCharmDecayer::
 halfHalfScalarCoupling(int imode,Energy,Energy,Energy,
			Complex& A,Complex&B) const {
  useMe();
  A = A1_[imode];
  B = B1_[imode];
}

// couplings for spin-1/2 to spin-1/2 spin-1
void KornerKramerCharmDecayer::halfHalfVectorCoupling(int imode,
						      Energy m0,Energy m1,Energy,
						      Complex& A1,Complex& A2,
						      Complex& B1,Complex& B2) const {
  useMe();
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A1 = B1_[imode];
  B1 = A1_[imode];
  A2 = B2_[imode]*(m0+m1);
  B2 = A2_[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-0
void KornerKramerCharmDecayer::halfThreeHalfScalarCoupling(int imode, Energy m0,
							   Energy m1,Energy,Complex& A,
							   Complex& B) const {
  useMe();
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A = B2_[imode]*(m0+m1);
  B = A2_[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-1
void KornerKramerCharmDecayer::
halfThreeHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex& A1, Complex& A2, Complex& A3,
			    Complex& B1, Complex& B2, Complex& B3) const {
  useMe();
  A1 = A1_[imode];
  A2 = A2_[imode]*(m0+m1);
  A3 = A3_[imode]*(m0+m1)*(m0+m1);
  B1 = B1_[imode];
  B2 = B2_[imode]*(m0+m1);
  B3 = B3_[imode]*(m0+m1)*(m0+m1);
}
 
void KornerKramerCharmDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":OneOverNc " <<  oneNC_ << "\n";
  output << "newdef " << name() << ":Fpi " << fpi_/MeV << "\n";
  output << "newdef " << name() << ":FK " << fk_/MeV  << "\n";
  output << "newdef " << name() << ":Frho " << frho_ << "\n";
  output << "newdef " << name() << ":fKstar " << fKstar_ << "\n";
  output << "newdef " << name() << ":Mdcplus " << mdcplus_/GeV << "\n";
  output << "newdef " << name() << ":Mscplus " << mscplus_/GeV << "\n";
  output << "newdef " << name() << ":Mdcminus " << mdcminus_/GeV << "\n";
  output << "newdef " << name() << ":Mscminus " << mscminus_/GeV << "\n";
  output << "newdef " << name() << ":Cplus " << cplus_ << "\n";
  output << "newdef " << name() << ":Cminus " << cminus_ << "\n";
  output << "newdef " << name() << ":H2 " << H2_/GeV << "\n";
  output << "newdef " << name() << ":H3 " << H3_/GeV << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(ix<initsize_) {
      output << "newdef " << name() << ":I1 " 
	     << ix << " " << I1_[ix] << "\n";
      output << "newdef " << name() << ":I2 " 
	     << ix << " " << I2_[ix] << "\n";
      output << "newdef " << name() << ":I3 " 
	     << ix << " " << I3_[ix] << "\n";
      output << "newdef " << name() << ":I4 " 
	     << ix << " " << I4_[ix] << "\n";
      output << "newdef " << name() << ":I5 " 
	     << ix << " " << I5_[ix] << "\n";
      output << "newdef " << name() << ":Ihat3 " 
	     << ix << " " << Ihat3_[ix] << "\n";
      output << "newdef " << name() << ":Ihat4 " 
	     << ix << " " << Ihat4_[ix] << "\n";
      output << "newdef " << name() << ":Incoming " 
	     << ix << " " << incoming_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingB " 
	     << ix << " " << outgoingB_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingM " 
	     << ix << " " << outgoingM_[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " 
	     << ix << " " << maxweight_[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":I1 " 
	     << ix << " " << I1_[ix] << "\n";
      output << "insert " << name() << ":I2 " 
	     << ix << " " << I2_[ix] << "\n";
      output << "insert " << name() << ":I3 " 
	     << ix << " " << I3_[ix] << "\n";
      output << "insert " << name() << ":I4 " 
	     << ix << " " << I4_[ix] << "\n";
      output << "insert " << name() << ":I5 " 
	     << ix << " " << I5_[ix] << "\n";
      output << "insert " << name() << ":Ihat3 " 
	     << ix << " " << Ihat3_[ix] << "\n";
      output << "insert " << name() << ":Ihat4 " 
	     << ix << " " << Ihat4_[ix] << "\n";
      output << "insert " << name() << ":Incoming " 
	     << ix << " " << incoming_[ix] << "\n";
      output << "insert " << name() << ":OutgoingB " 
	     << ix << " " << outgoingB_[ix] << "\n";
      output << "insert " << name() << ":OutgoingM " 
		 << ix << " " << outgoingM_[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " 
		 << ix << " " << maxweight_[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
