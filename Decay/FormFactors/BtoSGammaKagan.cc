// -*- C++ -*-
//
// BtoSGammaKagan.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaKagan class.
//

#include "BtoSGammaKagan.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/GaussianIntegrator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using Herwig::Math::Li2;

DescribeClass<BtoSGammaKagan,BtoSGammaHadronicMass>
describeHerwigBtoSGammaKagan("Herwig::BtoSGammaKagan",
			     "HwFormFactors.so");

HERWIG_INTERPOLATOR_CLASSDESC(BtoSGammaKagan1,double,double)

HERWIG_INTERPOLATOR_CLASSDESC(BtoSGammaKagan2,InvEnergy,Energy)


BtoSGammaKagan::BtoSGammaKagan() 
  : _initialize(false),_mt(175.*GeV),_mb(4.8*GeV),
    _mc(1.392*GeV),_ms(ZERO),_msovermb(1./50.),_zratio(0.),
    _lambda2(0.12*GeV2),_mw(80.425*GeV),_mz(91.1876*GeV),
    _MB(5279.4*MeV),_c20(0.),_c70(0.),_c80(0.),
    _beta0(23./3.),_beta1(116./3.),_alpha(1./137.036),
    _alphaSZ(0.118),_mub(4.8*GeV),_alphaSM(0.),
    _ckm(0.976),_delta(0.),_spectmax(0.00025/GeV),_maxtry(100),
    _fermilambda(ZERO),_fermia(0.),_ferminorm(1./GeV),
    _fermilambda1(-0.3*GeV2),_ycut(0.9999999999),
    _y(0.),_deltacut(0.9),_nsfunct(100),_nspect(100),_iopt(9999) {
  Energy mHin[100]={0*GeV,0.0505907*GeV,0.101181*GeV,0.151772*GeV,0.202363*GeV,
		    0.252953*GeV,0.303544*GeV,0.354135*GeV,0.404726*GeV,0.455316*GeV,
		    0.505907*GeV,0.556498*GeV,0.607088*GeV,0.657679*GeV,0.70827*GeV,
		    0.75886*GeV,0.809451*GeV,0.860042*GeV,0.910632*GeV,0.961223*GeV,
		    1.01181*GeV,1.0624*GeV,1.113*GeV,1.16359*GeV,1.21418*GeV,
		    1.26477*GeV,1.31536*GeV,1.36595*GeV,1.41654*GeV,1.46713*GeV,
		    1.51772*GeV,1.56831*GeV,1.6189*GeV,1.66949*GeV,1.72008*GeV,
		    1.77067*GeV,1.82126*GeV,1.87186*GeV,1.92245*GeV,1.97304*GeV,
		    2.02363*GeV,2.07422*GeV,2.12481*GeV,2.1754*GeV,2.22599*GeV,
		    2.27658*GeV,2.32717*GeV,2.37776*GeV,2.42835*GeV,2.47894*GeV,
		    2.52953*GeV,2.58013*GeV,2.63072*GeV,2.68131*GeV,2.7319*GeV,
		    2.78249*GeV,2.83308*GeV,2.88367*GeV,2.93426*GeV,2.98485*GeV,
		    3.03544*GeV,3.08603*GeV,3.13662*GeV,3.18721*GeV,3.2378*GeV,
		    3.2884*GeV,3.33899*GeV,3.38958*GeV,3.44017*GeV,3.49076*GeV,
		    3.54135*GeV,3.59194*GeV,3.64253*GeV,3.69312*GeV,3.74371*GeV,
		    3.7943*GeV,3.84489*GeV,3.89548*GeV,3.94607*GeV,3.99666*GeV,
		    4.04726*GeV,4.09785*GeV,4.14844*GeV,4.19903*GeV,4.24962*GeV,
		    4.30021*GeV,4.3508*GeV,4.40139*GeV,4.45198*GeV,4.50257*GeV,
		    4.55316*GeV,4.60375*GeV,4.65434*GeV,4.70493*GeV,4.75553*GeV,
		    4.80612*GeV,4.85671*GeV,4.9073*GeV,4.95789*GeV,5.00848*GeV};
  InvEnergy spin[100]={0./GeV,3.40885e-10/GeV,1.06258e-08/GeV,7.30539e-08/GeV,
		       2.75462e-07/GeV,7.49796e-07/GeV,1.66662e-06/GeV,3.22116e-06/GeV,
		       5.62241e-06/GeV,9.06604e-06/GeV,1.37419e-05/GeV,1.98035e-05/GeV,
		       2.73352e-05/GeV,3.6404e-05/GeV,4.69896e-05/GeV,5.90282e-05/GeV,
		       7.23334e-05/GeV,8.67468e-05/GeV,0.000101999/GeV,0.000117792/GeV,
		       0.000133824/GeV,0.0001497/GeV,0.000165089/GeV,0.000179618/GeV,
		       0.000193034/GeV,0.000205017/GeV,0.000215324/GeV,0.000223761/GeV,
		       0.000230224/GeV,0.00023456/GeV,0.000236774/GeV,0.000236858/GeV,
		       0.000234965/GeV,0.000231179/GeV,0.00022566/GeV,0.000218597/GeV,
		       0.000210199/GeV,0.000200691/GeV,0.000190323/GeV,0.000179277/GeV,
		       0.000167797/GeV,0.000156088/GeV,0.000144322/GeV,0.000132705/GeV,
		       0.000121364/GeV,0.00011042/GeV,9.99745e-05/GeV,9.01017e-05/GeV,
		       8.08564e-05/GeV,7.22729e-05/GeV,6.43679e-05/GeV,5.71471e-05/GeV,
		       5.05892e-05/GeV,4.46771e-05/GeV,3.93802e-05/GeV,3.46595e-05/GeV,
		       3.04797e-05/GeV,2.67948e-05/GeV,2.3561e-05/GeV,2.07344e-05/GeV,
		       1.82726e-05/GeV,1.61351e-05/GeV,1.42839e-05/GeV,1.26838e-05/GeV,
		       1.13031e-05/GeV,1.01119e-05/GeV,9.08476e-06/GeV,8.19854e-06/GeV,
		       7.43307e-06/GeV,6.77062e-06/GeV,6.19614e-06/GeV,5.69634e-06/GeV,
		       5.26005e-06/GeV,4.87775e-06/GeV,4.54138e-06/GeV,4.24416e-06/GeV,
		       3.9805e-06/GeV,3.74572e-06/GeV,3.53602e-06/GeV,3.34835e-06/GeV,
		       3.18023e-06/GeV,3.02972e-06/GeV,2.89566e-06/GeV,2.77756e-06/GeV,
		       2.67492e-06/GeV,2.58796e-06/GeV,2.52423e-06/GeV,2.51738e-06/GeV,
		       2.58615e-06/GeV,2.72927e-06/GeV,2.94112e-06/GeV,3.22002e-06/GeV,
		       3.57015e-06/GeV,4.00189e-06/GeV,4.53288e-06/GeV,5.19051e-06/GeV,
		       6.0167e-06/GeV,7.0773e-06/GeV,8.4804e-06/GeV,1.04156e-05/GeV};
  _mHinter=vector<Energy>(mHin,mHin+100);
  _spectrum=vector<InvEnergy>(spin,spin+100);
}

void BtoSGammaKagan::doinitrun() {
  BtoSGammaHadronicMass::doinitrun();
  _pmHinter = new_ptr(Interpolator<InvEnergy,Energy>(_spectrum,_mHinter,3));
}

void BtoSGammaKagan::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mt,GeV) << ounit(_mb,GeV) << ounit(_mc,GeV) << ounit(_ms,GeV) 
     << _msovermb << _zratio << ounit(_lambda2,GeV2) << ounit(_mw,GeV) << ounit(_mz,GeV) 
     << ounit(_MB,GeV) << _c20 << _c70 << _c80 << _beta0 << _beta1 << _alpha 
     << _alphaSZ << ounit(_mub,GeV) << _alphaSM << _ckm << _delta 
     << ounit(_mHinter,GeV) << ounit(_spectrum,1./GeV) 
     << ounit(_spectmax,1./GeV) << ounit(_fermilambda,GeV)
     << _fermia << ounit(_ferminorm,1./GeV) << ounit(_fermilambda1,GeV2)
     <<_ycut << _deltacut << _nsfunct 
     << _nspect << _maxtry << _initialize
     << _s22inter << _s27inter << _pmHinter;
}

void BtoSGammaKagan::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mt,GeV) >> iunit(_mb,GeV) >> iunit(_mc,GeV) >> iunit(_ms,GeV) 
     >> _msovermb >> _zratio >> iunit(_lambda2,GeV2) >> iunit(_mw,GeV) >> iunit(_mz,GeV) 
     >> iunit(_MB,GeV) >> _c20 >> _c70 >> _c80 >> _beta0 >> _beta1 >> _alpha 
     >> _alphaSZ >> iunit(_mub,GeV) >> _alphaSM >> _ckm >> _delta 
     >> iunit(_mHinter,GeV) >> iunit(_spectrum,1./GeV)
     >> iunit(_spectmax,1./GeV) >> iunit(_fermilambda,GeV) 
     >> _fermia >> iunit(_ferminorm,1./GeV) >> iunit(_fermilambda1,GeV2)
     >>_ycut >> _deltacut >> _nsfunct 
     >> _nspect >> _maxtry >> _initialize
     >> _s22inter >> _s27inter >> _pmHinter;
}

void BtoSGammaKagan::Init() {

  static ClassDocumentation<BtoSGammaKagan> documentation
    ("The BtoSGammaKagan class implements the calculation of hep-ph/9805303 for the"
     " hadronic mass spectrum in b to s gamma decays.",
     "The decay $b\\to s\\gamma$ was simulated using the hadronic mass spectrum from"
     "\\cite{Kagan:1998ym}.\n",
     "\\bibitem{Kagan:1998ym} A.~L.~Kagan and M.~Neubert,\n"
     "Eur.\\ Phys.\\ J.\\  C {\\bf 7} (1999) 5 [arXiv:hep-ph/9805303].\n"
     "%%CITATION = EPHJA,C7,5;%%\n");

  static Switch<BtoSGammaKagan,bool> interfaceInitialize
    ("Initialize",
     "Initialize the interpolation tables for the hadronic mass",
     &BtoSGammaKagan::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialize
    (interfaceInitialize,
     "Yes",
     "Perform initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "No",
     "No initialization is performed",
     false);

  static Parameter<BtoSGammaKagan,Energy> interfaceTopMass
    ("TopMass",
     "The mass of the top quark",
     &BtoSGammaKagan::_mt, GeV, 175.*GeV, 100.0*GeV, 200.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark",
     &BtoSGammaKagan::_mb, GeV, 4.8*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark",
     &BtoSGammaKagan::_mc, GeV, 1.392*GeV, 1.0*GeV, 2.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceStrangeMassRatio
    ("StrangeMassRatio",
     "The ratio of the strange quark mass to that of the bottom",
     &BtoSGammaKagan::_msovermb, 1./50., 0.001, 0.1,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceWMass
    ("WMass",
     "The mass of the W boson",
     &BtoSGammaKagan::_mw, GeV, 80.425*GeV, 75.*GeV, 85.*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceZMass
    ("ZMass",
     "The mass of the Z boson",
     &BtoSGammaKagan::_mz, GeV, 91.1876*GeV, 85.*GeV, 95.*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy2> interfaceLambda2
    ("Lambda2",
     "Hadronic parameter from hep-ph/9805303",
     &BtoSGammaKagan::_lambda2, GeV2, 0.12*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceBMesonMass
    ("BMesonMass",
     "The mass of the decaying B meson",
     &BtoSGammaKagan::_MB, GeV, 5.2794*GeV, 5.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceMu
    ("Mu",
     "The renormalisation scale",
     &BtoSGammaKagan::_mub, GeV, 4.8*GeV, 3.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceDelta
    ("Delta",
     "The cut off on the photon energy",
     &BtoSGammaKagan::_deltacut, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy2> interfaceLambda1
    ("Lambda1",
     "Hadronic scale for the fermi motion function",
     &BtoSGammaKagan::_fermilambda1, GeV2, -0.3*GeV2, -10.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfacealpha
    ("alpha",
     "The fine structure constant",
     &BtoSGammaKagan::_alpha, 1./137.036, 0.005, 0.01,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfacealphaSZ
    ("alphaSZ",
     "The strong coupling at the Z mass",
     &BtoSGammaKagan::_alphaSZ, 0.118, 0.1, 0.3,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceCKM
    ("CKM",
     "The CKM prefactor for the decay",
     &BtoSGammaKagan::_ckm, 0.976, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<BtoSGammaKagan,Energy> interfacemHValues
    ("mHValues",
     "The mH values for the interpolation of the spectrum",
     &BtoSGammaKagan::_mHinter, GeV, -1, 1.*GeV, ZERO, 10.*GeV,
     false, false, Interface::limited);

  static ParVector<BtoSGammaKagan,InvEnergy> interfaceSpectrum
    ("Spectrum",
     "Values of the spectrum for interpolation",
     &BtoSGammaKagan::_spectrum, 1./GeV, -1, ZERO, ZERO, 1./GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,InvEnergy> interfaceFermiNormalisation
    ("FermiNormalisation",
     "The normalisation for the fermi motion function",
     &BtoSGammaKagan::_ferminorm, 1./GeV, 1./GeV, 0*1./GeV, 0*1./GeV,
     false, false, Interface::nolimits);

  static Parameter<BtoSGammaKagan,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight the spectrum",
     &BtoSGammaKagan::_maxtry, 100, 10, 10000,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,InvEnergy> interfaceSpectrumMaximum
    ("SpectrumMaximum",
     "The maximum value of the spectrum for unweighting",
     &BtoSGammaKagan::_spectmax, 1./GeV, 1./GeV, ZERO, 10000.0/GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceycut
    ("ycut",
     "Limit of the value of y to avoid singualarities in integrals",
     &BtoSGammaKagan::_ycut, 0.9999999999, 0.999, 1.,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,unsigned int> interfaceNSpectrum
    ("NSpectrum",
     "Number of spectrum points to calculate for interpolation",
     &BtoSGammaKagan::_nspect, 100, 10, 1000,
     false, false, Interface::limited);
}

void BtoSGammaKagan::calculateWilsonCoefficients() {
  // strong coupling at various scales
  double alphaSMW(alphaS(_mw)),alphaSMT(alphaS(_mt));
  _alphaSM=alphaS(_mub);
  // coupling running top mass and ratio to W mass
  Energy mtbar=_mt*pow(alphaSMW/alphaSMT,12./23.)*
    (1.+12./23.*(253./72.-58./56.)/pi*(alphaSMW-alphaSMT)-4./3.*alphaSMT/pi);
  double xt(mtbar/_mw);xt*=xt;
  // leading order coefficients
  double eta(alphaSMW/_alphaSM);
  _c20=0.5*(pow(eta,-12./23.)+pow(eta,6./23.));
  double c7w=1./pow(xt-1.,3)*(0.25*xt*xt*(3.*xt-2.)/(xt-1)*log(xt)
			      +xt/24.*(7.-5.*xt-8.*xt*xt));
  double c8w=0.25*xt/pow(xt-1.,3)*(-3.*xt/(xt-1.)*log(xt)+0.5*(2.+5.*xt-xt*xt));
  // corrections
  double c7c=626126./272277.*pow(eta,14./23.)-56281./51730.*pow(eta,16./23.)
    -3./7.*pow(eta,6./23.)-1./14.*pow(eta,-12./23.)-0.6494*pow(eta,0.4086)
    -0.0380*pow(eta,-.4230)-0.0186*pow(eta,-0.8994)-0.0057*pow(eta,0.1456);
  double c8c=313063./363036.*pow(eta,14./23.)-0.9135*pow(eta,0.4086)
    +0.0873*pow(eta,-0.4230)-0.0571*pow(eta,-0.8994)+0.0209*pow(eta,0.1456);
  // total
  _c70=pow(eta,16./23.)*c7w+8./3.*(pow(eta,14./23.)-pow(eta,16./23.))*c8w+c7c;
  _c80=pow(eta,14./23.)*c8w+c8c;
  // corrections
  double sp(real(Li2(Complex(1.-1./xt))));
  double xt2(xt*xt),xt3(xt2*xt),xt4(xt3*xt),xt5(xt4*xt);
  double c71w=1./pow(xt-1.,4)*( sp/9.*xt*(-8.+80.*xt-122.*xt2-16.*xt3)
				+pow(log(xt),2)/3./(xt-1.)*xt2*(-28.+46.*xt+6.*xt2)
				+log(xt)/81./(xt-1.)*(208.-1364.*xt+3244.*xt2
						      -2262.*xt3-588.*xt4-102.*xt5)
				+1./486.*(-436.+2509.*xt-10740.*xt2+12205.*xt3+1646.*xt4)
				);
  double c81w=1./pow(xt-1.,4)*(sp/6.*(xt+41.*xt2+40.*xt3-4.*xt4)
			       +0.5*pow(log(xt),2)/(xt-1.)*(-17.*xt3-31.*xt2)
			       +log(xt)/216./(xt-1.)*(280.-1994.*xt+2857.*xt2
						      +4893.*xt3+1086.*xt4-210.*xt5)
			       +1./1296.*(-508.+610.*xt-28209.*xt2-14102.*xt3+737.*xt4));
  // coefficients for the series in eta
  double ai[8]={14./23.,16./23.,6./23.,-12./23.,0.4086,-0.4230,-0.8994,0.1456};
  double ei[8]={4661194./816831.,-8516./2217.,0.,0.,-1.9043,-0.1008,0.1216,0.0183};
  double fi[8]={-17.3023,8.5027,4.5508,0.7519,2.0040,0.7476,-0.5385,0.0914};
  double gi[8]={14.8088,-10.8090,-0.8740,0.4218,-2.9347,0.3971,0.1600,0.0225};
  double Ex=-2./3.*log(xt)+log(xt)/6./pow(1.-xt,4)*xt2*(15.-16.*xt+4.*xt2)
    +xt/12./pow(1.-xt,3)*(18.-11.*xt-xt2);
  double c71c=(297664./14283.*pow(eta,16./23.)-7164416./357075.*pow(eta,14./23.)
	       +256868./14283.*pow(eta,37./23.)-6698884./357075.*pow(eta,39./23.))*c8w
    +37208./4761.*(pow(eta,39./23.)-pow(eta,16./23.))*c7w;
  for(unsigned int ix=0;ix<8;++ix) {
    c71c+=(ei[ix]*eta*Ex+fi[ix]+gi[ix]*eta)*pow(eta,ai[ix]);
  }
  // total correction
  double c71(pow(eta,39./23.)*c71w+8./3.*(pow(eta,37./23.)-pow(eta,39./23.))*c81w+c71c);
  // electromagnetic corrections
  double c7em=
    (32./75.*pow(eta,-9./23.)-40./69.*pow(eta,-7./23.)+88./575.*pow(eta,16./23.))*c7w
    +(-32./575.*pow(eta,-9./23.)+32./1449.*pow(eta,-7./23.)+640./1449.*pow(eta,14./23.)
      -704./1725.*pow(eta,16./23.))*c8w
    -190./8073.*pow(eta,-35./23.)-359./3105.*pow(eta,-17./23.)
    +4276./121095.*pow(eta,-12./23.)+350531./1009125.*pow(eta,-9./23.)
    +2./4347.*pow(eta,-7./23.)-5956./15525.*pow(eta,6./23.)
    +38380./169533.*pow(eta,14./23.)-748./8625.*pow(eta,16./23.);
  // nlo correction to semi-leptonic decay
  double kSL(12./23.*(1./eta -1.));
  // corrections 
  double r7(-10./3.-8.*pi*pi/9.),gam27(416./81.),gam77(32./3.),gam87(-32./9.);
  double fz(semiLeptonicf());
  double kappa(3.382-4.14*(sqrt(_zratio)-0.29));
  double realr2(-4.092+12.78*(sqrt(_zratio)-0.29));
  double realr8(44./9.-8.*pi*pi/27.);
  // correction to the various terms (get multiplied by Delta(y))
  double delta77=1.+0.5*_alphaSM/pi*(r7+gam77*log(_mb/_mub)-16./3.)
    +(pow(1.-_zratio,4)/fz-1.)*6.*_lambda2/_mb/_mb
    +0.5*_alphaSM/pi*kappa;
  double delta27=0.5*_alphaSM/pi*(realr2+gam27*log(_mb/_mub))-_lambda2/9./_mc/_mc;
  double delta78=0.5*_alphaSM/pi*(realr8+gam87*log(_mb/_mub));
  // combination for the genuine NLO terms
  _delta = delta77*_c70*_c70+delta27*_c20*_c70+delta78*_c70*_c80
    +0.5*_alphaSM/pi*c71*_c70
    +_alpha/_alphaSM*(2.*c7em*_c70-kSL*_c70*_c70);
}

void BtoSGammaKagan::doinit() {
  BtoSGammaHadronicMass::doinit();
  // quark masses
  _ms=_msovermb*_mb;
  _zratio=sqr(_mc/_mb);
  // parameters for the fermi motion function
  _fermilambda=_MB-_mb;
  _fermia=-3.*sqr(_fermilambda)/_fermilambda1-1.;
  if(_initialize) {
    // calculate the wilson coefficients etc.
    calculateWilsonCoefficients();
    // calculate the interpolation tables for the s functions
    vector<double> sfunct,yvalues;
    // s_22 function
    sfunct.push_back(0.),yvalues.push_back(0.);
    double step(1./_nsfunct);
    _y=-0.5*step;
    // perform the integrals
    GaussianIntegrator integrator;
    _iopt=0;
    for(unsigned int ix=0;ix<_nsfunct;++ix) {
      _y+=step;
      yvalues.push_back(_y);
      sfunct.push_back(integrator.value(*this,0.,_y));
    }
    _s22inter = new_ptr(Interpolator<double,double>(sfunct,yvalues,3));
    // s_27 function;
    _iopt=1;
    sfunct.clear();yvalues.clear();
    sfunct.push_back(0.),yvalues.push_back(0.);
    _y=-0.5*step;
    // perform integrals
    for(unsigned int ix=0;ix<_nsfunct;++ix) {
      _y+=step;
      yvalues.push_back(_y);
      sfunct.push_back(integrator.value(*this,0.,_y));
    }
    _s27inter = new_ptr(Interpolator<double,double>(sfunct,yvalues,3));
    // compute the normalisation constant
    KaganIntegrand integrand(this);
    _iopt=0;
    _ferminorm*=1./integrator.value(integrand,_MB*(1.-_deltacut)-_mb,_MB-_mb);
    // now for the spectrum
    _mHinter.clear();
    _spectrum.clear();
    _spectmax=ZERO;
    // limits on the mass
    Energy minegamma(0.5*_MB*(1. - _deltacut)),maxegamma(0.5*_MB);
    Energy minhadronmass(max(minMass(),sqrt(_MB*_MB-2.*_MB*maxegamma)));
    Energy maxhadronmass(min(maxMass(),sqrt(_MB*_MB-2.*_MB*minegamma)));
    Energy hstep=(maxhadronmass-minhadronmass)/(_nspect-1);
    Energy mhadron(minhadronmass);
    // function to be integrated
    _iopt=1;
    // prefactor
    InvEnergy2 pre(6.*0.105*2./_MB/_MB*_alpha/pi/semiLeptonicf()*_ckm*_ckm);
    // compute the table
    for(unsigned int ix=0;ix<_nspect;++ix) {
      // calculate y
      _y=1.-mhadron*mhadron/_MB/_MB;
      // perform the integral
      _spectrum.push_back(pre*mhadron*integrator.value(integrand,_MB*_y-_mb,_MB-_mb));
      _spectmax=max(_spectmax,_spectrum.back());
      _mHinter.push_back(mhadron);
      // increment the loop
      mhadron+=hstep;
    }
  }
}

Energy BtoSGammaKagan::hadronicMass(Energy mb,Energy mquark) {
  useMe();
  Energy minmass(max(minMass(),mquark)),maxmass(min(maxMass(),mb)),mass;
  Energy minegamma(0.5*_MB*(1. - _deltacut));//,maxegamma(0.5*_MB);
  minmass=max(minmass,ZERO); // ZERO==sqrt(_MB*_MB-2.*_MB*maxegamma));
  maxmass=min(maxmass,sqrt(_MB*_MB-2.*_MB*minegamma));
  unsigned int ntry(0);
  double rand;
  do {
    rand=UseRandom::rnd();
    mass = minmass*(1.-rand)+rand*maxmass;
    ++ntry;
  }
  while(ntry<_maxtry&&(*_pmHinter)(mass)<UseRandom::rnd()*_spectmax);
  if(ntry>=_maxtry) throw Exception() 
    << "Unweighting failed in BtoSGammaKagan::hadronicMass()" 
    << Exception::eventerror;
  return mass;
}

void BtoSGammaKagan::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BtoSGammaKagan " << name() << " \n";
  output << "newdef " << name() << ":TopMass "    << _mt/GeV << " \n";
  output << "newdef " << name() << ":BottomMass " << _mb/GeV << " \n";
  output << "newdef " << name() << ":CharmMass "  << _mc/GeV << " \n";
  output << "newdef " << name() << ":StrangeMassRatio " << _msovermb << " \n";
  output << "newdef " << name() << ":WMass " << _mw/GeV << " \n";
  output << "newdef " << name() << ":ZMass " << _mz/GeV << " \n";
  output << "newdef " << name() << ":Lambda2 " << _lambda2/GeV2 << " \n";
  output << "newdef " << name() << ":BMesonMass " << _MB/GeV << " \n";
  output << "newdef " << name() << ":Mu " << _mub/GeV << " \n";
  output << "newdef " << name() << ":Delta " << _deltacut << " \n";
  output << "newdef " << name() << ":Lambda1 " << _fermilambda1/GeV2 << " \n";
  output << "newdef " << name() << ":alpha " << _alpha << " \n";
  output << "newdef " << name() << ":CKM " << _ckm << " \n";
  output << "newdef " << name() << ":FermiNormalisation " << _ferminorm*GeV << " \n";
  output << "newdef " << name() << ":MaximumTries " << _maxtry << " \n";
  output << "newdef " << name() << ":ycut " << _ycut << " \n";
  output << "newdef " << name() << ":NSpectrum " <<  _nspect << " \n";
  for(unsigned int ix=0;ix<_mHinter.size();++ix) {
    if(ix<100) output << "newdef " << name() << ":mHValues " << ix << " " 
		      << _mHinter[ix]/GeV << " \n";
    else       output << "insert " << name() << ":mHValues " << ix << " " 
		      << _mHinter[ix]/GeV << " \n";
  }
  for(unsigned int ix=0;ix<_spectrum.size();++ix) {
    if(ix<100) output << "newdef " << name() << ":Spectrum " << ix << " " 
		      << _spectrum[ix]*GeV << " \n";
    else       output << "insert " << name() << ":Spectrum " << ix << " " 
		      << _spectrum[ix]*GeV << " \n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
