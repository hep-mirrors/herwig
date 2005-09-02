// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaKagan class.
//

#include "BtoSGammaKagan.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BtoSGammaKagan.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

BtoSGammaKagan::~BtoSGammaKagan() {}

void BtoSGammaKagan::persistentOutput(PersistentOStream & os) const {
  os << _mt << _mb << _mc << _ms << _msovermb << _zratio << _lambda2 << _mw << _mz 
     << _MB << _c20 << _c70 << _c80 << _beta0 << _beta1 << _alpha << _alphaSZ << _mub 
     << _alphaSM << _ckm << _delta << _yinter << _spectrum << _spectmax << _fermilambda 
     << _fermia << _ferminorm << _fermilambda1 <<_ycut << _deltacut << _nsfunct 
     << _nspect << _maxtry << _initialize;
}

void BtoSGammaKagan::persistentInput(PersistentIStream & is, int) {
  is >> _mt >> _mb >> _mc >> _ms >> _msovermb >> _zratio >> _lambda2 >> _mw >> _mz 
     >> _MB >> _c20 >> _c70 >> _c80 >> _beta0 >> _beta1 >> _alpha >> _alphaSZ >> _mub 
     >> _alphaSM >> _ckm >> _delta >> _yinter >> _spectrum >> _spectmax >> _fermilambda 
     >> _fermia >> _ferminorm >> _fermilambda1 >>_ycut >> _deltacut >> _nsfunct 
     >> _nspect >> _maxtry >> _initialize;
}

ClassDescription<BtoSGammaKagan> BtoSGammaKagan::initBtoSGammaKagan;
// Definition of the static class description member.

void BtoSGammaKagan::Init() {

  static ClassDocumentation<BtoSGammaKagan> documentation
    ("The BtoSGammaKagan class implements the calculation of hep-ph/9805303 for the hadronic"
     " mass spectrum in b to s gamma decays");

  static Switch<BtoSGammaKagan,bool> interfaceInitialize
    ("Initialize",
     "Initialize the interpolation tables for the hadronic mass",
     &BtoSGammaKagan::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialize
    (interfaceInitialize,
     "Initialize",
     "Perform initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "NoInitialization",
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
     &BtoSGammaKagan::_lambda2, GeV2, 0.12*GeV2, 0.0*GeV2, 10.0*GeV2,
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

  static ParVector<BtoSGammaKagan,double> interfaceyValues
    ("yValues",
     "The y values for the interpolation of the spectrum",
     &BtoSGammaKagan::_yinter, -1,0., 0.0, 1.,
     false, false, Interface::limited);

  static ParVector<BtoSGammaKagan,double> interfaceSpectrum
    ("Spectrum",
     "Values of the spectrum for interpolation",
     &BtoSGammaKagan::_spectrum,-1, 0., 0., 0,
     false, false, Interface::upperlim);

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

  static Parameter<BtoSGammaKagan,double> interfaceSpectrumMaximum
    ("SpectrumMaximum",
     "The maximum value of the spectrum for unweighting",
     &BtoSGammaKagan::_spectmax, 1.0, 0., 0,
     false, false, Interface::lowerlim);

  static Parameter<BtoSGammaKagan,double> interfaceycut
    ("ycut",
     "Limit of the value of y to avoid singualarities in integrals",
     &BtoSGammaKagan::_ycut, 0.999999999, 0.999, 1.,
     false, false, Interface::limited);

}

void BtoSGammaKagan::calculateWilsonCoefficients()
{
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
  for(unsigned int ix=0;ix<8;++ix)
    {c71c+=(ei[ix]*eta*Ex+fi[ix]+gi[ix]*eta)*pow(eta,ai[ix]);}
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

void BtoSGammaKagan::doinit() throw(InitException) {
  BtoSGammaHadronicMass::doinit();
  cout << "testing start of do init " << endl;
  // quark masses
  _ms=_msovermb*_mb;
  _zratio=pow(_mc/_mb,2);
  // parameters for the fermi motion function
  _fermilambda=_MB-_mb;
  _fermia=-3.*pow(_fermilambda,2.)/_fermilambda1-1.;
  if(_initialize)
    {
      cout << "testing calculating coefficients" << endl;
      // calculate the wilson coefficients etc.
      calculateWilsonCoefficients();
      // calculate the interpolation tables for the s functions
      vector<double> sfunct,yvalues;
      // s_22 function
      sfunct.push_back(0.),yvalues.push_back(0.);
      double step(1./_nsfunct);
      _y=-0.5*step;
      // perform the integrals
      Genfun::AbsFunction * integrand= new KaganIntegrand(this,0);
      GaussianIntegral *integral= new GaussianIntegral(0.,1.); 
      for(unsigned int ix=0;ix<_nsfunct;++ix)
	{
	  _y+=step;
	  integral->resetLimits(0.,_y);
	  yvalues.push_back(_y);
	  sfunct.push_back((*integral)[*integrand]);
	  cout << "testing first integral " << _y << " " << sfunct.back() << endl;
	}
      _s22inter = new Interpolator(sfunct,yvalues,3);
      // s_27 function
      delete integrand;
      integrand=new KaganIntegrand(this,1);
      sfunct.resize(0);yvalues.resize(0);
      sfunct.push_back(0.),yvalues.push_back(0.);
      _y=-0.5*step;
      // perform integrals
      for(unsigned int ix=0;ix<_nsfunct;++ix)
	{
	  _y+=step;
	  integral->resetLimits(0.,_y);
	  yvalues.push_back(_y);
	  sfunct.push_back((*integral)[*integrand]);
	  cout << "testing second integral " << _y << " " << sfunct.back() << endl;
	}
      _s27inter = new Interpolator(sfunct,yvalues,3);
      delete integrand;
      // compute the normalisation constant
      integrand=new KaganIntegrand(this,3);
      integral->resetLimits(_MB*(1.-_deltacut)-_mb,_MB-_mb);
      _ferminorm*=1./(*integral)[*integrand];
      cout << "testing normalisation " << _ferminorm << endl;
      delete integrand;


      // now for the spectrum
      _yinter.resize(0);
      _spectrum.resize(0);
      _spectmax=0.;
      // limits on the mass
      double minegamma(0.5*_MB*(1. - _deltacut)),maxegamma(0.5*_MB);
      double minhadronmass(sqrt(_MB*_MB-2.*_MB*maxegamma));
      double maxhadronmass(sqrt(_MB*_MB-2.*_MB*minegamma));
      double hstep=(maxhadronmass-minhadronmass)/_nspect;
      double mhadron(minhadronmass);
      // function to be integrated      
      integrand=new KaganIntegrand(this,2);
      // prefactor
      double pre(6.*0.105*2./_MB/_MB*_alpha/pi/semiLeptonicf()*_ckm*_ckm);
      cout << "testing prefactor " << _alpha << " " << pre << endl;
      // compute the table
      for(unsigned int ix=0;ix<_nspect+1;++ix)
	{
	  // calculate y
	  _y=1.-mhadron*mhadron/_MB/_MB;
	  // perform the integral
	  integral->resetLimits(_MB*_y-_mb,_MB-_mb);
	  _spectrum.push_back(pre*mhadron*(*integral)[*integrand]);
	  cout << "testing main integral " << mhadron/GeV << "  " << _spectrum.back()*GeV << endl;
	  _spectmax=max(_spectmax,_spectrum.back());
	  _yinter.push_back(_y);
	  // increment the loop
	  mhadron+=hstep;
	}
      // delete the interpolators for the s functions
      delete _s22inter;
      delete _s27inter;    
    }
}

Energy BtoSGammaKagan::hadronicMass(Energy mb,Energy mquark)
{
  Energy minmass(max(minMass(),mquark)),maxmass(min(maxMass(),mb)),mass;
  unsigned int ntry(0);
  double rand;
  do
    {
      rand=generator()->rnd();
      mass = minmass*(1.-rand)+rand*maxmass;
      ++ntry;
    }
  while(ntry<_maxtry&&(*_pyinter)(mass)<generator()->rnd()*_spectmax);
  if(ntry>=_maxtry)
    {throw Exception() << "Unweighting failed in BtoSGammaKagan::hadronicMass()" 
		       << Exception::eventerror;}
  return mass;
}

// function for the integral
namespace Herwig {
using namespace Genfun;
using namespace ThePEG;

FUNCTION_OBJECT_IMP(KaganIntegrand)

  KaganIntegrand::KaganIntegrand(BtoSGammaKaganPtr in,unsigned int ioptin)
{_kagan=in,_iopt=ioptin;}

KaganIntegrand::~KaganIntegrand() {}

KaganIntegrand::KaganIntegrand(const KaganIntegrand & right)
  : _kagan(right._kagan),_iopt(right._iopt) {}

// calculate the integrand  
double KaganIntegrand::operator() (double x) const 
{
  if(_iopt==0)
    {return _kagan->integrands22(x);}
  else if(_iopt==1)
    {return _kagan->integrands27(x);}
  else if(_iopt==2)
    {return _kagan->integrandPy(x);}
  else
    {return _kagan->fermiFunction(x);}
}

}
 
