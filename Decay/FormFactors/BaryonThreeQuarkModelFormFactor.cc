// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonThreeQuarkModelFormFactor class.
//

#include "BaryonThreeQuarkModelFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonThreeQuarkModelFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void BaryonThreeQuarkModelFormFactor::doinit() throw(InitException) {
  BaryonFormFactor::doinit();
  // initialization in needed
  if(_initialize)
    {
      Genfun::AbsFunction * integrand= new BaryonCFunction(this);      
      _C0.resize(0);_C1.resize(0);_C2.resize(0);
      GaussianIntegral *integral= new GaussianIntegral(0.,1.);
      double pre(0.),root(2.*sqrt(6.));
      double gamma1(1),gamma2(1),gamma3(sqrt(acos(-1.)));
      unsigned int ix,iy;
      for(iy=0;iy<2;++iy)
	{
	  if(iy==0){_mu2=_mlight*_mlight/_LambdaQ/_LambdaQ;}
	  else{_mu2=_mstrange*_mstrange/_LambdaQ/_LambdaQ;}
	  for(ix=0;ix<=_order;++ix)
	    {
	      if(ix>0){gamma1*=ix;}
	      if(ix%2==1){gamma2*=(ix+1)/2.0;gamma3*=ix/2.0;}
	      if(ix%2==0){pre=pow(root,ix)/12.*gamma2/gamma1;}
	      else{pre=pow(root,ix)/12.*gamma3/gamma1;}
	      // for the xi_0 function
	      _a=_mu2;_b=2.;_N=ix;
	      _C0.push_back(pre*(*integral)[*integrand]);
	      // for the xi_1 function
	      _a=_mu2;_b=1.;
	      _C1.push_back(pre*(*integral)[*integrand]);
	      // for the xi_2 function
	      _a=0.;
	      _b=0.;
	      _C2.push_back(pre*(*integral)[*integrand]);
	    }
	}
      // tidy up
      delete integrand;
    }
}

BaryonThreeQuarkModelFormFactor::~BaryonThreeQuarkModelFormFactor() {}

void BaryonThreeQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _order << _mlight << _mstrange << _LambdaQ << _Lambdaqq 
     << _Lambdasq << _Lambdass << _C0 << _C1 << _C2;}

void BaryonThreeQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _order >> _mlight >> _mstrange >> _LambdaQ >> _Lambdaqq
     >> _Lambdasq >> _Lambdass >> _C0 >> _C1 >> _C2;}

ClassDescription<BaryonThreeQuarkModelFormFactor> BaryonThreeQuarkModelFormFactor::initBaryonThreeQuarkModelFormFactor;
// Definition of the static class description member.

void BaryonThreeQuarkModelFormFactor::Init() {

  static ClassDocumentation<BaryonThreeQuarkModelFormFactor> documentation
    ("The BaryonThreeQuarkModelFormFactor class implements"
     " the form-factors for semi-leptonic decay of baryon containing a"
     " heavy quark from PRD56, 348.");

  static Parameter<BaryonThreeQuarkModelFormFactor,unsigned int> interfaceOrder
    ("Order",
     "The order of terms to include in the series expansion of the form-factor.",
     &BaryonThreeQuarkModelFormFactor::_order, 10, 0, 1000,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLightMass
    ("LightMass",
     "The mass of the light quark",
     &BaryonThreeQuarkModelFormFactor::_mlight, GeV, .42*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark",
     &BaryonThreeQuarkModelFormFactor::_mstrange, GeV, .57*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaQ
    ("LambdaQ",
     "Heavy Baryon Size Parameter",
     &BaryonThreeQuarkModelFormFactor::_LambdaQ, GeV, 2.5*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaqq
    ("Lambdaqq",
     "The size parameter for light quarks",
     &BaryonThreeQuarkModelFormFactor::_Lambdaqq, GeV, 0.71*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdasq
    ("Lambdasq",
     "The size parameter for one strange quark",
     &BaryonThreeQuarkModelFormFactor::_Lambdasq, GeV, 0.85*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdass
    ("Lambdass",
     "The size parameter with two strange quarks.",
     &BaryonThreeQuarkModelFormFactor::_Lambdass, GeV, 1.0*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC0
    ("C0",
     "The coefficient of the expansion for xi_0.",
     &BaryonThreeQuarkModelFormFactor::_C0,
     0, 0, 0, -1E20, 1E20, false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC1
    ("C1",
     "The coefficient of the expansion for xi_1.",
     &BaryonThreeQuarkModelFormFactor::_C1,
     0, 0, 0, -1E20, 1E20, false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC2
    ("C2",
     "The coefficient of the expansion for xi_2.",
     &BaryonThreeQuarkModelFormFactor::_C2,
     0, 0, 0, -1E20, 1E20, false, false, true);


  static Switch<BaryonThreeQuarkModelFormFactor,bool> interfaceInitialize
    ("Initialize",
     "Initialize the coefficient for the expansion of the form-factor",
     &BaryonThreeQuarkModelFormFactor::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialize
    (interfaceInitialize,
     "Initialize",
     "Perform the initialize",
     true);
  static SwitchOption interfaceInitializeNoInitialize
    (interfaceInitialize,
     "NoInitialize",
     "No initialization",
     false);
}

void BaryonThreeQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int,int id0,int id1,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  // this model is based on heavy quark theory
  // therefore most of the factors are zero
  f2v=0.;f2a=0.;f3v=0.;f3a=0.;
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4332)
    {lambdabar=_Lambdass/_LambdaQ;ioff=_order+1;}
  else if(abs(id1)==4232||abs(id1)==4132||abs(id1)==3322)
    {lambdabar=_Lambdasq/_LambdaQ;ioff=_order+1;}
  else{lambdabar=_Lambdaqq/_LambdaQ;}
  // the omega value
  double omega=.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi(phiFunction(omega));
  // use the xi0 functions
  if(abs(id0)==5122||abs(id0)==5232||abs(id0)==5132)
    {
      double power(1.),numer(0.),denom(0.);
      for(unsigned int ix=0;ix<=_order;++ix)
	{
	  numer+=phi[ix]*power*_C0[ix+ioff];
	  denom+=power*_C0[ix+ioff];
	  power*=lambdabar;
	}
      f1v=numer/denom;
      f1a=-f1v;
    }
  else
    {
      double power(1.),numer[2]={0.,0.},denom(0.);
      for(unsigned int ix=0;ix<=_order;++ix)
	{
	  numer[0]+=phi[ix]*power*_C1[ix+ioff];
	  denom+=power*_C1[ix+ioff];
	  numer[1]+=power*_C2[ix+ioff]*(phi[ix]-phi[ix+2]);
	  power*=lambdabar;
	}
      numer[1]/=(omega-1.);
      double xi1(numer[0]/denom),xi2(numer[1]/denom);
      // the couplings in the velocity form
      Complex g1v,g1a,g2v,g2a,g3a,g3v;
      g1v = -(omega*xi1-(omega*omega-1.)*xi2)/3.;
      g1a = g1v;
      g2v = 2./3.*(xi1-(omega-1.)*xi2);
      g3v = g2v;
      g2a = 2./3.*(xi1-(omega+1.)*xi2);
      g3a =-g2a;
      // convert to our form
      f1v = g1v+0.5*(m0+m1)*(g2v/m0+g3v/m1);
      f1a =-g1a+0.5*(m0-m1)*(g2a/m0+g3a/m1);
      f2v = 0.5*(m0+m1)*( g2v/m0+g3v/m1);
      f3v = 0.5*(m0+m1)*( g2v/m0-g3v/m1);
      f2a =-0.5*(m0+m1)*( g2a/m0+g3a/m1);
      f3a = 0.5*(m0+m1)*(-g2a/m0+g3a/m1);
    }
}

void  BaryonThreeQuarkModelFormFactor::
 SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int,int,int id1,Energy m0,
				 Energy m1, Complex & f1v,Complex & f2v,
				 Complex & f3v,Complex & f4v,Complex & f1a,
				 Complex & f2a,Complex & f3a,Complex & f4a )
{
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4334)
    {lambdabar=_Lambdass/_LambdaQ;ioff=_order+1;}
  else if(abs(id1)==4234||abs(id1)==4134||abs(id1)==3324)
    {lambdabar=_Lambdasq/_LambdaQ;ioff=_order+1;}
  else{lambdabar=_Lambdaqq/_LambdaQ;}
  // the omega value
  double omega=.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi=phiFunction(omega);
  double power(1.),numer[2]={0.,0.},denom(0.);
  for(unsigned int ix=0;ix<=_order;++ix)
    {
      numer[0]+=phi[ix]*power*_C1[ix+ioff];
      denom+=power*_C1[ix+ioff];
      numer[1]+=power*_C2[ix+ioff]*(phi[ix]-phi[ix+2]);
      power*=lambdabar;
    }
  numer[1]/=(omega-1.);
  double xi1=numer[0]/denom;
  double xi2=numer[1]/denom;
  // the couplings in the velocity form
  Complex N1,N2,N3,N4,K1,K2,K3,K4; 
  double orr(1./sqrt(3.));
  Energy msum(m0+m1);Energy2 msum2(msum*msum);
  // in the form of the heavy quark papers
  N1 = orr*(xi1-(omega-1.)*xi2);
  K1 = orr*(xi1-(omega+1.)*xi2);
  N2 = 0.;
  K2 = 0.;
  N3 = -2.*orr*xi2;
  K3 = -N3;
  N4 = 2.*orr*xi1;
  K4 = -2*orr*xi1;
  // convert to our form
  f1v = N4;
  f1a =-K4;
  f2v = N1*msum/m0;
  f2a =-K1*msum/m0;
  f3v = msum2/m0*(N2/m0+N3/m1);
  f3a =-msum2/m0*(K2/m0+K3/m1);
  f4v = msum2/m0/m0*N2;
  f4a =-msum2/m0/m0*K2;
}

void BaryonThreeQuarkModelFormFactor::dataBaseOutput(ofstream & output,bool header,
						     bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create /Herwig++/BaryonThreeQuarkModelFormFactor " 
	    << fullName() << " \n";}
  output << "set " << fullName() << ":Order       " << _order        << " \n";
  output << "set " << fullName() << ":LightMass   " << _mlight/GeV   << " \n";
  output << "set " << fullName() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "set " << fullName() << ":LambdaQ     " << _LambdaQ/GeV  << " \n";
  output << "set " << fullName() << ":Lambdaqq    " << _Lambdaqq/GeV << " \n";
  output << "set " << fullName() << ":Lambdasq    " << _Lambdasq/GeV << " \n";
  output << "set " << fullName() << ":Lambdass    " << _Lambdass/GeV << " \n";
  // the number of terms to include in the sum for the form-factors
  for(unsigned int ix=0;ix<_C0.size();++ix)
    {output << "insert " << fullName() << ":C0 " << ix << "   " << _C0[ix] << " \n";}
  for(unsigned int ix=0;ix<_C1.size();++ix)
    {output << "insert " << fullName() << ":C1 " << ix << "   " << _C1[ix] << " \n";}
  for(unsigned int ix=0;ix<_C2.size();++ix)
    {output << "insert " << fullName() << ":C2 " << ix << "   " << _C2[ix] << " \n";}
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}


// function for the integral
namespace Herwig {
using namespace Genfun;
using namespace ThePEG;

FUNCTION_OBJECT_IMP(BaryonCFunction)

BaryonCFunction::BaryonCFunction(BaryonThreeQuarkModelFormFactorPtr in)
{_formFactor=in;}

// calculate the integrand  
double BaryonCFunction::operator() (double x) const 
{return _formFactor->integrandC(x);}

}


/*

// from .h


  //
   *  The integrand for the semi-analytic calculation of the semi-leptonic width.
   *  This is mainly included for testing purposes.
   * @param omega The \f$\omega\f$ parameter of the heavy quark form-factors.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param type The type of the decay 
   * @param imass The baryonic mass parameter to use.
   //
  inline double widthIntegrand(double omega,Energy m0, Energy m1, int type, int imass);

// from .icc



  cout << "testing the paper form of the me" << endl;
  // first matrix element
  Energy m0=getParticleData(5122)->mass();
  Energy m1=getParticleData(4122)->mass();
  double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  Genfun::AbsFunction *me1=new BaryonCMatrixElement(this,m0,m1,1,0);
  GaussianIntegral *integral= new GaussianIntegral(1.,omegamax);
  cout << "testing lambda decay" << (*integral)[*me1]/6.582119E-22 << endl;
  delete me1; delete integral;
  // second matrix element
  m0=getParticleData(5222)->mass();
  m1=getParticleData(4222)->mass();
  omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  me1=new BaryonCMatrixElement(this,m0,m1,2,0);
  integral= new GaussianIntegral(1.,omegamax);
  cout << "testing sigma decay" << (*integral)[*me1]/6.582119E-22 << endl;
  delete me1; delete integral;
  // third matrix element
  m0=getParticleData(5232)->mass();
  m1=getParticleData(4232)->mass();
  omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  me1=new BaryonCMatrixElement(this,m0,m1,1,1);
  integral= new GaussianIntegral(1.,omegamax);
  cout << "testing xi decay" << (*integral)[*me1]/6.582119E-22 << endl;
  delete me1; delete integral;
  // fourth matrix element
  m0=getParticleData(5332)->mass();
  m1=getParticleData(4332)->mass();
  omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  me1=new BaryonCMatrixElement(this,m0,m1,2,2);
  integral= new GaussianIntegral(1.,omegamax);
  cout << "testing omega decay" << (*integral)[*me1]/6.582119E-22 << endl;
  delete me1; delete integral;
  // fifth matrix element
  m0=getParticleData(5222)->mass();
  m1=getParticleData(4224)->mass();
  omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  me1=new BaryonCMatrixElement(this,m0,m1,3,0);
  integral= new GaussianIntegral(1.,omegamax);
  cout << "testing sigma decay" << (*integral)[*me1]/6.582119E-22 << endl;
  delete me1; delete integral;
  // fourth matrix element
  m0=getParticleData(5332)->mass();
  m1=getParticleData(4334)->mass();
  omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  me1=new BaryonCMatrixElement(this,m0,m1,3,2);
  integral= new GaussianIntegral(1.,omegamax);
  cout << "testing omega decay" << (*integral)[*me1]/6.582119E-22 << endl;



inline double BaryonThreeQuarkModelFormFactor::widthIntegrand(double omega,Energy m0,
							      Energy m1, int type,
							      int imass)
 {
   InvEnergy2 GF=1.16639E-5/GeV2;
   double vcb=1.;
   // the IW function
   vector<double> phi=phiFunction(omega);
   double Hpp(0.),Hmm(0.),Hp0(0.),Hm0(0.),Hppb(0.),Hmmb(0.);
   double HA3(0.),HV3(0.),HA1(0.),HV1(0.),HA0(0.),HV0(0.);
   double lambdabar,xit;
   if(imass==0){lambdabar=_Lambdaqq/_LambdaQ;}
   else if(imass==1){lambdabar=_Lambdasq/_LambdaQ;}
   else{lambdabar=_Lambdass/_LambdaQ;}
   if(type==1)
     {
       double power(1.),numer(0.),denom(0.);
       unsigned int iy;
       for(unsigned int ix=0;ix<=_order;++ix)
	 {
	   if(imass>0){iy=ix+_order+1;}
	   else{iy=ix;}
	   numer+=phi[ix]*power*_C0[iy];
	   denom+=power*_C0[iy];
	   power*=lambdabar;
	 }
       double xi0(numer/denom);
       Hpp = xi0;
       Hmm = xi0;
       Hp0 = xi0*((m0+m1)*sqrt(omega-1.)-(m0-m1)*sqrt(omega+1.));
       Hm0 = xi0*((m0+m1)*sqrt(omega-1.)+(m0-m1)*sqrt(omega+1.));
     }
   else if(type==2)
     {
       double power(1.),numer[2]={0.,0.},denom(0.);
       unsigned int iy;
       for(unsigned int ix=0;ix<=_order;++ix)
	 {
	   if(imass>0){iy=ix+_order+1;}
	   else{iy=ix;}
	   numer[0]+=phi[ix]*power*_C1[iy];
	   denom+=power*_C1[iy];
	   numer[1]+=power*_C2[iy]*(phi[ix]-phi[ix+2]);
	   power*=lambdabar;
	 }
       numer[1]/=(omega-1.);
       double xi1=numer[0]/denom;
       double xi2=numer[1]/denom;
       xit = xi1*omega-xi2*(omega*omega-1.);
       double xiLp = xi1*(omega+2.)-xi2*(omega*omega-1.);
       double xiLm = xi1*(omega-2.)-xi2*(omega*omega-1.);
       //double xiLp(xit),xiLm(xit);
       Hpp = 1./3.*xit;
       Hmm = 1./3.*xit;
       Hp0 = 1./3.*((m0+m1)*sqrt(omega-1.)*xiLp-(m0-m1)*xiLm*sqrt(omega+1.));
       Hm0 = 1./3.*((m0+m1)*sqrt(omega-1.)*xiLp+(m0-m1)*xiLm*sqrt(omega+1.));
     }
   else
     {
       double power(1.),numer[2]={0.,0.},denom(0.);
       unsigned int iy;
       for(unsigned int ix=0;ix<=_order;++ix)
	 {
	   if(imass>0){iy=ix+_order+1;}
	   else{iy=ix;}
	   numer[0]+=phi[ix]*power*_C1[iy];
	   denom+=power*_C1[iy];
	   numer[1]+=power*_C2[iy]*(phi[ix]-phi[ix+2]);
	   power*=lambdabar;
	 }
       numer[1]/=(omega-1.);
       double xi1=numer[0]/denom;
       double xi2=numer[1]/denom;
       // the form factors
       double N1(0.),N2(0.),N3(0.),N4(0.),K1(0.),K2(0.),K3(0.),K4(0.),orr(1./sqrt(3.));
       N1 = orr*(xi1-(omega-1.)*xi2);
       K1 = orr*(xi1-(omega+1.)*xi2);
       N2 = 0.;
       K2 = 0.;
       N3 =-2.*orr*xi2;
       K3 =-N3;
       N4 = 2.*orr*xi1;
       K4 =-2.*orr*xi1;
       // the coefficients
       xit  = xi1*omega-xi2*(omega*omega-1.);
       double xiLps = xi1*(omega-1.)-xi2*(omega*omega-1.);
       double xiLms = xi1*(omega+1.)-xi2*(omega*omega-1.);
       Hpp  = sqrt(2.)/3.*xit;
       Hmm  =-sqrt(2.)/3.*xit;
       Hppb = sqrt(2./3.)*xi1;
       Hmmb =-sqrt(2./3.)*xi1;
       Hp0  = sqrt(2.)/3.*((m0+m1)*sqrt(omega-1.)*xiLps-(m0-m1)*xiLms*sqrt(omega+1.));
       Hm0  = sqrt(2.)/3.*((m0+m1)*sqrt(omega-1.)*xiLps+(m0-m1)*xiLms*sqrt(omega+1.));
       HA3 = sqrt(2.*m0*m1*(omega+1.))*K4;
       HV3 =-sqrt(2.*m0*m1*(omega-1.))*N4;
       HV1 = sqrt(2./3.)*sqrt(m0*m1*(omega-1.))*(N4-2.*(omega+1.)*N1);
       HA1 = sqrt(2./3.)*sqrt(m0*m1*(omega+1.))*(K4-2.*(omega-1.)*K1);
       HV0 = -2./sqrt(3.)*sqrt(m0*m1*(omega-1.))*((m0*omega-m1)*N4-(m0-m1)*(omega+1.)*N1
						  +m1*(omega*omega-1.)*N2
						  +m0*(omega*omega-1.)*N3);
       HA0 =  2./sqrt(3.)*sqrt(m0*m1*(omega+1.))*((m0*omega-m1)*K4+(m0+m1)*(omega-1.)*K1
						  +m1*(omega*omega-1.)*K2
						  +m0*(omega*omega-1.)*K3);
     }
   // prefactors
   double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
   Hpp  *= -2.*sqrt(m0*m1)*(sqrt(omega-1.)-sqrt(omega+1.));
   Hmm  *= -2.*sqrt(m0*m1)*(sqrt(omega-1.)+sqrt(omega+1.));
   Hppb *= -2.*sqrt(m0*m1)*(sqrt(omega-1.)-sqrt(omega+1.));
   Hmmb *= -2.*sqrt(m0*m1)*(sqrt(omega-1.)+sqrt(omega+1.));
   Hp0  *= 1./sqrt(omegamax-omega);
   Hm0  *= 1./sqrt(omegamax-omega);
   HV0 *=1./sqrt(2.*m0*m1*(omegamax-omega));
   HA0 *=1./sqrt(2.*m0*m1*(omegamax-omega));
   double pi=acos(-1.);
   double kw=GF*GF*vcb*vcb/8./pi/pi/pi*m1*m1*m1/6.*
     (omegamax-omega)*sqrt(omega*omega-1.);
   double output;
   if(type<=2)
     {output=kw*(Hpp*Hpp+Hmm*Hmm+Hp0*Hp0+Hm0*Hm0+Hppb*Hppb+Hmmb*Hmmb);}
   else
     {output=2.*kw*(HA3*HA3 +HV3*HV3+HA0*HA0+HV0*HV0+HA1*HA1+HV1*HV1);}
   return output;
}


// from .cc

  // output some plots for testing
  double lambdabar;
  ofstream output("ThreeQuark.top");
  output << "set font duplex" << endl;
  output << "title top \"Figure 3 from paper \"" << endl;
  output << "set limits x 1 1.44 y 0.5 1" << endl;
  for(unsigned int ix=0;ix<5;++ix)
    {
      double omegamin(1.),omegamax(1.44),step((omegamax-omegamin)/100.),omega(1.),xi;
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step)
	{
	  vector<double> phi(phiFunction(omega));
	  double power(1.),numer(0.),denom(0.);
	  for(unsigned int iy=0;iy<=_order;++iy)
	    {
	      numer+=phi[iy]*power*_C0[iy+ioff];
	      denom+=power*_C0[iy+ioff];
	      power*=lambdabar;
	    }
	  xi=numer/denom;
	  output << omega << "   " << xi << endl; 
	}
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
  output << "new frame " << endl;
  output << "set font duplex" << endl;
  output << "title top \"Figure 6 from paper \"" << endl;
  output << "set limits x 1 1.4 y 0.5 1" << endl;
  for(unsigned int ix=0;ix<5;++ix)
    {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step)
	{
	  vector<double> phi(phiFunction(omega));
	  double power(1.),numer[2]={0.,0.},denom(0.);
	  for(unsigned int iy=0;iy<=_order;++iy)
	    {
	      numer[0]+=phi[iy]*power*_C1[iy+ioff];
	      denom+=power*_C1[iy+ioff];
	      numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	      power*=lambdabar;
	    }
	  numer[1]/=(omega-1.);
	  double xi1(numer[0]/denom),xi2(numer[1]/denom);
	  output << omega << "   " << xi1 << endl; 
	}
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
  output << "new frame " << endl;
  output << "set font duplex" << endl;
  output << "title top \"Figure 7 from paper \"" << endl;
  output << "set limits x 1 1.33 y 0.4 1" << endl;
  for(unsigned int ix=0;ix<5;++ix)
    {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(_order+1);
      if(ix==0){lambdabar=800*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=900*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=1000*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=1050*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=1100*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step)
	{
	  vector<double> phi(phiFunction(omega));
	  double power(1.),numer[2]={0.,0.},denom(0.);
	  for(unsigned int iy=0;iy<=_order;++iy)
	    {
	      numer[0]+=phi[iy]*power*_C1[iy+ioff];
	      denom+=power*_C1[iy+ioff];
	      numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	      power*=lambdabar;
	    }
	  numer[1]/=(omega-1.);
	  double xi1(numer[0]/denom),xi2(numer[1]/denom);
	  output << omega << "   " << xi1 << endl; 
	}
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }


// function for the integral
namespace Herwig {
using namespace Genfun;
using namespace ThePEG;

FUNCTION_OBJECT_IMP(BaryonCMatrixElement)

  BaryonCMatrixElement::BaryonCMatrixElement(BaryonThreeQuarkModelFormFactorPtr in,
					     Energy m0, Energy m1, int type, int mass)
 {
   _formFactor=in;
   _m0=m0;
   _m1=m1;
   _type=type;
   _mass=mass;
 }

BaryonCMatrixElement::~BaryonCMatrixElement() {}

BaryonCMatrixElement::BaryonCMatrixElement(const BaryonCMatrixElement & right)
  : _formFactor(right._formFactor), _m0(right._m0), _m1(right._m1), _type(right._type),
    _mass(right._mass) {}

// calculate the integrand  
double BaryonCMatrixElement::operator() (double x) const 
 {return _formFactor->widthIntegrand(x,_m0,_m1,_type,_mass);}

}




 */
