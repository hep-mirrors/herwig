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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"

using namespace Herwig;
using namespace ThePEG;

BaryonThreeQuarkModelFormFactor::BaryonThreeQuarkModelFormFactor() 
  : _initialize(false), _order(50),_mlight(420*MeV),_mstrange(570*MeV),
    _LambdaQ(2.5*GeV),_Lambdaqq(0.71*GeV),_Lambdasq(850*MeV),_Lambdass(1.0*GeV) 
{
  // modes handled by this form factor
  // lambda_b
  addFormFactor(5122,4122,2,2,1,2,5,4);
  // xi_b
  addFormFactor(5232,4232,2,2,2,3,5,4);
  addFormFactor(5132,4132,2,2,1,3,5,4);
  // sigma_b
  addFormFactor(5112,4112,2,2,1,1,5,4);
  addFormFactor(5212,4212,2,2,2,1,5,4);
  addFormFactor(5222,4222,2,2,2,2,5,4);
  // omega_b
  addFormFactor(5332,4332,2,2,3,3,5,4);
  // sigma_b-> sigma_c*
  addFormFactor(5112,4114,2,4,1,1,5,4);
  addFormFactor(5212,4214,2,4,2,1,5,4);
  addFormFactor(5222,4224,2,4,2,2,5,4);
  // omega_b -> omega_c*
  addFormFactor(5332,4334,2,4,3,3,5,4);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void BaryonThreeQuarkModelFormFactor::doinit() throw(InitException) {
  BaryonFormFactor::doinit();
  // initialization in needed
  if(_initialize) {
    GaussianIntegrator integrator;
    BaryonCFunction integrand(this);
    _C0.resize(0);_C1.resize(0);_C2.resize(0);
    double pre(0.),root(2.*sqrt(6.));
    double gamma1(1),gamma2(1),gamma3(sqrt(acos(-1.)));
    unsigned int ix,iy;
    for(iy=0;iy<2;++iy) {
      _mu2 = iy==0 ? sqr(_mlight  /_LambdaQ) : sqr(_mstrange/_LambdaQ);
      for(ix=0;ix<=_order;++ix) {
	if(ix>0)    gamma1*=ix;
	if(ix%2==1) gamma2*=(ix+1)/2.0;gamma3*=ix/2.0;
	if(ix%2==0) pre=pow(root,double(ix))/12.*gamma2/gamma1;
	else        pre=pow(root,double(ix))/12.*gamma3/gamma1;
	// for the xi_0 function
	_a=_mu2;_b=2.;_N=ix;
	_C0.push_back(pre*integrator.value(integrand,0.,1.));
	// for the xi_1 function
	_a=_mu2;_b=1.;
	_C1.push_back(pre*integrator.value(integrand,0.,1.));
	// for the xi_2 function
	_a=0.;
	_b=0.;
	_C2.push_back(pre*integrator.value(integrand,0.,1.));
      }
    }
  }
}

void BaryonThreeQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _order << ounit(_mlight,MeV) << ounit(_mstrange,MeV) << ounit(_LambdaQ,MeV) << ounit(_Lambdaqq,MeV) 
     << ounit(_Lambdasq,MeV) << ounit(_Lambdass,MeV) << _C0 << _C1 << _C2;}

void BaryonThreeQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _order >> iunit(_mlight,MeV) >> iunit(_mstrange,MeV) >> iunit(_LambdaQ,MeV) >> iunit(_Lambdaqq,MeV)
     >> iunit(_Lambdasq,MeV) >> iunit(_Lambdass,MeV) >> _C0 >> _C1 >> _C2;}

ClassDescription<BaryonThreeQuarkModelFormFactor> 
BaryonThreeQuarkModelFormFactor::initBaryonThreeQuarkModelFormFactor;
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
			   Complex & f1a,Complex & f2a,Complex & f3a) {
  // this model is based on heavy quark theory
  // therefore most of the factors are zero
  f2v=0.;f2a=0.;f3v=0.;f3a=0.;
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4332) {
    lambdabar=_Lambdass/_LambdaQ;
    ioff=_order+1;
  }
  else if(abs(id1)==4232||abs(id1)==4132||abs(id1)==3322) {
    lambdabar=_Lambdasq/_LambdaQ;
    ioff=_order+1;
  }
  else {
    lambdabar=_Lambdaqq/_LambdaQ;
  }
  // the omega value
  double omega=.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi(phiFunction(omega));
  // use the xi0 functions
  if(abs(id0)==5122||abs(id0)==5232||abs(id0)==5132) {
    double power(1.),numer(0.),denom(0.);
    for(unsigned int ix=0;ix<=_order;++ix) {
      numer+=phi[ix]*power*_C0[ix+ioff];
      denom+=power*_C0[ix+ioff];
      power*=lambdabar;
    }
    f1v=numer/denom;
    f1a=-f1v;
  }
  else {
    double power(1.),numer[2]={0.,0.},denom(0.);
    for(unsigned int ix=0;ix<=_order;++ix) {
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
    f1v = g1v+Complex(0.5*(m0+m1)*(g2v/m0+g3v/m1));
    f1a =-g1a+Complex(0.5*(m0-m1)*(g2a/m0+g3a/m1));
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
				 Complex & f2a,Complex & f3a,Complex & f4a ) {
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4334) {
    lambdabar=_Lambdass/_LambdaQ;
    ioff=_order+1;
  }
  else if(abs(id1)==4234||abs(id1)==4134||abs(id1)==3324) {
    lambdabar=_Lambdasq/_LambdaQ;
    ioff=_order+1;
  }
  else {
    lambdabar=_Lambdaqq/_LambdaQ;
  }
  // the omega value
  double omega=.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi=phiFunction(omega);
  double power(1.),numer[2]={0.,0.},denom(0.);
  for(unsigned int ix=0;ix<=_order;++ix) {
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
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonThreeQuarkModelFormFactor " 
		    << fullName() << " \n";
  output << "set " << fullName() << ":Order       " << _order        << " \n";
  output << "set " << fullName() << ":LightMass   " << _mlight/GeV   << " \n";
  output << "set " << fullName() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "set " << fullName() << ":LambdaQ     " << _LambdaQ/GeV  << " \n";
  output << "set " << fullName() << ":Lambdaqq    " << _Lambdaqq/GeV << " \n";
  output << "set " << fullName() << ":Lambdasq    " << _Lambdasq/GeV << " \n";
  output << "set " << fullName() << ":Lambdass    " << _Lambdass/GeV << " \n";
  // the number of terms to include in the sum for the form-factors
  for(unsigned int ix=0;ix<_C0.size();++ix)
    output << "insert " << fullName() << ":C0 " << ix << "   " << _C0[ix] << " \n";
  for(unsigned int ix=0;ix<_C1.size();++ix)
    output << "insert " << fullName() << ":C1 " << ix << "   " << _C1[ix] << " \n";
  for(unsigned int ix=0;ix<_C2.size();++ix)
    output << "insert " << fullName() << ":C2 " << ix << "   " << _C2[ix] << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
