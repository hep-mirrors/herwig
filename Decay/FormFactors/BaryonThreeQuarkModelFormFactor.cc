// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonThreeQuarkModelFormFactor class.
//

#include "BaryonThreeQuarkModelFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/GaussianIntegrator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

// function for the integral of the partial width
namespace {
using namespace Herwig;

struct BaryonCMatrixElement {

typedef ThePEG::Ptr<BaryonThreeQuarkModelFormFactor>::pointer 
BaryonThreeQuarkModelFormFactorPtr;

BaryonCMatrixElement(BaryonThreeQuarkModelFormFactorPtr in,
		     Energy m0, Energy m1, int type, int mass,
		     int id0,int id1) {
  _formFactor=in;
  _m0=m0;
  _m1=m1;
  _type=type;
  _mass=mass;
  _id0=id0;
  _id1=id1;
}
  
// calculate the integrand  
Energy operator() (double x) const {
  return _formFactor->widthIntegrand(x,_m0,_m1,_type,_mass,_id0,_id1);
} 
  /** Argument type for GaussianIntegrator */
  typedef double ArgType;
  /** Return type for GaussianIntegrator */
  typedef Energy ValType;

private:
  // private variables
  BaryonThreeQuarkModelFormFactorPtr _formFactor;
  Energy _m0,_m1;
  int _type,_mass;
  int _id0,_id1;
};

}

BaryonThreeQuarkModelFormFactor::BaryonThreeQuarkModelFormFactor() 
  : _initialize(false), _order(50),_mlight(420*MeV),_mstrange(570*MeV),
    _LambdaQ(2.5*GeV),_Lambdaqq(0.71*GeV),_Lambdasq(850*MeV),_Lambdass(1.0*GeV) {
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

void BaryonThreeQuarkModelFormFactor::doinit() {
  BaryonFormFactor::doinit();
  // initialization in needed
  if(_initialize) {
    GaussianIntegrator integrator;
    _C0.clear();_C1.clear();_C2.clear();
    double pre(0.),root(2.*sqrt(6.));
    double gamma1(1),gamma2(1),gamma3(sqrt(acos(-1.)));
    unsigned int ix,iy;
    for(iy=0;iy<2;++iy) {
      _mu2 = iy==0 ? sqr(_mlight  /_LambdaQ) : sqr(_mstrange/_LambdaQ);
      for(ix=0;ix<=_order;++ix) {
	if(ix>0)    gamma1*=ix;
	if(ix%2==1) { gamma2*=(ix+1)/2.0; gamma3*=ix/2.0; }
	if(ix%2==0) pre=pow(root,double(ix))/12.*gamma2/gamma1;
	else        pre=pow(root,double(ix))/12.*gamma3/gamma1;
	// for the xi_0 function
	_a=_mu2;_b=2.;_N=ix;
	_C0.push_back(pre*integrator.value(*this,0.,1.));
	// for the xi_1 function
	_a=_mu2;_b=1.;
	_C1.push_back(pre*integrator.value(*this,0.,1.));
	// for the xi_2 function
	_a=0.;
	_b=0.;
	_C2.push_back(pre*integrator.value(*this,0.,1.));
      }
    }
    generator()->log() << "Checking results of BaryonThreeQuarkModelFormFactor"
		       << "vs results from the original paper\n";
    // first matrix element
    Energy m0=getParticleData(5122)->mass();
    Energy m1=getParticleData(4122)->mass();
    double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    BaryonCMatrixElement int1(this,m0,m1,1,0,5122,4122);
    Energy width = integrator.value(int1,1.,omegamax); 
    generator()->log() << "Lambda_b0->Lambda_c+ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // second matrix element
    m0=getParticleData(5222)->mass();
    m1=getParticleData(4222)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,2,0,5222,4222);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Sigma_b+->Sigma_c++ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // third matrix element
    m0=getParticleData(5232)->mass();
    m1=getParticleData(4232)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,1,1,5232,4232);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Xi_b0->Xi_c+ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fourth matrix element
    m0=getParticleData(5332)->mass();
    m1=getParticleData(4332)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,2,2,5332,4332);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Omega_b-->Omega_c0 decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fifth matrix element
    m0=getParticleData(5222)->mass();
    m1=getParticleData(4224)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,3,0,5222,4224);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Sigma_vb+->Sigma_c*++ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fourth matrix element
    m0=getParticleData(5332)->mass();
    m1=getParticleData(4334)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,3,2,5332,4334);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Omega_b-->Omega_c*0 decay" 
		       << width/6.582119E-22/MeV << "\n";
    // output some plots for testing
    double lambdabar = -999.999;
    ofstream output("ThreeQuark.top");
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 3 from paper \"" << endl;
    output << "newdef limits x 1 1.44 y 0.5 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.44),
	step((omegamax-omegamin)/100.),omega(1.),xi;
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer(0.),denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
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
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 6 from paper \"" << endl;
    output << "newdef limits x 1 1.4 y 0.5 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer[2]={0.,0.},denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
	  numer[0]+=phi[iy]*power*_C1[iy+ioff];
	  denom+=power*_C1[iy+ioff];
	  numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	  power*=lambdabar;
	}
	numer[1]/=(omega-1.);
	double xi1(numer[0]/denom);
	output << omega << "   " << xi1 << endl; 
      }
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
    output << "new frame " << endl;
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 7 from paper \"" << endl;
    output << "newdef limits x 1 1.33 y 0.4 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(_order+1);
      if(ix==0){lambdabar=800*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=900*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=1000*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=1050*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=1100*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer[2]={0.,0.},denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
	  numer[0]+=phi[iy]*power*_C1[iy+ioff];
	  denom+=power*_C1[iy+ioff];
	  numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	  power*=lambdabar;
	}
	numer[1]/=(omega-1.);
	double xi1(numer[0]/denom);
	output << omega << "   " << xi1 << endl; 
      }
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
  }
}

void BaryonThreeQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _order << ounit(_mlight,MeV) << ounit(_mstrange,MeV) 
     << ounit(_LambdaQ,MeV) << ounit(_Lambdaqq,MeV) 
     << ounit(_Lambdasq,MeV) << ounit(_Lambdass,MeV) << _C0 << _C1 << _C2;}

void BaryonThreeQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _order >> iunit(_mlight,MeV) >> iunit(_mstrange,MeV) 
     >> iunit(_LambdaQ,MeV) >> iunit(_Lambdaqq,MeV)
     >> iunit(_Lambdasq,MeV) >> iunit(_Lambdass,MeV) >> _C0 >> _C1 >> _C2;}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BaryonThreeQuarkModelFormFactor,BaryonFormFactor>
describeHerwigBaryonThreeQuarkModelFormFactor("Herwig::BaryonThreeQuarkModelFormFactor", "HwFormFactors.so");

void BaryonThreeQuarkModelFormFactor::Init() {

  static ClassDocumentation<BaryonThreeQuarkModelFormFactor> documentation
    ("The BaryonThreeQuarkModelFormFactor class implements"
     " the form-factors for semi-leptonic decay of baryon containing a"
     " heavy quark from PRD56, 348.",
     "The form factors from \\cite{Ivanov:1996fj} were used.",
     "%\\cite{Ivanov:1996fj}\n"
     "\\bibitem{Ivanov:1996fj}\n"
     "  M.~A.~Ivanov, V.~E.~Lyubovitskij, J.~G.~Korner and P.~Kroll,\n"
     "  ``Heavy baryon transitions in a relativistic three-quark model,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 56} (1997) 348\n"
     "  [arXiv:hep-ph/9612463].\n"
     "  %%CITATION = PHRVA,D56,348;%%\n"
     );

  static Parameter<BaryonThreeQuarkModelFormFactor,unsigned int> interfaceOrder
    ("Order",
     "The order of terms to include in the series expansion of the form-factor.",
     &BaryonThreeQuarkModelFormFactor::_order, 10, 0, 1000,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLightMass
    ("LightMass",
     "The mass of the light quark",
     &BaryonThreeQuarkModelFormFactor::_mlight, GeV, .42*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark",
     &BaryonThreeQuarkModelFormFactor::_mstrange, GeV, .57*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaQ
    ("LambdaQ",
     "Heavy Baryon Size Parameter",
     &BaryonThreeQuarkModelFormFactor::_LambdaQ, GeV, 2.5*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaqq
    ("Lambdaqq",
     "The size parameter for light quarks",
     &BaryonThreeQuarkModelFormFactor::_Lambdaqq, GeV, 0.71*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdasq
    ("Lambdasq",
     "The size parameter for one strange quark",
     &BaryonThreeQuarkModelFormFactor::_Lambdasq, GeV, 0.85*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdass
    ("Lambdass",
     "The size parameter with two strange quarks.",
     &BaryonThreeQuarkModelFormFactor::_Lambdass, GeV, 1.0*GeV, ZERO, 2.0*GeV,
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
     "Yes",
     "Perform the initialize",
     true);
  static SwitchOption interfaceInitializeNoInitialize
    (interfaceInitialize,
     "No",
     "No initialization",
     false);
}

void BaryonThreeQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int,int id0,int id1,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a) {
  useMe();
  // this model is based on heavy quark theory
  // therefore most of the factors are zero
  Complex g1v(0.),g1a(0.),g2v(0.),g2a(0.),g3a(0.),g3v(0.);
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4332) {
    lambdabar=_Lambdass/_LambdaQ;
    ioff=_order+1;
  }
  else if(abs(id1)==4232||abs(id1)==4132) {
    lambdabar=_Lambdasq/_LambdaQ;
    ioff=_order+1;
  }
  else {
    lambdabar=_Lambdaqq/_LambdaQ;
  }
  // the omega value
  double omega = 0.5/m0/m1*(m0*m0+m1*m1-q2);
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
    g1v=numer/denom;
    g1a=g1v;
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
    g1v = -(omega*xi1-(sqr(omega)-1.)*xi2)/3.;
    g1a = g1v;
    g2v = 2./3.*(xi1-(omega-1.)*xi2);
    g3v = g2v;
    g2a = 2./3.*(xi1-(omega+1.)*xi2);
    g3a =-g2a;
  }
  // convert to our form
  f1v = g1v+Complex(0.5*(m0+m1)*(g2v/m0+g3v/m1));
  f1a =-g1a+Complex(0.5*(m0-m1)*(g2a/m0+g3a/m1));
  f2v = Complex(0.5*(m0+m1)*( g2v/m0+g3v/m1));
  f3v = Complex(0.5*(m0+m1)*( g2v/m0-g3v/m1));
  f2a =-Complex(0.5*(m0+m1)*( g2a/m0+g3a/m1));
  f3a = Complex(0.5*(m0+m1)*(-g2a/m0+g3a/m1));
}

void  BaryonThreeQuarkModelFormFactor::
SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int,int,int id1,Energy m0,
				Energy m1, Complex & f1v,Complex & f2v,
				Complex & f3v,Complex & f4v,Complex & f1a,
				Complex & f2a,Complex & f3a,Complex & f4a ) {
  useMe();
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
  double omega=0.5/m0/m1*(m0*m0+m1*m1-q2);
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
  // couplings in the velocity form
  Complex g1v,g2v,g3v,g4v,g1a,g2a,g3a,g4a;
  double orr(1./sqrt(3.));
  Energy msum(m0+m1);
  Energy2 msum2(msum*msum);
  g1v = 2.*orr*xi1;
  g1a = -g1v;
  g2v = orr*(xi1-(omega-1)*xi2);
  g2a = orr*(xi1-(omega+1)*xi2);
  g3v = 0.;
  g3a = 0.;
  g4v = -2.*orr*xi2;
  g4a = -g4v;
  // convert to our form
  f1v = g1v;
  f1a =-g1a;
  f2v = g2v*msum/m0;
  f2a =-g2a*msum/m0;
  f3v = Complex(msum2/m0*(g3v/m0+g4v/m1));
  f3a =-Complex(msum2/m0*(g3a/m0+g4a/m1));
  f4v = Complex(msum2/m0/m0*g3v);
  f4a =-Complex(msum2/m0/m0*g3a);
}

void BaryonThreeQuarkModelFormFactor::dataBaseOutput(ofstream & output,bool header,
						     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonThreeQuarkModelFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":Order       " << _order        << " \n";
  output << "newdef " << name() << ":LightMass   " << _mlight/GeV   << " \n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "newdef " << name() << ":LambdaQ     " << _LambdaQ/GeV  << " \n";
  output << "newdef " << name() << ":Lambdaqq    " << _Lambdaqq/GeV << " \n";
  output << "newdef " << name() << ":Lambdasq    " << _Lambdasq/GeV << " \n";
  output << "newdef " << name() << ":Lambdass    " << _Lambdass/GeV << " \n";
  // the number of terms to include in the sum for the form-factors
  for(unsigned int ix=0;ix<_C0.size();++ix)
    output << "insert " << name() << ":C0 " << ix << "   " << _C0[ix] << " \n";
  for(unsigned int ix=0;ix<_C1.size();++ix)
    output << "insert " << name() << ":C1 " << ix << "   " << _C1[ix] << " \n";
  for(unsigned int ix=0;ix<_C2.size();++ix)
    output << "insert " << name() << ":C2 " << ix << "   " << _C2[ix] << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

double BaryonThreeQuarkModelFormFactor::operator ()(double x) const {
  // convert the integration variable
  double y=(1.-x)/x;
  double output =exp(-24.*_mu2*y)*y;
  // the integrals
  double I,Im2;
  SN(y,_N,Im2,I);
  double Nfact=0.5*_N+1.0;
  output *= (_a+Nfact/24./(1.+y))*I+1./12./(1.+y)*(_b-Nfact)*Im2;
  return output;
}

void BaryonThreeQuarkModelFormFactor::SN(double y, int N, double & SNm2,
					 double & SN) const {
  // special cases for the low lying values
  if(N==0) {
    double root=sqrt((1.+y)*(3.+y));
    SN   = 0.5/y*sqrt((1.+y)/(3.+y))*log((root+y)/(root-y));
    SNm2 = 0.5/(y+3.)*(SN+(1.+y)/(3.+4.*y)); 
  }
  else if(N==1) {
    SN   = sqrt(1.+y)/y*asin(y/sqrt(1.+y)/sqrt(3.+y));
    SNm2 = 1./(3.+y)*sqrt((1.+y)/(3.+4.*y)); 
  }
  else if(N==2) {
    double root=sqrt((1.+y)*(3.+y));
    SN   = 1.;
    SNm2 = 0.5/y*sqrt((1.+y)/(3.+y))*log((root+y)/(root-y));
  }
  // the general case
  else {
    int ix; double root;
    if(N%2==0) {
      SN=1.;
      ix=2;
      root=1.;
    }
    else {
      SN = sqrt(1.+y)/y*asin(y/sqrt(1.+y)/sqrt(3.+y));
      ix=1;
      root=sqrt((1.+y)/(3.+4.*y));
    }
    do {
      SNm2=SN;
      ix+=2;
      root*=(3.+4.*y)/(1.+y);
      SN=1.0/(ix-1.0)*(root+(ix-2.0)*(y+3.)*SNm2);
    }
    while(ix<N);
  }
}

// return the phi_N functions calculated using recursion
vector<double> BaryonThreeQuarkModelFormFactor::phiFunction(double omega) {
  vector<double> output;
  double root(sqrt(omega*omega-1.));
  output.push_back(1./root*log(omega+root));
  if(omega<1.00001) output.back()=1.;
  if(_order>0) output.push_back(2./(omega+1.));
  if(_order<2) return output;
  for(unsigned int ix=2;ix<=_order+2;++ix) {
    output.push_back(2./ix/(omega+1.)*(1.+(ix-1)*output[ix-2]));
  }
  return output;
}

Energy BaryonThreeQuarkModelFormFactor::widthIntegrand(double omega,Energy m0,
						       Energy m1, int type,
						       int ,int id0,
						       int id1) {
  // prefactors
  double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  double pi=acos(-1.);
  InvEnergy kw=sqr(generator()->standardModel()->fermiConstant())
    /8./pi/pi/pi*m1*m1*m1/6.*(omegamax-omega)*sqrt(omega*omega-1.);
  Energy2 q2 = sqr(m0)+sqr(m1)-2.*m0*m1*omega;
  if(type<=2) {
    Complex f1v,f2v,f3v,f1a,f2a,f3a;
    SpinHalfSpinHalfFormFactor(q2,0,id0,id1,m0,m1,f1v,f2v,f3v,f1a,f2a,f3a);
    Complex left  =f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
    Complex right =f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
    double g1v = 0.5*( left+right).real();
    double g2v = m0*(f2v+f3v).real()/(m0+m1);
    double g3v = m1*(f2v-f3v).real()/(m0+m1);
    double g1a = -0.5*(+right-left).real();
    double g2a = -m0*(f2a+f3a).real()/(m0+m1);
    double g3a = -m1*(f2a-f3a).real()/(m0+m1);
    Energy Hpp = -2.*sqrt(m0*m1*(omega-1.))*g1v+2.*sqrt(m0*m1*(omega+1))*g1a;
    Energy Hmm = -2.*sqrt(m0*m1*(omega-1.))*g1v-2.*sqrt(m0*m1*(omega+1))*g1a;
    Energy Hp0 = 
      (sqrt(2.*m0*m1*(omega-1))*((m0+m1)*g1v+m1*(omega+1)*g2v+m0*(omega+1)*g3v)-
       sqrt(2.*m0*m1*(omega+1))*((m0-m1)*g1a-m1*(omega-1)*g2a-m0*(omega-1)*g3a))/sqrt(q2);
    Energy Hm0 = 
      (sqrt(2.*m0*m1*(omega-1))*((m0+m1)*g1v+m1*(omega+1)*g2v+m0*(omega+1)*g3v)+
       sqrt(2.*m0*m1*(omega+1))*((m0-m1)*g1a-m1*(omega-1)*g2a-m0*(omega-1)*g3a))/sqrt(q2);
    return kw*sqr(0.04)*(sqr(Hpp)+sqr(Hmm)+sqr(Hp0)+sqr(Hm0));
  }
  else {
    Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
    double  g1v,g2v,g3v,g4v,g1a,g2a,g3a,g4a;
    SpinHalfSpinThreeHalfFormFactor(q2,0,id0,id1,m0,m1,
				    f1v,f2v,f3v,f4v,
				    f1a,f2a,f3a,f4a);
    g1v = f1v.real();
    g1a = -f1a.real();
    g2v = m0/(m0+m1)*f2v.real();
    g2a =-m0/(m0+m1)*f2a.real();
    g3v = sqr(m0/(m0+m1))*f4v.real();
    g3a =-sqr(m0/(m0+m1))*f4a.real();
    g4v = m0*m1/sqr(m0+m1)*(f3v.real()-f4v.real());
    g4a =-m0*m1/sqr(m0+m1)*(f3a.real()-f4a.real());
    Energy HppC = 
      +sqrt(2./3.)*sqrt(m0*m1*(omega-1))*(g1v-2.*(omega+1)*g2v)
      -sqrt(2./3.)*sqrt(m0*m1*(omega+1))*(g1a-2.*(omega-1)*g2a);
    Energy HmmC = 
      -sqrt(2./3.)*sqrt(m0*m1*(omega-1))*(g1v-2.*(omega+1)*g2v)
      -sqrt(2./3.)*sqrt(m0*m1*(omega+1))*(g1a-2.*(omega-1)*g2a);
    Energy HppbC = 
      -sqrt(2.*m0*m1*(omega-1))*g1v
      -sqrt(2.*m0*m1*(omega+1))*g1a;
    Energy HmmbC = 
      sqrt(2.*m0*m1*(omega-1))*g1v
      -sqrt(2.*m0*m1*(omega+1))*g1a;
    Energy Hp0C = (
		   -2./sqrt(3.)*sqrt(m0*m1*(omega-1))*
		   ((m0*omega-m1)*g1v-(m0-m1)*(omega+1)*g2v+m1*(sqr(omega)-1)*g3v+m0*(sqr(omega)-1)*g4v)
		   -2./sqrt(3.)*sqrt(m0*m1*(omega+1))*
		   ((m0*omega-m1)*g1a+(m0+m1)*(omega-1)*g2a+m1*(sqr(omega)-1)*g3a+m0*(sqr(omega)-1)*g4a))/sqrt(q2);
    Energy Hm0C = (
		   +2./sqrt(3.)*sqrt(m0*m1*(omega-1))*
		   ((m0*omega-m1)*g1v-(m0-m1)*(omega+1)*g2v+m1*(sqr(omega)-1)*g3v+m0*(sqr(omega)-1)*g4v)
		   -2./sqrt(3.)*sqrt(m0*m1*(omega+1))*
		   ((m0*omega-m1)*g1a+(m0+m1)*(omega-1)*g2a+m1*(sqr(omega)-1)*g3a+m0*(sqr(omega)-1)*g4a))/sqrt(q2);
    return kw*sqr(0.04)*(sqr(HppC)+sqr(HmmC)+sqr(HppbC)+sqr(HmmbC)+sqr(Hp0C)+sqr(Hm0C));
  }
}
