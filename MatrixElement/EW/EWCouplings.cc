// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EWCouplings class.
//

#include "EWCouplings.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <boost/numeric/ublas/operation.hpp>

using namespace Herwig;

namespace {

Complex trace(boost::numeric::ublas::matrix<Complex> M) {
  assert(M.size1()==M.size2());
  Complex output(0.);
  for(unsigned int ix=0;ix<M.size1();++ix)
    output += M(ix,ix);
  return output;
}
}


EWCouplings::EWCouplings(unsigned int loops, unsigned int steps, Energy highScale,
			 Energy lowScale) 
  : ewScale_(91.1876*GeV), highScale_(highScale), lowScale_(lowScale),
    includeSU3_(true), includeEW_(true), initialized_(false), massChoice_(false),
    mZ_(91.1876*GeV), mW_(80.399*GeV), 
    mT_(173.1*GeV), // 179.08045 (should be this?)
    loops_(loops), highSteps_(steps), lowSteps_(steps)
{}

void EWCouplings::initialize() {
  using Constants::pi;
  if(initialized_) return;
  // set the particle masses
  if(massChoice_) {
    mZ_ = getParticleData(ParticleID::Z0   )->mass();
    mW_ = getParticleData(ParticleID::Wplus)->mass();
    mT_ = getParticleData(ParticleID::t    )->mass();
    mH_ = getParticleData(ParticleID::h0   )->mass();
  }
  // logs of scales
  double logEWScale   = log(ewScale_/GeV);
  double logHighScale = log(highScale_/GeV);
  // step size
  double stepsize = (logHighScale - logEWScale)/double(highSteps_);
  // Initialize parameters at the ewScale
  // 32 parameters, mostly zero due massless quarks
  unsigned int N = 32;
  vector<Complex> y(N,0.), dydx(N,0.), yout(N,0.);
  initializeCouplings(y);
  double x = logEWScale;
  derivatives(x,y,dydx);
  // energy scale + 6 parameters: g1,g2,g3,y_t,lambda,vev
  table_ = boost::numeric::ublas::matrix<double>(highSteps_+1,7);
  table_(0,0) = logEWScale;
  for (unsigned int i=1; i<N; i++) {
    if (i <=3)  table_(0,i   ) = (y[i-1].real()*y[i-1].real())/(4.0*pi);
    if (i >=29) table_(0,i-25) = y[i].real();
  }
  int counter = 1;
   
  // Use 4th order runge-kutta to integrate to highScale
  int steps = highSteps_;
  while (steps > 0) {
//      _RK4(y,dydx,x,stepsize,yout,_loops);

    // Advance x and calculate derivatives at new starting point
    for(unsigned int j=0; j<N; j++) {
      y[j] = yout[j];  
    }
    x += stepsize;
    derivatives(x,y,dydx);
     
    table_(counter,0) = x;
    for (unsigned int i=1; i<N; i++) {
      if (i<=3 ) table_(counter,i   ) = (y[i-1].real()*y[i-1].real())/(4.0*pi);
      if (i>=29) table_(counter,i-25) = y[i].real();
    }
     
    steps--;
    counter++;
  }
   
//    _LowInit(); // Initialize couplings at mu < 91.1876 GeV
}

EWCouplings::~EWCouplings() {}

IBPtr EWCouplings::clone() const {
  return new_ptr(*this);
}

IBPtr EWCouplings::fullclone() const {
  return new_ptr(*this);
}

void EWCouplings::persistentOutput(PersistentOStream & os) const {
  os << ounit(ewScale_,GeV) <<  ounit(highScale_,GeV) <<  ounit(lowScale_,GeV)
     << includeSU3_ << includeEW_ << ounit(mZ_,GeV) << ounit(mW_,GeV)
     << ounit(mT_,GeV) << ounit(mH_,GeV) << massChoice_ << initialized_;
}

void EWCouplings::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ewScale_,GeV) >>  iunit(highScale_,GeV) >> iunit(lowScale_,GeV)
     >> includeSU3_ >> includeEW_ >> iunit(mZ_,GeV) >> iunit(mW_,GeV) 
     >> iunit(mT_,GeV) >> iunit(mH_,GeV) >> massChoice_ >> initialized_;
}


//The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EWCouplings,Interfaced>
describeHerwigEWCouplings("Herwig::EWCouplings", "HwMEEW.so");

void EWCouplings::Init() {

  static ClassDocumentation<EWCouplings> documentation
    ("The EWCouplings class implements");

  static Switch<EWCouplings,bool> interfaceMassChoice
    ("MassChoice",
     "Where to get the SM particle masses from",
     &EWCouplings::massChoice_, false, false, false);
  static SwitchOption interfaceMassChoiceLocal
    (interfaceMassChoice,
     "Local",
     "Use local values",
     false);
  static SwitchOption interfaceMassChoiceParticleData
    (interfaceMassChoice,
     "ParticleData",
     "Get the values from the ParticleData object",
     true);

}



void EWCouplings::initializeLow() {
  using Constants::pi;
  // For scales less than ewScale, the only couplings calculated here are those of
  // alpha_EW and alpha3:
   
  double logEWScale  = log(ewScale_ /GeV);
  double loglowScale = log(lowScale_/GeV);
   
  double stepsize = (loglowScale - logEWScale)/double(lowSteps_);
  int steps = lowSteps_;
   
  // Initialize parameters at the ewScale
  unsigned int N=2; // Total # of parameters = 2
  vector<Complex> y(N), dydx(N), yout(N);
   for (unsigned int i=0; i<N; i++) {
      y[i] = 0.0;
      dydx[i] = 0.0;
      yout[i] = 0.0;
   }
   double x = logEWScale;
   
   // Initialize Couplings at the ewScale, including the One-Loop Threshold Matching:
   double a1 = (3.0/5.0)*table_(0,1);
   double a2 = table_(0,2);
   double a3 = table_(0,3);
   double aEM_ewScale = 1./(1./a1+1./a2-1./(6.*pi)*(1.-21.*log(mW()/ewScale_)+16./3.*log(mTatmZ()/ewScale_)));
   double aS_ewScale = 1./(1./a3+1./(3.*pi)*log(ewScale_/mTatmZ()));
   y[0] = sqrt(4.0*pi*aEM_ewScale);
   y[1] = sqrt(4.0*pi*aS_ewScale);
   
   derivatives(x,y,dydx);
   // energy scale + 2 parameters: gEM,g3
   lowTable_.resize(lowSteps_+1,3);
   lowTable_(0,0) = logEWScale;
   for (unsigned int i=1; i<=N; i++) {
      lowTable_(0,i) = (y[i-1].real()*y[i-1].real())/(4.0*pi);
   }
   int counter = 1;
   
   // Use 4th order runge-kutta to integrate to highScale
   while (steps > 0) {
     RK4(y,dydx,x,stepsize,yout);
     // Advance x and calculate derivatives at new starting point
     for(unsigned int j=0; j<N; j++) {
       y[j] = yout[j];  
     }
     x += stepsize;
     derivatives(x,y,dydx);
		
     lowTable_(counter,0) = x;
     for (unsigned int i=1; i<=N; i++) {
       lowTable_(counter,i) = (y[i-1].real()*y[i-1].real())/(4.0*pi);
     }

      steps--;
      counter++;
   }
}

void EWCouplings::RK4(vector<Complex> &y, vector<Complex> &dydx,
		      const double x, const double h, vector<Complex> &yout) {
   
  unsigned int n = y.size();
  std::vector<Complex> dym(n), dyt(n), yt(n);
  double hh = h*0.5;
  double h6 = h/6.0;
  double xh = x + hh;
  const Complex I(0,1.0);
	
  for(unsigned int i=0; i<n; i++) yt[i] = y[i] + hh*dydx[i];
  derivatives(xh, yt, dyt);
  for(unsigned int i=0; i<n; i++) yt[i] = y[i] + hh*dyt[i];
  derivatives(xh, yt, dym);
  for(unsigned int i=0; i<n; i++) {
    yt[i] = y[i] + h*dym[i];
    dym[i] += dyt[i];
  }
  derivatives(x+h, yt, dyt);
  for(unsigned int i=0; i<n; i++) {
    yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);
  }
}


void EWCouplings::initializeCouplings(vector<Complex> & y) {
  // \todo make these values parameters so they can be changed
  InvEnergy2 gFermi = 1.16637*pow(10.0,-5)/GeV2;
  Energy vev = 1.0/(sqrt(sqrt(2.0)*gFermi)); // vev = 246.221
  
  y[0] = 0.461531463;     // g1 = Sqrt[5/3] * Sqrt[4*pi*a1] with a1(Mz) = 0.01017054
  y[1] = 0.651547066;     // g2 = Sqrt[4*pi*a2] with a2(Mz) = 0.03378168
  y[2] = 1.215650108;     // g3 = Sqrt[4*pi*as] with as(Mz) = 0.1176
  
  // Note lambda_t = sqrt(2.0)*mt/vev only valid for mt(mt); need mt(mZ) here
  // Top Yukawa lambda from Manohar
  //Complex lambda_t = 1.02858; 
  // Top Yukawa lambda from Sascha          
  double lambda_t = 0.991172;  
  // Quartic coupling lambda (need to multiply by a factor of 2 when accessing the quartic coupling)
  double lambda = (mH_/vev)*(mH_/vev);
  y[29] = lambda_t;
  y[30] = lambda;
  y[31] = vev/GeV;
}

void EWCouplings::derivatives(const double x, vector<Complex> & y, 
			      vector<Complex> &dydx) {
  // zero the output
  for (unsigned int i=0; i<dydx.size(); i++) dydx[i]=0.0;
  // low scale
  if (y.size()==2) {
    lowBetaGauge(x,y,dydx);
  }
  // high scale
  else {
    betaGauge(x,y,dydx);
    betaYukawa(x,y,dydx);
    betaHiggs(x,y,dydx);
  }
}

void EWCouplings::betaGauge(const double x, vector<Complex> &y, vector<Complex> & dydx) {
  using Constants::pi;
  using boost::numeric::ublas::axpy_prod;
  using boost::numeric::ublas::herm;
  const Complex I(0,1.0);
  // Yukawa
  boost::numeric::ublas::matrix<Complex> Yuk_u(3,3), Yuk_d(3,3), Yuk_e(3,3);
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Yuk_u(ix,iy) = y[21+3*ix+iy];
      Yuk_d(ix,iy) = y[12+3*ix+iy];
      Yuk_e(ix,iy) = y[ 3+3*ix+iy];
    }
  }
  // gauge
  boost::numeric::ublas::vector<Complex> gauge(3);
  for(unsigned int l=0; l<3; l++) gauge[l] = y[l];
  // Evaluate beta functions for gauge couplings
  double Ng = 0.5*numberOfFlavours(x);
  boost::numeric::ublas::vector<Complex>  b(3), g2(3), gc(3), Cu(3), Cd(3), Ce(3);
  boost::numeric::ublas::matrix<Complex> B1(3,3),B2(3,3), B3(3,3), B(3,3);
  b[0] = -4.0/3.0*Ng - 1.0/10.0;
  b[1] = 22.0/3.0 - 4.0/3.0*Ng - 1.0/6.0;
  b[2] = 11.0 - 4.0/3.0*Ng;
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      B1(ix,iy) = 0.;
      B2(ix,iy) = 0.;
      B3(ix,iy) = 0.;
    }
  }
  B1(1,1) = 136.0/3.0;
  B1(2,2) = 102.0;
  B2(0,0) = 19.0/15.0;
  B2(0,1) = 1.0/5.0;
  B2(0,2) = 11.0/30.0;
  B2(1,0) = 3.0/5.0;
  B2(1,1) = 49.0/3.0;
  B2(1,2) = 3.0/2.0 ;
  B2(2,0) = 44.0/15.0;
  B2(2,1) = 4.0;
  B2(2,2) = 76.0/3.0;
  B3(0,0) = 9.0/50.0;
  B3(0,1) = 3.0/10.0;
  B3(1,0) = 9.0/10.0;
  B3(1,1) = 13.0/6.0;
  B = B1 - Ng*B2 - B3;
  Cu[0] = 17.0/10.0;
  Cu[1] = 3.0/2.0;
  Cu[1] = 2.0;
  Cd[0] = 1.0/2.0;
  Cd[1] = 3.0/2.0;
  Cd[1] = 2.0;
  Ce[0] = 3.0/2.0;
  Ce[1] = 1.0/2.0;
  Ce[1] = 0.0;
  for (int i=0; i<3; i++) {
    g2(i) = pow(gauge(i),2);
  }
  // gc = trans(g2) * B
  axpy_prod(B, g2, gc, true);
  // compute the answer
  if(loops_ >= 1) {
    for(int l=0; l<3; l++) {
      dydx[l] += -b(l)*pow(gauge(l),3)/(16.0*pow(pi,2));
    }
    if (loops_ >= 2) {
      boost::numeric::ublas::matrix<Complex> temp;
      axpy_prod(herm(Yuk_u),Yuk_u,temp);
      Complex tr1 = trace(temp);
      axpy_prod(herm(Yuk_d),Yuk_d,temp);
      Complex tr2 = trace(temp);
      axpy_prod(herm(Yuk_e),Yuk_e,temp);
      Complex tr3 = trace(temp);
      for(int l=0; l<3; l++) {
	dydx[l] += -pow(gauge(l),3)/pow(16.0*pow(pi,2),2)*
	  (gc(l) + Cu(l)*tr1 + Cd(l)*tr2 + Ce(l)*tr3 );
      }
    }
  }
}

void EWCouplings::lowBetaGauge(const double, vector<Complex> &y,
			       vector<Complex> &dydx) {
  using Constants::pi;
  const Complex I(0,1.0);
  Complex e = y[0],  g3 = y[1];
  // Evaluate beta functions for gauge couplings
  double Nu = 2.0, Nd = 3.0, Nl = 3.0;
  if(loops_ >=1) {
    dydx[0] += (16.0/9.0*Nu + 4.0/9.0*Nd + 4.0/3.0*Nl)*pow(e,3)/pow(4.0*pi,2);
    dydx[1] += (2.0/3.0*(Nu+Nd)-11.0)*pow(g3,3)/pow(4.0*pi,2);
    // Note this also includes the three-loop contribution for alpha_3
    if (loops_ >= 2) { 
      dydx[0] += (64.0/27.0*Nu+4.0/27.0*Nd+4.0*Nl)*pow(e,5)/pow(4.0*pi,4) +
	(64.0/9.0*Nu+16.0/9.0*Nd)*pow(e,3)*pow(g3,2)/pow(4.0*pi,4);
      dydx[1] += (38.0/3.0*(Nu+Nd)-102.0)*pow(g3,5)/pow(4.0*pi,4) + 
	(8.0/9.0*Nu+2.0/9.0*Nd)*pow(g3,3)*pow(e,2)/pow(4.0*pi,4) + 
	(5033.0/18.0*(Nu+Nd)-325.0/54.0*(Nu+Nd)*(Nu+Nd)-2857.0/2.0)*pow(g3,7)/pow(4.0*pi,6);
    }
  }
}

void EWCouplings::betaYukawa(const double x, vector< Complex > &y, vector<Complex > &dydx) {
  using Constants::pi;
  const Complex I(0,1.0);
  boost::numeric::ublas::identity_matrix<Complex> Id(3,3);
  // Yukawa
  boost::numeric::ublas::matrix<Complex> Yuk_u(3,3), Yuk_d(3,3), Yuk_e(3,3);
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Yuk_u(ix,iy) = y[21+3*ix+iy];
      Yuk_d(ix,iy) = y[12+3*ix+iy];
      Yuk_e(ix,iy) = y[ 3+3*ix+iy];
    }
  }
  // gauge
  double Ng = 0.5*numberOfFlavours(x);
  boost::numeric::ublas::vector<Complex> gauge(3);
  for(unsigned int l=0; l<3; l++) gauge[l] = y[l];
  Complex lambda = y[30];
  // traces of yukawa matrices
  boost::numeric::ublas::matrix<Complex> mTemp,MUU,MDD,MLL,MUU2,MDD2,MLL2,MUUDD,MDDUU;
  axpy_prod(herm(Yuk_u),Yuk_u,MUU);
  Complex trU  = trace( MUU);
  axpy_prod(MUU,MUU,MUU2);
  Complex trUU = trace(MUU2);
  axpy_prod(herm(Yuk_d),Yuk_d,MDD);
  Complex trD  = trace( MUU);
  axpy_prod(MDD,MDD,MDD2);
  Complex trDD = trace(MDD2);
  axpy_prod(MUU,MDD,MUUDD);
  Complex trUD = trace(MUUDD);
  axpy_prod(MDD,MUU,MDDUU);
  axpy_prod(herm(Yuk_e),Yuk_e,MLL);
  Complex trL  = trace( MLL);
  axpy_prod(MLL,MLL,MLL2);
  Complex trLL = trace(MLL2);
  Complex g02 = sqr(gauge[0]);
  Complex g12 = sqr(gauge[1]);
  Complex g22 = sqr(gauge[2]);
  // Evaluate beta functions for yukawa couplings
  boost::numeric::ublas::zero_matrix<Complex> zero3x3(3,3);
  boost::numeric::ublas::matrix<Complex> dYuk_udx(zero3x3), dYuk_ddx(zero3x3), dYuk_edx(zero3x3), beta1_u(zero3x3),
    beta1_d(zero3x3), beta1_e(zero3x3), beta2_u(zero3x3), beta2_d(zero3x3), beta2_e(zero3x3);
  Complex Y2 = 3.0*trU+3.0*trD + trL;
  Complex Y4 = (17.0/20.0*g02 + 9.0/4.0*g12 + 8.0*g22)*trU +
    (1.0/4.0*g02 + 9.0/4.0*g12 + 8.0*g22)*trD + 3.0/4.0*(g02 + g12)*trL;
  Complex chi4 = 27.0/4.0*trUU + 27.0/4.0*trDD + 9.0/4.0*trLL - 6.0/4.0*trUD;
  if(loops_ >= 1) {
    beta1_u = 3.0/2.0*(MUU - MDD) + (Y2 - 17.0/20.0*g02 - 9.0/4.0*g12 - 8.0*g22)*Id;
    beta1_d = 3.0/2.0*(MDD - MUU) + (Y2 - 1.0/4.0*g02 - 9.0/4.0*g12 - 8.0*g22)*Id;
    beta1_e = 3.0/2.0*(MLL) + (Y2 - 9.0/4.0*g02 - 9.0/4.0*g12)*Id;
    
    axpy_prod(Yuk_u,beta1_u,mTemp);
    dYuk_udx += (1.0/(16.0*pow(pi,2)))*mTemp;
    axpy_prod(Yuk_d,beta1_d,mTemp);
    dYuk_ddx += (1.0/(16.0*pow(pi,2)))*mTemp;
    axpy_prod(Yuk_e,beta1_e,mTemp);
    dYuk_edx += (1.0/(16.0*pow(pi,2)))*mTemp;
      
    if (loops_ >= 2) {
      beta2_u = 3.0/2.0*MUU2 - MUUDD - 1.0/4.0*MDDUU + 11.0/4.0*MDD2 + Y2*(5.0/4.0*MDD - 9.0/4.0*MUU) +
	(-chi4 + 3.0/2.0*pow(lambda,2))*Id - 2.0*lambda*(3.0*MUU + MDD) + (223.0/80.0*g02 + 135.0/16.0*g12 + 16.0*g22)*(MUU) 
	- (43.0/80.0*g02 - 9.0/16.0*g12 + 16.0*g22)*(MDD) + 5.0/2.0*Y4*Id + 
	((9.0/200.0 + 29.0/45.0*Ng)*pow(gauge[0],4) - 9.0/20.0*pow(gauge[0]*gauge[1],2) + 19.0/15.0*pow(gauge[0]*gauge[2],2) - (35.0/4.0 - Ng)*pow(gauge[1],4) + 9.0*pow(gauge[1]*gauge[2],2) - (404.0/3.0 - 80.0/9.0*Ng)*pow(gauge[2],4))*Id;
      beta2_d = 3.0/2.0*MDD2 - MDDUU - 1.0/4.0*MUUDD + 11.0/4.0*MUU2 + Y2*(5.0/4.0*MUU - 9.0/4.0*MDD) + (-chi4 + 3.0/2.0*pow(lambda,2))*Id - 2.0*lambda*(3.0*MDD + MUU) + (187.0/80.0*g02 + 135.0/16.0*g12 + 16.0*g22)*(MDD) - (79.0/80.0*g02 - 9.0/16.0*g12 + 16.0*g22)*(MUU) + 5.0/2.0*Y4*Id - ((29.0/200.0 + 1.0/45.0*Ng)*pow(gauge[0],4) - 27.0/20.0*pow(gauge[0]*gauge[1],2) + 31.0/15.0*pow(gauge[0]*gauge[2],2) - (35.0/4.0 - Ng)*pow(gauge[1],4) + 9.0*pow(gauge[1]*gauge[2],2) - (404.0/3.0 - 80.0/9.0*Ng)*pow(gauge[2],4))*Id;
      beta2_e = 3.0/2.0*MLL2 - 9.0/4.0*Y2*MLL + (-chi4 + 3.0/2.0*pow(lambda,2))*Id - 6.0*lambda*(MLL) + (387.0/80.0*g02 + 135.0/15.0*g12)*(MLL) + 5.0/2.0*Y4*Id + ((51.0/200.0 + 11.0/5.0*Ng)*pow(gauge[0],4) + 27.0/20.0*pow(gauge[0]*gauge[1],2) - (35.0/4.0 - Ng)*pow(gauge[1],4))*Id;
         
      axpy_prod(Yuk_u,beta2_u,mTemp);
      dYuk_udx += (1.0/pow(16.0*pow(pi,2),2))*mTemp;
      axpy_prod(Yuk_d,beta2_d,mTemp);
      dYuk_ddx += (1.0/pow(16.0*pow(pi,2),2))*mTemp;
      axpy_prod(Yuk_e,beta2_e,mTemp);
      dYuk_edx += (1.0/pow(16.0*pow(pi,2),2))*mTemp;
    }
  }
  
  boost::numeric::ublas::vector<Complex> temp(27);

  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      dydx[ 3+ix+3*iy] = dYuk_edx(ix,iy);
      dydx[12+ix+3*iy] = dYuk_ddx(ix,iy);
      dydx[21+ix+3*iy] = dYuk_udx(ix,iy);
    }
  }
}

void EWCouplings::betaHiggs(const double x, vector<Complex> &y,
			    vector<Complex> &dydx) {
  using Constants::pi;
  const Complex I(0,1.0);
  double Ng = 0.5*numberOfFlavours(x);
  // Yukawa
  boost::numeric::ublas::matrix<Complex> Yuk_u(3,3), Yuk_d(3,3), Yuk_e(3,3);
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Yuk_u(ix,iy) = y[21+3*ix+iy];
      Yuk_d(ix,iy) = y[12+3*ix+iy];
      Yuk_e(ix,iy) = y[ 3+3*ix+iy];
    }
  }
  // gauge
  boost::numeric::ublas::vector<Complex> gauge(3);
  for(int l=0; l<3; l++) gauge(l) = y[l];
  Complex lambda = y[30];
  complex<Energy> vev = y[31]*GeV;
  // Evaluate beta functions for higgs coupling
  Complex beta1_lambda(0.), beta2_lambda(0.), gamma1_vev(0.), gamma2_vev(0.),
    Y2(0.), H(0.), Y4(0.), chi4(0.);
  // traces of yukawa matrices
  boost::numeric::ublas::matrix<Complex> temp,temp2,MUU,MDD,MLL;
  axpy_prod(herm(Yuk_u),Yuk_u,MUU);
  Complex trU  = trace( MUU);
  axpy_prod(MUU,MUU,temp);
  Complex trUU = trace(temp);
  axpy_prod(MUU,temp,temp2);
  Complex trUUU = trace(temp2);
  axpy_prod(herm(Yuk_d),Yuk_d,MDD);
  Complex trD  = trace( MUU);
  axpy_prod(MDD,MDD,temp);
  Complex trDD = trace(temp);
  axpy_prod(MDD,temp,temp2);
  Complex trDDD = trace(temp2);
  axpy_prod(MUU,MDD,temp);
  Complex trUD = trace(temp);
  axpy_prod(herm(Yuk_e),Yuk_e,MLL);
  Complex trL  = trace( MLL);
  axpy_prod(MLL,MLL,temp);
  Complex trLL = trace(temp);
  axpy_prod(MLL,temp,temp2);
  Complex trLLL = trace(temp2);
  axpy_prod(MUU+MDD,MDD,temp);
  axpy_prod(MUU,temp,temp2);
  Complex trUUDD = trace(temp2);

  Complex g02 = sqr(gauge[0]);
  Complex g12 = sqr(gauge[1]);
  Complex g22 = sqr(gauge[2]);
  Y2 = 3.0*trU+3.0*trD + trL;
  Y4 = (17.0/20.0*g02 + 9.0/4.0*g12 + 8.0*g22)*(trU) + (1.0/4.0*g02 + 9.0/4.0*g12 + 8.0*g22)*(trD) + 3.0/4.0*(g02 + g12)*(trL);
  chi4 = 27.0/4.0*trUU + 27.0/4.0*trDD + 9.0/4.0*trLL - 6.0/4.0*trUD;
  H = 3.0*trUU + 3.0*trDD + trLL;
   
  if(loops_ >= 1) {
    beta1_lambda = 12.0*pow(lambda,2) - (9.0/5.0*g02 + 9.0*g12)*lambda + 9.0/4.0*(3.0/25.0*pow(gauge[0],4)+2.0/5.0*pow(gauge[0]*gauge[1],2) + pow(gauge[1],4)) + 4.0*Y2*lambda - 4.0*H;
    gamma1_vev = 9.0/4.0*(1.0/5.0*g02+g12)-Y2;
      
    dydx[30] += 1.0/(16.0*pow(pi,2))*beta1_lambda;
    dydx[31] += vev/(16.0*pow(pi,2))*gamma1_vev/GeV;
    
    if (loops_ >= 2) {
      beta2_lambda = -78.0*pow(lambda,3) + 18.0*(3.0/5.0*g02 + 3.0*g12)*pow(lambda,2) - ( (313.0/8.0 - 10.0*Ng)*pow(gauge[1],4) - 117.0/20.0*pow(gauge[0]*gauge[1],2) + 9.0/25.0*(229.0/4.0+50.0/9.0*Ng)*pow(gauge[0],4) )*lambda + (497.0/8.0 - 8.0*Ng)*pow(gauge[1],6) - 3.0/5.0*(97.0/24.0 + 8.0/3.0*Ng)*g02*pow(gauge[1],4) - 9.0/25.0*(239.0/24.0 + 40.0/9.0*Ng)*pow(gauge[0],4)*g12 - 27.0/125.0*(59.0/24.0 + 40.0/9.0*Ng)*pow(gauge[0],6) - 64.0*g22*(trUU + trDD) - 8.0/5.0*g02*(2.0*trUU - trDD + 3.0*trLL) - 3.0/2.0*pow(gauge[1],4)*Y4 + 10.0*lambda*( (17.0/20.0*g02 + 9.0/4.0*g12 + 8.0*g22)*trU + (1.0/4.0*g02 + 9.0/4.0*g12 + 8.0*g22)*trD + 3.0/4.0*(g02 + g12)*trL ) + 3.0/5.0*g02*( (-57.0/10.0*g02 + 21.0*g12 )*trU + (3.0/2.0*g02 + 9.0*g12)*trD + (-15.0/2.0*g02 + 11.0*g12)*trL ) - 24.0*pow(lambda,2)*Y2 - lambda*H + 6.0*lambda*trUD + 20.0*(3.0*trUUU + 3.0*trDDD + trLLL) - 12.0*trUUDD;
      gamma2_vev = -3.0/2.0*pow(lambda,2) - 5.0/2.0*Y4 + chi4 - 27.0/80.0*pow(gauge[0]*gauge[1],2) - (93.0/800.0 + 1.0/2.0*Ng)*pow(gauge[0],4) + (511.0/32.0 - 5.0/2.0*Ng)*pow(gauge[1],4);
         
      dydx[30] += 1.0/pow(16.0*pow(pi,2),2)*beta2_lambda;
      dydx[31] += vev/pow(16.0*pow(pi,2),2)*gamma2_vev/GeV;
    }
  }
}

double EWCouplings::interpolate(double t, int paramIndex) {
  double stepsize = table_(1,0)-table_(0,0);
  double tol = 0.001*stepsize;
	
   double logewScale   = log(ewScale_/GeV);
   double loghighScale = log(highScale_/GeV);
   
   if (t<logewScale-tol || t>loghighScale+tol) {
     cerr << "Stepsize: " << stepsize << std::endl;
     cerr << "tol: " << tol << std::endl;
     cerr << "logewScale: " << logewScale << std::endl;
     cerr << "loghighScale: " << loghighScale << std::endl;
     cerr << "paramIndex: " << paramIndex << std::endl;
     cerr << "t: " << t << std::endl;
     
     cerr << "Couplings::_Interp(double t, int parmIndex) trying to obtain parameter ";
     cerr << "value outside the available range. Returning zero." << std::endl;
     assert(false);
   }
   
   // return value at EW scale
   if (abs(t-logewScale)<tol) return table_(0,paramIndex);
   
   // return value at high scale
   if (abs(t-loghighScale)<tol) {
     return table_(table_.size1()-1,paramIndex);
   }
	
   unsigned int numSteps = int((t-table_(0,0))/stepsize);
   
   // Linear Interpolation:
   double x1 = table_(numSteps,0);
   double y1 = table_(numSteps,paramIndex);
   double x2 = table_(numSteps+1,0);
   double y2 = table_(numSteps+1,paramIndex);
   return y1+((y2-y1)/(x2-x1))*(t-x1);
}

double EWCouplings::interpolateLow(double t, int paramIndex) {
  double stepsize = lowTable_(0,0)-lowTable_(1,0);
  double tol = 0.00001*stepsize;
	
  double logewScale  = log(ewScale_ /GeV);
  double loglowScale = log(lowScale_/GeV);
   
  if (t<loglowScale-tol || t>logewScale+tol) {
    cerr<< "Couplings::_LowInterp(double t, int parmIndex) trying to obtain parameter ";
    cerr << "value outside the available range. Returning zero." << std::endl;
    assert(false);
  }
   
  if (abs(t-logewScale)<tol) {
    return lowTable_(0,paramIndex);
  }
   
  if (std::abs(t-loglowScale)<tol) {
    return lowTable_(lowTable_.size1()-1,paramIndex);
  }
	
  int numSteps = (int)((lowTable_(0,0)-t)/stepsize);
   
  // Linear Interpolation:
  double x1 = lowTable_(numSteps,0);
  double y1 = lowTable_(numSteps,paramIndex);
  double x2 = lowTable_(numSteps+1,0);
  double y2 = lowTable_(numSteps+1,paramIndex);
  
  return y1+((y2-y1)/(x2-x1))*(t-x1);
}














