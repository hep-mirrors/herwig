// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopPOWHEGDecayer class.
//

#include "SMTopPOWHEGDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMTopPOWHEGDecayer::SMTopPOWHEGDecayer() {}

IBPtr SMTopPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMTopPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}

void SMTopPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void SMTopPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SMTopPOWHEGDecayer,SMTopDecayer>
describeHerwigSMTopPOWHEGDecayer("Herwig::SMTopPOWHEGDecayer", "HwPerturbativeDecay.so");

void SMTopPOWHEGDecayer::Init() {

  static ClassDocumentation<SMTopPOWHEGDecayer> documentation
    ("There is no documentation for the SMTopPOWHEGDecayer class");

}

HardTreePtr SMTopPOWHEGDecayer::generateHardest(ShowerTreePtr) {
  using Constants::pi;
  using Constants::twopi;
  // test of the phase-space
  Energy m0 = getParticleData(ParticleID::t    )->mass();
  Energy m1 = getParticleData(ParticleID::b    )->mass();
  Energy m2 = getParticleData(ParticleID::Wplus)->mass();
  Energy m3 = ZERO;
  Energy2 m0s(sqr(m0)),m1s(sqr(m1)),m2s(sqr(m2)),m3s(ZERO);
  double mu1  = m1/m0   , mu2  = m2/m0   ;
  double mu12 = sqr(mu1), mu22 = sqr(mu2);
  double lambda = sqrt(1.+sqr(mu12)+sqr(mu22)-2.*mu12-2.*mu22-2.*mu12*mu22);
  unsigned long npoint(1000000);
  Energy2 wgtsum(ZERO);
  Energy4 wgtsq (ZERO);
  ofstream file("dalitz.top");
  file << "SET FONT DUPLEX\n";
  file << "SET ORDER X Y\n";
  file << "SET LIMITS X 0 0.9 Y 0.9 1.3\n";
  unsigned int nplot(0);
  Energy2 wgtmax = sqr(m0)/sqr(pi)/lambda/16.;
  for(unsigned long ix=0;ix<npoint;++ix) {
    Energy2 m122max = sqr(m0);
    Energy2 m122min = sqr(m1+m2);
    Energy2 m122 = m122min+UseRandom::rnd()*(m122max-m122min);
    Energy m12 = sqrt(m122);
    Energy E2 = 0.5*(m122-m1s+m2s)/m12;
    Energy E3 = 0.5*(m0s-m122-m3s)/m12;
    Energy p2 = sqrt(sqr(E2)-m2s);
    Energy p3 = sqrt(sqr(E3)-m3s);
    Energy2 m232max = sqr(E2+E3)-sqr(p2-p3);
    Energy2 m232min = sqr(E2+E3)-sqr(p2+p3);
    Energy2 m232 =  m232min+UseRandom::rnd()*(m232max-m232min);
    Energy2 wgt = (m122max-m122min)*(m232max-m232min)/sqr(pi)/sqr(m0)/lambda/16.;

    double x1 = 1.-(m232-m1s)/m0s;
    double xg = 1.-(m122)/m0s;
    double x2=2.-x1-xg;

    if(nplot<100000) {
      if(wgt>wgtmax)
	cerr << "weight problem " << wgt/wgtmax << "\n";
      if(wgt>UseRandom::rnd()*wgtmax) {
	++nplot;
	file << xg << "\t" << x2 << "\n";
	if(nplot%50000==0) file << "PLOT\n";
      }
    }

    wgtsum += wgt;
    wgtsq  += sqr(wgt);
  }
  file << "PLOT\nNEW FRAME\n";
  file << "SET ORDER X Y\n";
  file << "SET LIMITS X 0 0.9 Y 0.9 1.3\n";
  Energy2 vol = wgtsum/double(npoint);
  Energy2 err = sqrt((wgtsq /double(npoint)-sqr(vol))/double(npoint));

  Energy2 avol = 1./sqr(2.*pi)/4./lambda
    *m0s * (- mu22  * (1. - mu12) * log((1. - mu12 + mu22 + lambda) / 
					(1. - mu12 + mu22 - lambda)) 
	    - mu12  * (1. - mu22) * log((1. + mu12 - mu22 + lambda) / 
					(1. + mu12 - mu22 - lambda)) 
	    + lambda * (1. + mu12 + mu22) / 2.);
  cerr << "testing analytic " << vol/GeV2 << "\n";
  cerr << "phase-space volume is " << vol/GeV2 << " +/- " << err/GeV2 << "\n";
  cerr << "ratio              is " << vol/avol << " +/- " << err/avol << "\n";



  Energy pTmax = 0.5*(m0-m2-sqr(m1)/(m0-m2));
  double ymin=-10.,ymax=10.;
  wgtsum = ZERO;
  wgtsq  = ZERO;
  nplot=0;
  Energy maximumPt(ZERO);
  wgtmax *=500.;
  for(unsigned long ix=0;ix<npoint;++ix) {
    Energy pT = UseRandom::rnd()*pTmax;
    double phi = twopi*UseRandom::rnd();
    double y = ymin+UseRandom::rnd()*(ymax-ymin);
    double xT = 2.*pT/m0, xT2 = sqr(xT);
    double x3 = 2.*pT*cosh(y)/m0;
    // solve for x1
    double disc = 16.*(sqr(x3)-xT2)*((1.-x3-sqr(mu1+mu2))*(1.-x3-sqr(mu1-mu2))-sqr(mu2)*xT2);
    double a =  4.-4.*x3+xT2;
    double b = -2.*(2.-x3)*(2.*(1.-x3+sqr(mu1)-sqr(mu2))+xT2);
    if(disc<0.) continue;
    disc = sqrt(disc);
    double x1[2] = {0.5*(-b+disc)/a,0.5*(-b-disc)/a};
    double x2[2];
    Energy2 wgt(ZERO);
    for(unsigned int iy=0;iy<2;++iy) {
      // check x1 within limits
      if(x1[iy]<2.*mu1 || x1[iy]>1.+mu12-mu22 ) continue;
      // calc x2 and check limits
      x2[iy] = 2-x3-x1[iy];
      double root = sqrt(max(0.,sqr(x1[iy])-4.*mu12));



      double x2Plus  = 1.-mu12+mu22-0.5*(1.-x1[iy]+mu12-mu22)/(1.-x1[iy]+mu12)
	*(x1[iy]-2.*mu12-root);
      double x2Minus = 1.-mu12+mu22-0.5*(1.-x1[iy]+mu12-mu22)/(1.-x1[iy]+mu12)
	*(x1[iy]-2.*mu12+root);


      if(x2[iy]>x2Plus || x2[iy]<x2Minus) continue;
      static double tolerance = 1e-6; 
      bool isMomentaReconstructed = false;
      Lorentz5Momentum p3(pT*cos(phi),pT*sin(phi),pT*sinh(y),pT*cosh(y),ZERO);
      Lorentz5Momentum p1,p2;
      double z1 = sqrt(sqr(x1[iy])-xT2-4.*mu12);
      double z2 = sqr(x2[iy])-4.*mu22;
      if(z2>0.) z2 = sqrt(z2);
      else {
	cerr << "problem with z2 " << z2 << "\n";
	z2=0.;
      }
      double z3 = sqrt(sqr(x3)-xT2);
      if(pT*sinh(y)>ZERO) {
	if( abs(-z2+z3+z1) <= tolerance) {
	  p1 = Lorentz5Momentum(-pT*cos(phi),-pT*sin(phi),
				0.5*m0*z1,0.5*m0*x1[iy]);
	  p2 = Lorentz5Momentum(ZERO,ZERO,-0.5*m0*z2,
				0.5*m0*x2[iy],m2);
	  isMomentaReconstructed=true;
	}
	else if( abs(-z2+z3-z1) <= tolerance) {
	  p1 = Lorentz5Momentum(-pT*cos(phi),-pT*sin(phi),
				-0.5*m0*z1,0.5*m0*x1[iy]);
	  p2 = Lorentz5Momentum(ZERO,ZERO,-0.5*m0*z2,
				0.5*m0*x2[iy],m2);
	  isMomentaReconstructed=true;
	}
      }
      else if(pT*sinh(y) < ZERO){
	if( abs(-z2-z3+z1) <= tolerance) {
	  p1 = Lorentz5Momentum(-pT*cos(phi),-pT*sin(phi),
				0.5*m0*z1,0.5*m0*x1[iy]);
	  p2 = Lorentz5Momentum(ZERO,ZERO,-0.5*m0*z2,
				0.5*m0*x2[iy],m2);
	  isMomentaReconstructed=true;
	}
	else if ( abs(-z2-z3-z1) <= tolerance) {
	  p1 = Lorentz5Momentum(-pT*cos(phi),-pT*sin(phi),
				0.5*m0*z1,0.5*m0*x1[iy]);
	  p2 = Lorentz5Momentum(ZERO,ZERO,-0.5*m0*z2,
				0.5*m0*x2[iy],m2);
	  isMomentaReconstructed=true;
	}
      }
      else if(abs(-z2+z1) <= tolerance) {
	p1 = Lorentz5Momentum(-pT*cos(phi),-pT*sin(phi),
			      0.5*m0*z1,0.5*m0*x1[iy]);
	p2 = Lorentz5Momentum(ZERO,ZERO,-0.5*m0*z2,
			      0.5*m0*x2[iy],m2);
	isMomentaReconstructed=true;
      }
      if(!isMomentaReconstructed) continue;


//       cerr << "testing masses " << m1/GeV << " " << m2/GeV << "\n";


//       cerr << "testing p1 z " << p1.m2()/GeV2 << " "
// 	   << (-sqr(pT)-0.25*sqr(m0)*(sqr(x1[iy])-xT2-4.*mu12)
// 	       +0.25*sqr(m0)*sqr(x1[iy]))/GeV2<< " "
// 	   << (-sqr(pT)-0.25*sqr(m0)*(-xT2-4.*mu12))/GeV2
// 	   << "\n";
//       cerr << "testing pieces "
// 	   << (-sqr(pT)+0.25*sqr(m0)*xT2)/GeV2
// 	   << " "
// 	   << -0.25*sqr(m0)*(-4.*mu12)/GeV2
// 	   << "\n";

//       cerr << "testing bottom " << p1/GeV << " " << p1.m()/GeV << " " 
// 	   << p1.mass()/GeV << "\n";
//       cerr << "testing W bos  " << p2/GeV << " " << p2.m()/GeV << " " 
// 	   << p2.mass()/GeV << "\n";
//       cerr << "testing gluon  " << p3/GeV << " " << p3.m()/GeV << " " 
// 	   << p3.mass()/GeV << "\n";

//       cerr << "testing x sum " << x1[iy]+x2[iy]+x3 << "\n";

//       cerr << "testing p sum " << (p1+p2+p3)/GeV << "\n";


      Energy ew=p2.e();
      Energy pw=abs(p2.z());

      wgt += abs(pT*pTmax*(ymax-ymin)/sqr(twopi)/lambda*sqr(pw)/(pw*(m0-pT*cosh(y))-ew*pT*sinh(y)));

//       wgt += pT*pTmax*(ymax-ymin)/sqr(twopi)/lambda/2.*pow(sqr(x2[iy])-4*mu22,1.5)/
// 	fabs(2.*x1[iy]*mu22+mu22*x2[iy]+x2[iy]+mu12*x2[iy]-x2[iy]*x1[iy]-sqr(x2[iy]));

//       cerr << "testing weights "
// 	   <<  pT*pTmax*(ymax-ymin)/sqr(twopi)/lambda*sqr(pw)/(pw*(m0-pT*cosh(y))-ew*pT*sinh(y))/GeV2 << " "
// 	   <<pT*pTmax*(ymax-ymin)/sqr(twopi)/lambda/2.*pow(sqr(x2[iy])-4*mu22,1.5)/
// 	fabs(2.*x1[iy]*mu22+mu22*x2[iy]+x2[iy]+mu12*x2[iy]-x2[iy]*x1[iy]-sqr(x2[iy]))/GeV2
// 	   << "\n";

      
      if(nplot<100000) {
	if(wgt>wgtmax)
	  cerr << "weight problem " << wgt/wgtmax << " "
	       << pw*(m0-pT*cosh(y))/GeV2 << " " << ew*pT*sinh(y)/GeV2 << " "
	       << pT/GeV << " " << pw/GeV << " " << y << "\n";
	if(wgt>UseRandom::rnd()*wgtmax) {
	  ++nplot;
	  file << x3 << "\t" << x2[iy] << "\n";
	  if(nplot%50000==0) file << "PLOT\n";
	}
      }
      if(pT>maximumPt) maximumPt = pT;
    }
    wgtsum += wgt;
    wgtsq += sqr(wgt);
  }
  file << "PLOT\n";
  double xw=2.*mu2;
  while(xw<=1.+mu22-mu12) {
    double xg = 1.+(xw*(mu22-mu12)-2.*mu22)/(xw-2.*mu22);
    file << xg << "\t" << xw << "\n";
    xw+=0.001;
  }
  file << "JOIN RED\n";
  xw=2.*mu2;
  while(xw<=1.+mu22-mu12) {
    double xg =(2.*(-1.+xw+mu12-mu22))/(-2.+xw);
    file << xg << "\t" << xw << "\n";
    xw+=0.001;
  }
  file << "JOIN BLUE\n";
  xw=2.*mu2;
  while(xw<=1.+mu22-mu12) {
    double xg = 1-mu12+mu22-2*mu22*(1+mu12-mu22)/(xw-2*mu22);
    file << xg << "\t" << xw << "\n";
    xw+=0.001;
  }
  file << "JOIN GREEN\n";




  file.close();
  vol = wgtsum/double(npoint);
  err = sqrt((wgtsq /double(npoint)-sqr(vol))/double(npoint));

  cerr << "phase-space volume B is " << vol/GeV2 << " +/- " << err/GeV2 << "\n";
  cerr << "ratio              B is " << vol/avol << " +/- " << err/avol << "\n";

  cerr << "testing max pT " << maximumPt/GeV << " "
       << 0.5*(m0-m2-sqr(m1)/(m0-m2))/GeV << "\n";

  exit(0);
  return HardTreePtr();
}
