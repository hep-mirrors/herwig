#include "QQbarG.h"
#include <math.h>
#include <stdlib.h>

#define sqr(x) ((x)*(x))
#define EPS 0.00000001

using namespace Herwig;

QQbarG::QQbarG(Energy Q, Energy m) {
  d_Q = Q;
  d_m = m;
  setRho(sqr(d_m/d_Q));
  setKtildeSymm();
  std::srand(1); 
}

QQbarG::QQbarG(Energy Q, Energy m, double k) {
  d_Q = Q;
  d_m = m;
  setRho(sqr(d_m/d_Q));
  setKtilde(k);
}

QQbarG::~QQbarG() {}

void QQbarG::setRho(double r) { 
  d_rho = r; d_v = sqrt(1.-4.*d_rho);
}

void QQbarG::setKtilde(double k) { 
  d_kt1 = k; 
  setKtilde2(); 
}

void QQbarG::setKtilde2() { 
   double num = d_rho * d_kt1 + 0.25 * d_v *(1.+d_v)*(1.+d_v);
   double den = d_kt1 - d_rho;
   d_kt2 = num/den;
}

void QQbarG::setKtildeSymm() { 
  d_kt1 = (1. + sqrt(1. - 4.*d_rho))/2.;
  setKtilde2();
}

void QQbarG::setKtildeLargest() { 
  d_kt1 = 4.*(1. - 2.*sqrt(d_rho));
  setKtilde2();
}

void QQbarG::setKtildeSmallest() { 
  setKtildeLargest();
  swap(d_kt1, d_kt2);
}

void QQbarG::getXXbar(double kti, double z, double &x, double &xbar) {
  x = (1. + sqr(d_v)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
       + z*sqrt(sqr(d_v) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z))
       - kti*(-1. + z)*z*(2. + z*(-2 
       + sqrt(sqr(d_v)+ kti*(-1. + z)*z*(2. + kti*(-1. + z)*z))
				  )))
    /(1. - kti*(-1. + z)*z 
      + sqrt(sqr(d_v) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z)));
  xbar = 1. + kti*(-1. + z)*z;
}
 
double QQbarG::u(double x2) {
  return 0.5*(1. + d_rho/(1.-x2+d_rho));
}

double QQbarG::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho);
  return uval + num/den;
}

double QQbarG::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double QQbarG::PS(double x, double xbar) {
   double u = 0.5*(1. + d_rho / (1.-xbar+d_rho));
   double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho);
   double brack = (1.+z*z)/(1.-z)- 2.*d_rho/(1-xbar);
   // interesting: the splitting function without the subtraction
   // term. Actually gives a much worse approximation in the collinear
   // limit.  double brack = (1.+z*z)/(1.-z);
   double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho);
   return brack/den;
}

double QQbarG::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho)*(x1+2.*d_rho) + (x2+2.*d_rho)*(x2+2.*d_rho) 
    - 8.*d_rho*(1.+2.*d_rho);
  double den = (1.+2.*d_rho)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho/((1.-x1)*(1.-x1)) 
	  - 2*d_rho/((1.-x2)*(1.-x2)))/d_v;
}

double QQbarG::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho)*(x1+2.*d_rho) + (x2+2.*d_rho)*(x2+2.*d_rho) 
    + 2.*d_rho*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho);
  double den = d_v*d_v*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho/((1.-x1)*(1.-x1)) 
	  - 2*d_rho/((1.-x2)*(1.-x2)))/d_v;
}

double QQbarG::getRatioV(double x, double xbar) {
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if (k1 < d_kt1 && k2 < d_kt2) return 0.0;
  if (k1 < d_kt1) return MEV(x, xbar)/PS(x, xbar);
  // No...Is it in the anti-quark emission zone?
  if(k2 < d_kt2) return MEV(x, xbar)/PS(xbar, x);
  return 0.0;
}

double QQbarG::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) {
    //    std::cout << "oops" << std::endl; 
    return 0.0;
  }
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
//   std::cout << "(x, xbar) = " << x << ", " << xbar << ")" << std::endl;
//   std::cout << "qWeight: (k1, k2, kt1m, kt2m) = (" 
// 	    << k1 << ", " 
// 	    << k2 << ", " 
// 	    << d_kt1 << ", " 
// 	    << d_kt2 << ")" << std::endl;
  // Is it in the quark emission zone?
  if(k1 < d_kt1) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2) rval *= 0.5;
    return rval;
  }
  //  std::cout << "passed through! " << std::endl; 
  return 1.0;
}

double QQbarG::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) {
    //std::cout << "oops" << std::endl; 
    return 0.0;
  }
  //   double k1 = getKfromX(x, xbar);
  //   double k2 = getKfromX(xbar, x);
  double k1 = getKfromX(xbar, x);
  double k2 = getKfromX(x, xbar);
//   std::cout << "qbarWeight: (k1, k2, kt1m, kt2m) = (" 
// 	    << k1 << ", " 
// 	    << k2 << ", " 
// 	    << d_kt1 << ", " 
// 	    << d_kt2 << ")" << std::endl;
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the anti-quark emission zone?
    if(k1 < d_kt1) rval *= 0.5;
    return rval;
  }
  //  std::cout << "passed through!" << std::endl; 
  return 1.0;
}

double QQbarG::qWeightX(Energy qtilde, double z) {
  double x, xb; 
//   std::cout << "qWeightX called with (" 
// 	    << sqr(qtilde/d_Q) << ", " << z 
// 	    << ")";
  getXXbar(sqr(qtilde/d_Q), z, x, xb);
  //std::cout << " -> (" << x << ", " << xb << ")" << std::endl;
  return qWeight(x, xb); 
}

double QQbarG::qbarWeightX(Energy qtilde, double z) {
  double x, xb; 
//   std::cout << "qbarWeightX called with (" 
// 	    << sqr(qtilde/d_Q) << ", " << z
// 	    << ")";
  getXXbar(sqr(qtilde/d_Q), z, x, xb);
  //  std::cout << " -> (" << x << ", " << xb << ")" << std::endl;
  return qbarWeight(x, xb); 
}

double QQbarG::getHard(double &x1, double &x2) {
  //  x1 = drand48(); 
  //  x2 = drand48(); 
  double w = 0.0;
  // map soft corner to origin:  
  // double y1 = 1.-x1; 
  // double y2 = 1.-x2; 
  // this is equivalent to previous...
  double y1 = drand48(); 
  double y2 = drand48(); 
  // simply double MC efficiency 
  // -> weight has to be divided by two (Jacobian)
  if (y1 + y2 > 1) {
    y1 = 1.-y1; 
    y2 = 1.-y2;
  }
  bool inSoft = false; 
  if (y1 < 0.25) { 
    if (y2 < 0.25) {
      inSoft = true; 
      if (y1 < y2) {
	y1 = 0.25-y1;
	y2 = y1*(1.5 - 2.*y2);
      }	else {
	y2 = 0.25 - y2;
	y1 = y2*(1.5 - 2.*y1);
      }
    } else {
      if (y2 < y1 + 2.*sqr(y1)) return w;
    }
  } else {
    if (y2 < 0.25) {
      if (y1 < y2 + 2.*sqr(y2)) return w;
    }
  } 

  // inside PS?
  x1 = 1.-y1;
  x2 = 1.-y2;
  if(y1*y2*(1.-y1-y2) < d_rho*sqr(y1+y2)) return w;
  double k1 = getKfromX(x1, x2);
  double k2 = getKfromX(x2, x1);
  // Is it in the quark emission zone?
  if (k1 < d_kt1) return 0.0;
  // No...is it in the anti-quark emission zone?
  if (k2 < d_kt2) return 0.0;  
  // Point is in dead zone: compute q qbar g weight
  w = MEV(x1, x2); 
  // for axial: 
  //  w = MEA(x1, x2); 
  // Reweight soft region
  if (inSoft) { 
    if (y1 < y2) w *= 2.*y1;
    else w *= 2.*y2;
  }
  // w *= 2./3./pi*(proper alphas value...)
  w *= 1./3./3.14159265*0.117997; 
  // = CF as/(2 pi) div by two due to 'Jacobian'
  return w; 
}


vector<Lorentz5Momentum> QQbarG::applyHard(const PVector &p) {
  double x, xbar; 
  vector<Lorentz5Momentum> fs; 
  if (getHard(x, xbar) < drand48() || p.size() != 2) {
    return fs; 
  } else {
    Lorentz5Momentum pcm = p[0]->momentum() + p[1]->momentum(); 
    Lorentz5Momentum pq, pa, pg;
    if (abs(p[0]->id()) < 7) {
      if (p[0]->id() > 0) {
	pq = p[0]->momentum(); 
	pa = p[1]->momentum(); 
      } else {
	pa = p[0]->momentum(); 
	pq = p[1]->momentum(); 
      }
    }
    Vector3 beta = (pcm.findBoostToCM()); 
    pq.boost(beta);    
    pa.boost(beta);
//     cout << "applyHard, (x, xbar) = (" 
// 	 << x << ", " << xbar << ")" << endl
// 	 << "  pq   = " << pq
// 	 << ", uq = " << pq.vect().unit() << endl
// 	 << "  pa   = " << pa 
// 	 << ", ua = " << pa.vect().unit() << endl 
// 	 << "  pcm  = " << pcm << endl 
// 	 << "  pcm' = " << pq+pa << endl
// 	 << "  beta = " << beta << endl; 
    Vector3 u1, u2, u3;
    double xg = 2.-x-xbar; 
    if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return fs;

    // moduli of momenta in units of Q and cos theta
    // stick to q direction?
    // p1 is the one that is kept, p2 is the other fermion, p3 the gluon.
    Energy e1, e2, e3; 
    double ct2, ct3, pp1, pp2, pp3;
    bool keepq = true; 
    if (drand48() > sqr(x)/(sqr(x)+sqr(xbar))) 
      keepq = false; 
    if (keepq) {
      pp1 = d_Q*sqrt(sqr(x)-4.*d_rho)/2.;
      pp2 = d_Q*sqrt(sqr(xbar)-4.*d_rho)/2.;
      e1 = d_Q*x/2.; 
      e2 = d_Q*xbar/2.; 
      u1 = pq.vect().unit();
    } else {
      pp2 = d_Q*sqrt(sqr(x)-4.*d_rho)/2.;
      pp1 = d_Q*sqrt(sqr(xbar)-4.*d_rho)/2.;
      e2 = d_Q*x/2.; 
      e1 = d_Q*xbar/2.; 
      u1 = pa.vect().unit();
    }
    pp3 = d_Q*xg/2.;       
    e3 = pp3; 
    u2 = u1.orthogonal();
    u2 /= u2.mag();
    u3 = u1.cross(u2);
    u3 /= u3.mag();
    if (pp1 == 0 || pp2 == 0 || pp3 == 0) {
      if (pp1 == 0) {
	ct2 = 1; 
	ct3 = -1; 
      } 
      if (pp2 == 0 || pp3 == 0) {
	ct2 = 1; 
	ct3 = 1; 
      }
    } else {
      ct3 = (sqr(pp1)+sqr(pp3)-sqr(pp2))/(2.*pp1*pp3);
      ct2 = (sqr(pp1)+sqr(pp2)-sqr(pp3))/(2.*pp1*pp2);
    }
//     cout << "  p1 = " << pp1 << ", e1 = " << e1 
// 	 << ", e12-p12 = " << sqrt(sqr(e1)-sqr(pp1)) << endl
// 	 << "  p2 = " << pp2 << ", e2 = " << e2 
// 	 << ", e22-p22 = " << sqrt(sqr(e2)-sqr(pp2)) << endl
// 	 << "  p3 = " << pp3 << ", e3 = " << e3 
// 	 << ", e32-p32 = " << sqrt(sqr(e3)-sqr(pp3)) << endl; 
    double phi = 2.*pi*drand48();
    double cphi = cos(phi);
    double sphi = sin(phi); 
    double st2 = sqrt(1.-sqr(ct2));
    double st3 = sqrt(1.-sqr(ct3));
    Vector3 pv1, pv2, pv3; 
    pv1 = pp1*u1;
//     pv2 = pp2*(-ct2*u1 + st2*(cphi*u2 + sphi*u3));
//     pv3 = pp3*(-ct3*u1 - st3*(cphi*u2 + sphi*u3));    
    pv2 = -ct2*pp2*u1 + st2*cphi*pp2*u2 + st2*sphi*pp2*u3;
    pv3 = -ct3*pp3*u1 - st3*cphi*pp3*u2 - st3*sphi*pp3*u3;
//     cout << "  pv1 = " << pv1 << ", mag = " << pv1.mag() 
// 	 << ", e1 = " << e1 << endl
// 	 << "  pv2 = " << pv2 << ", mag = " << pv2.mag() 
// 	 << ", e2 = " << e2 << endl
// 	 << "  pv3 = " << pv3 << ", mag = " << pv3.mag() 
// 	 << ", e3 = " << e3 << endl;
    if (keepq) {
      pq = Lorentz5Momentum(pv1, e1);
      pa = Lorentz5Momentum(pv2, e2);
    } else {
      pa = Lorentz5Momentum(pv1, e1);
      pq = Lorentz5Momentum(pv2, e2);
    }
    pg = Lorentz5Momentum(pv3, e3);
//     cout << "  getting FS with " << endl
// 	 << "  (ct2, ct3, phi) = (" 
// 	 << ct2 << ", " << ct3 << ", " << phi << ")" << endl
// 	 << "  keepq? " << (keepq ? "y":"n") << endl
// 	 << "  u1 = " << u1 << ", u12 = " << u1.mag() << endl
// 	 << "  u2 = " << u2 << ", u22 = " << u2.mag() << endl
// 	 << "  u3 = " << u3 << ", u32 = " << u3.mag() << endl
// 	 << "  (u1.u2, u1.u3, u2.u3) = ("
// 	 << u1*u2 << ", " << u1*u3 << ", " << u2*u3 << ")" << endl
// 	 << "  pq   = " << pq << ", m = " << pq.m() << endl
// 	 << "    uq = " << pq.vect().unit() << endl
// 	 << "  pa   = " << pa << ", m = " << pa.m() << endl 
// 	 << "    ua = " << pa.vect().unit() << endl 
// 	 << "  pg   = " << pg << ", m = " << pg.m() << endl 
// 	 << "  pcm  = " << pq+pa+pg << endl;
    pq.boost(-beta);
    pa.boost(-beta);
    pg.boost(-beta);
    fs.push_back(pq); 
    fs.push_back(pa); 
    fs.push_back(pg);
    //    cout << "  pcm' = " << pq+pa+pg << endl; 
  }
  return fs;   
}



#undef EPS
