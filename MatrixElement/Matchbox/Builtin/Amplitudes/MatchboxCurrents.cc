// -*- C++ -*-
//
// MatchboxCurrents.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxCurrents.h"
#include "Herwig++/Utilities/Maths.h"

using namespace Herwig;

using Constants::pi;

inline Complex csqr(Complex a) {
  return a*a;
}

inline double theta(double x) {
  if ( x >= 0. )
    return 1.;
  return 0.;
}

inline double sign(double x) {
  if ( x >= 0. )
    return 1.;
  return -1.;
}

// quick'n'dirty fix to template troubles

Complex operator*(const Complex& a, double b) {
  return Complex(a.real()*b,a.imag()*b);
}

Complex operator*(double b, const Complex& a) {
  return Complex(a.real()*b,a.imag()*b);
}

Complex operator+(const Complex& a, double b) {
  return Complex(a.real()+b,a.imag());
}

Complex operator+(double b, const Complex& a) {
  return Complex(a.real()+b,a.imag());
}

Complex operator-(const Complex& a, double b) {
  return Complex(a.real()-b,a.imag());
}

Complex operator-(double b, const Complex& a) {
  return Complex(b-a.real(),-a.imag());
}


// end fix, needs to be looked at in ThePEG/Config/

void MatchboxCurrents::setupLeptons(int l,    const Lorentz5Momentum& pl,
				    int lbar, const Lorentz5Momentum& plbar) {
  Lorentz5Momentum plFlat = pl;
  plFlat.setMass(ZERO); plFlat.rescaleEnergy();
  if ( pl.t() < ZERO )
    plFlat.setT(-plFlat.t());
  momentum(l,plFlat,true,pl.mass());
  Lorentz5Momentum plbarFlat = plbar;
  plbarFlat.setMass(ZERO); plbarFlat.rescaleEnergy();
  if ( plbar.t() < ZERO )
    plbarFlat.setT(-plbarFlat.t());
  momentum(lbar,plbarFlat,true,plbar.mass());
}

const LorentzVector<Complex>& MatchboxCurrents::llbarLeftCurrent(int l,    int lHel,
								 int lbar, int lbarHel) {

  if ( getCurrent(hash<0>(1,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent(Complex(0.,1.)*minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,2.)*mass(lbar)/plusProduct(l,lbar))*momentum(l));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,-2.)*mass(l)/minusProduct(l,lbar))*momentum(lbar));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,1.)*mass(l)*mass(lbar)/invariant(l,lbar))*minusCurrent(lbar,l));
  }
  return cachedCurrent();
    
}

const LorentzVector<Complex>& MatchboxCurrents::llbarRightCurrent(int l,    int lHel,
								  int lbar, int lbarHel) {
    
  if ( getCurrent(hash<0>(2,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,1.)*mass(l)*mass(lbar)/invariant(l,lbar))*minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,-2.)*mass(l)/plusProduct(l,lbar))*momentum(lbar));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,2.)*mass(lbar)/minusProduct(l,lbar))*momentum(l));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent(Complex(0.,1.)*minusCurrent(lbar,l));
  }
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftCurrent(int q,    int qHel,
								 int qbar, int qbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<1>(1,1,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(Complex(0.,1.)*minusCurrent(q,qbar));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  return cachedCurrent();
    
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightCurrent(int q,    int qHel,
								  int qbar, int qbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<1>(2,1,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(Complex(0.,1.)*minusCurrent(qbar,q));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftCurrent(int q,    int qHel,
								  int qbar, int qbarHel,
								  int g,    int gHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*minusCurrent(q,qbar)
		    +(1./minusProduct(g,qbar))*minusCurrent(q,g)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*minusCurrent(q,qbar)
		    +(1./plusProduct(q,g))*minusCurrent(g,qbar)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightCurrent(int q,    int qHel,
								   int qbar, int qbarHel,
								   int g,    int gHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*minusCurrent(qbar,q)
		    +(1./minusProduct(q,g))*minusCurrent(qbar,g)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*minusCurrent(qbar,q)
		    +(1./plusProduct(g,qbar))*minusCurrent(g,q)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralLeftCurrent(int i, int,
								   int j, int,
								   int k, int g1Hel,
								   int l, int g2Hel,
								   int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,k))/
      (jl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,k))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,l))/
      (kl*(jk + jl + kl)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(i,n))/
      (ik*jl*minusProduct(k,n)) + 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(k,j)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(l,j)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,k)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,l)*minusCurrent(i,j)*sqr(minusProduct(i,n)))/
      (ik*(ik + il + kl)*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(k,j)*sqr(minusProduct(i,n)))/
      (ik*(ik + il + kl)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,n)*minusProduct(j,n))/
      (ik*jl*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,j)*sqr(minusProduct(j,n)))/
      (jl*(jk + jl + kl)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(l,n))/
      (jl*(jk + jl + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*minusProduct(k,n));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*plusProduct(j,k)*plusProduct(j,n)*minusCurrent(i,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(i,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(k,j)*minusProduct(i,l)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(l,j)*minusProduct(i,l)*minusProduct(i,n))/
      (ik*(ik + il + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(i,n)*minusProduct(j,l))/
      (ik*jl*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(l,j)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(k,j)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,j)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(k,n)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,n)*minusCurrent(i,l)*minusProduct(j,l)*minusProduct(l,n))/
      (jl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(i,k)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(k,n)*minusCurrent(i,k)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(i,l)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*plusProduct(i,n)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(i,k))/
      (ik*jl*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,l)*minusProduct(i,k))/(ik*jl) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(j,n))/
      (ik*jl*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(j,n))/
      (ik*jl*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(k,j)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(j,n)*minusCurrent(i,l)*minusProduct(j,n)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,l)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(k,n))/
      (ik*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,n))/
      (ik*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(k,j)*minusProduct(i,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(i,k)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(l,n)*minusCurrent(i,l)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n));
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*sqr(plusProduct(i,n))*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(l,n)) + 
      (Complex(0,2)*sqr(plusProduct(j,n))*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(k,n)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(k,n)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)) + 
      (Complex(0,2)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*minusCurrent(k,j)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(l,n)*minusCurrent(l,j)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*sqr(plusProduct(j,n))*minusCurrent(i,l)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(i,k)*sqr(minusProduct(k,l)))/
      (kl*(jk + jl + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(i,l)*sqr(minusProduct(k,l)))/
      (kl*(jk + jl + kl)*plusProduct(l,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedLeftCurrent(int i, int,
								 int j, int,
								 int k, int g1Hel,
								 int l, int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,k))/
      (jl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,k))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,l))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,j))/
      (jl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(i,j))/
      (jl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,k)*minusProduct(i,j))/
      (jl*(jk + jl + kl)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(i,j)*sqr(minusProduct(i,j)))/
      (jl*(jk + jl + kl)*minusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(i,k)*minusProduct(i,k))/
      (kl*(jk + jl + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(i,l))/
      (jl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(i,l))/
      (kl*(jk + jl + kl)*minusProduct(i,k));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-1)*sqr(plusProduct(i,k))*minusCurrent(i,j)*minusProduct(i,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(l,j)*minusProduct(i,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,k)*minusCurrent(i,j)*sqr(minusProduct(i,l)))/
      (kl*(ik + il + kl)*minusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(k,j)*sqr(minusProduct(i,l)))/
      (kl*(ik + il + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,j)*sqr(minusProduct(i,l)))/
      (kl*(ik + il + kl)*minusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,k)*minusCurrent(i,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(i,j)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,k)*minusCurrent(i,l)*minusProduct(i,l)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(i,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(i,k)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(i,k)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(i,l)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1)*sqr(plusProduct(i,l))*minusCurrent(i,j)*minusProduct(i,k))/
      (kl*(ik + il + kl)*plusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(k,j)*minusProduct(i,k))/
      (kl*(ik + il + kl)*plusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,l)*minusProduct(i,k))/(ik*jl) + 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,j)*minusProduct(i,j)*minusProduct(i,k))/
      (ik*jl*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,l)*minusCurrent(i,j)*sqr(minusProduct(i,k)))/
      (ik*(ik + il + kl)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*minusCurrent(i,j)*sqr(minusProduct(i,k)))/
      (kl*(ik + il + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,j)*sqr(minusProduct(i,k)))/
      (ik*(ik + il + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,j)*sqr(minusProduct(i,k)))/
      (kl*(ik + il + kl)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(l,j)*sqr(minusProduct(i,k)))/
      (kl*(ik + il + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,k)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,j)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,l)*minusCurrent(i,l)*minusProduct(i,j)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(i,k)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,k)*minusCurrent(i,l)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(i,l)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l));
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*plusProduct(i,j)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(i,l)) + 
      (Complex(0,2)*sqr(plusProduct(i,j))*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(k,j)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,2)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)) + 
      (Complex(0,2)*minusCurrent(l,j)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*minusCurrent(k,j)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,l)*minusCurrent(l,j)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(i,j)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(i,j)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*sqr(plusProduct(i,j))*minusCurrent(i,l)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(i,k)*sqr(minusProduct(k,l)))/
      (kl*(jk + jl + kl)*plusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(i,l)*sqr(minusProduct(k,l)))/
      (kl*(jk + jl + kl)*plusProduct(i,l));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralRightCurrent(int i, int,
								    int j, int,
								    int k, int g1Hel,
								    int l, int g2Hel,
								    int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(j,k))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,l))/
      (ik*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,l))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,l)*minusCurrent(j,l)*minusProduct(i,n))/
      (ik*(ik + il + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(k,i)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,n))/
      (ik*(ik + il + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,n))/
      (kl*(ik + il + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,l)*minusCurrent(j,k)*minusProduct(i,n))/
      (ik*(ik + il + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(j,k)*minusProduct(j,n))/
      (ik*jl*minusProduct(l,n)) + 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(l,i)*minusProduct(j,n))/
      (kl*(jk + jl + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,l)*minusCurrent(j,i)*sqr(minusProduct(i,n)))/
      (ik*(ik + il + kl)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,n)*minusProduct(j,n))/
      (ik*jl*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(j,i)*sqr(minusProduct(j,n)))/
      (jl*(jk + jl + kl)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(l,i)*sqr(minusProduct(j,n)))/
      (jl*(jk + jl + kl)*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,k)*minusProduct(k,n))/
      (ik*(ik + il + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(j,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*minusProduct(k,n));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,n)*minusCurrent(j,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(l,k)*minusProduct(j,l))/(ik*jl) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,n)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(i,n)*minusProduct(j,l))/
      (ik*jl*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(l,i)*minusProduct(i,n)*minusProduct(j,l))/
      (ik*jl*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(i,n)*minusCurrent(j,k)*minusProduct(i,n)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(l,i)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(l,n))/
      (jl*(jk + jl + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(k,n)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(k,n)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(l,n))/
      (jl*(jk + jl + kl)*minusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(l,n))/
      (kl*(jk + jl + kl)*minusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(j,k)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,n)*minusCurrent(j,k)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,n)*minusCurrent(j,l)*minusProduct(k,l)*minusProduct(l,n))/
      (kl*(ik + il + kl)*plusProduct(l,n)*minusProduct(k,n));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*plusProduct(i,l)*plusProduct(i,n)*minusCurrent(j,l)*minusProduct(i,k))/
      (ik*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(j,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(i,n)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(i,n))/
      (ik*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(j,n))/
      (ik*jl*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(j,n))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(i,n)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(k,i)*minusProduct(j,n)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(i,n)*minusCurrent(j,k)*minusProduct(i,k)*minusProduct(k,n))/
      (ik*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,n)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(j,k)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(j,l)*plusProduct(l,n)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(k,l)*plusProduct(l,n)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(k,n))/
      (kl*(jk + jl + kl)*plusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(l,n)*minusCurrent(j,k)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(k,l)*minusCurrent(j,l)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(l,n)*minusCurrent(j,l)*minusProduct(k,l)*minusProduct(k,n))/
      (kl*(ik + il + kl)*plusProduct(k,n)*minusProduct(l,n));
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*sqr(plusProduct(i,n))*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,n)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(l,i)*minusProduct(i,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(k,n)) + 
      (Complex(0,2)*sqr(plusProduct(j,n))*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(l,n)) + 
      (Complex(0,2)*sqr(plusProduct(i,n))*minusCurrent(j,k)*minusProduct(i,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*plusProduct(k,n)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(l,n)) - 
      (Complex(0,2)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)) - 
      (Complex(0,2)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(l,n)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(l,n)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(k,n)) + 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(j,k)*sqr(minusProduct(k,l)))/
      (kl*(ik + il + kl)*plusProduct(k,n)) - 
      (Complex(0,2)*plusProduct(i,n)*minusCurrent(j,l)*sqr(minusProduct(k,l)))/
      (kl*(ik + il + kl)*plusProduct(l,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedRightCurrent(int i, int,
								  int j, int,
								  int k, int g1Hel,
								  int l, int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(j,k))/
      (kl*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,l))/
      (ik*(ik + il + kl)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,l))/
      (kl*(ik + il + kl)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(k,i)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(j,k)*minusProduct(i,j))/
      (ik*jl*minusProduct(i,l)) + 
      (Complex(0,2)*sqr(plusProduct(k,l))*minusCurrent(l,i)*minusProduct(i,j))/
      (kl*(jk + jl + kl)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(j,k)*plusProduct(j,l)*minusCurrent(j,i)*sqr(minusProduct(i,j)))/
      (jl*(jk + jl + kl)*minusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(j,l)*plusProduct(k,l)*minusCurrent(l,i)*sqr(minusProduct(i,j)))/
      (jl*(jk + jl + kl)*minusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,k)*minusProduct(i,k))/
      (ik*(ik + il + kl)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(j,k)*minusProduct(i,k))/
      (kl*(ik + il + kl)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(j,l)*minusProduct(i,l))/
      (kl*(ik + il + kl)*minusProduct(i,k));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*sqr(plusProduct(i,k))*minusCurrent(j,i)*minusProduct(i,l))/
      (ik*(ik + il + kl)*plusProduct(i,l)) - 
      (Complex(0,1)*sqr(plusProduct(i,k))*minusCurrent(j,i)*minusProduct(i,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,k)*minusCurrent(j,i)*sqr(minusProduct(i,l)))/
      (kl*(ik + il + kl)*minusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,l)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(i,k)*minusCurrent(j,k)*minusProduct(j,l))/
      (ik*jl*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(l,k)*minusProduct(j,l))/(ik*jl) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(i,j)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,j)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(k,i)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,l)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*minusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,l)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*minusProduct(i,k)) - 
      (Complex(0,2)*sqr(plusProduct(i,k))*minusCurrent(j,k)*minusProduct(k,l))/
      (ik*(ik + il + kl)*plusProduct(i,l)) - 
      (Complex(0,2)*sqr(plusProduct(i,k))*minusCurrent(j,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,k)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)*minusProduct(i,k)) + 
      (Complex(0,1)*plusProduct(i,k)*minusCurrent(j,k)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*minusProduct(i,k)) - 
      (Complex(0,1)*sqr(plusProduct(i,k))*minusCurrent(j,l)*minusProduct(i,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,l)*minusProduct(i,k));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1)*sqr(plusProduct(i,l))*minusCurrent(j,i)*minusProduct(i,k))/
      (kl*(ik + il + kl)*plusProduct(i,k)) - 
      (Complex(0,1)*plusProduct(i,l)*minusCurrent(j,i)*sqr(minusProduct(i,k)))/
      (kl*(ik + il + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(j,l)*minusCurrent(k,i)*minusProduct(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*minusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(l,i)*minusProduct(i,j)*minusProduct(j,k))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,j)*plusProduct(k,l)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(k,l)*minusCurrent(k,i)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*minusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(l,i)*minusProduct(i,k)*minusProduct(j,k))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,k)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,1)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(k,i)*minusProduct(i,k)*minusProduct(j,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*sqr(plusProduct(i,l))*minusCurrent(j,l)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,k)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(j,l)*minusCurrent(j,i)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,l)*plusProduct(k,l)*minusCurrent(k,i)*minusProduct(i,j)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,1)*sqr(plusProduct(i,l))*minusCurrent(j,k)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*plusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,1)*plusProduct(i,l)*minusCurrent(j,l)*minusProduct(i,k)*minusProduct(k,l))/
      (kl*(ik + il + kl)*minusProduct(i,l));
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*sqr(plusProduct(i,j))*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,l)) + 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(j,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(j,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)) - 
      (Complex(0,2)*plusProduct(i,k)*minusCurrent(k,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,l)) - 
      (Complex(0,2)*minusCurrent(l,i)*minusProduct(j,k)*minusProduct(k,l))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,j)*minusCurrent(j,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)) - 
      (Complex(0,2)*minusCurrent(k,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)) - 
      (Complex(0,2)*plusProduct(i,l)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(k,l))/
      (jl*(jk + jl + kl)*plusProduct(i,k)) - 
      (Complex(0,2)*plusProduct(i,l)*minusCurrent(l,i)*minusProduct(j,l)*minusProduct(k,l))/
      (kl*(jk + jl + kl)*plusProduct(i,k));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggLeftCurrent(int q,    int qHel,
								   int qbar, int qbarHel,
								   int g1,   int g1Hel,
								   int g2,   int g2Hel) {
  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<3>(1,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,q);
    LorentzVector<Complex> nj = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,qbar);
    LorentzVector<Complex> nl = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,0);
    LorentzVector<Complex> nlbar = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,1);
    LorentzVector<Complex> fixed = qqbarggFixedLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t()*conj(x1.t())) + real(x1.x()*conj(x1.x())) + real(x1.y()*conj(x1.y())) + real(x1.z()*conj(x1.z()));
    double c2 = 
      real(x2.t()*conj(x2.t())) + real(x2.x()*conj(x2.x())) + real(x2.y()*conj(x2.y())) + real(x2.z()*conj(x2.z()));
    double c3 = 
      real(x3.t()*conj(x3.t())) + real(x3.x()*conj(x3.x())) + real(x3.y()*conj(x3.y())) + real(x3.z()*conj(x3.z()));
    double c4 = 
      real(x4.t()*conj(x4.t())) + real(x4.x()*conj(x4.x())) + real(x4.y()*conj(x4.y())) + real(x4.z()*conj(x4.z()));
    ostream& ncheck = checkStream("qqbarggLeftCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(qqbarggFixedLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarggLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1)+momentum(g2));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggRightCurrent(int q,    int qHel,
								    int qbar, int qbarHel,
								    int g1,   int g1Hel,
								    int g2,   int g2Hel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<3>(2,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,q);
    LorentzVector<Complex> nj = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,qbar);
    LorentzVector<Complex> nl = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,0);
    LorentzVector<Complex> nlbar = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,1);
    LorentzVector<Complex> fixed = qqbarggFixedRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t()*conj(x1.t())) + real(x1.x()*conj(x1.x())) + real(x1.y()*conj(x1.y())) + real(x1.z()*conj(x1.z()));
    double c2 = 
      real(x2.t()*conj(x2.t())) + real(x2.x()*conj(x2.x())) + real(x2.y()*conj(x2.y())) + real(x2.z()*conj(x2.z()));
    double c3 = 
      real(x3.t()*conj(x3.t())) + real(x3.x()*conj(x3.x())) + real(x3.y()*conj(x3.y())) + real(x3.z()*conj(x3.z()));
    double c4 = 
      real(x4.t()*conj(x4.t())) + real(x4.x()*conj(x4.x())) + real(x4.y()*conj(x4.y())) + real(x4.z()*conj(x4.z()));
    ostream& ncheck = checkStream("qqbarggRightCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(qqbarggFixedRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarggRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1)+momentum(g2));
#endif
  return cachedCurrent();


}

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarLeftCurrent(int q,    int qHel,
								      int qbar, int qbarHel,
								      int k,    int kHel,
								      int kbar, int kbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);

  if ( qHel != 1 || qbarHel != 1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  int i = q; int j = qbar; int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusProduct(k,i)*plusProduct(i,l)*minusCurrent(i,j)+
		     minusProduct(i,k)*plusProduct(l,k)*minusCurrent(k,j))/
		    (kl+il+ik)-
		    (minusProduct(j,k)*plusProduct(l,j)*minusCurrent(i,j)+
		     minusProduct(l,k)*plusProduct(l,j)*minusCurrent(i,l))/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }


  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusProduct(l,i)*plusProduct(i,k)*minusCurrent(i,j)+
		     minusProduct(i,l)*plusProduct(k,l)*minusCurrent(l,j))/
		    (kl+il+ik)-
		    (minusProduct(j,l)*plusProduct(k,j)*minusCurrent(i,j)+
		     minusProduct(k,l)*plusProduct(k,j)*minusCurrent(i,k))/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarRightCurrent(int q,    int qHel,
								       int qbar, int qbarHel,
								       int k,    int kHel,
								       int kbar, int kbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);

  if ( qHel != -1 || qbarHel != -1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  int i = q; int j = qbar; int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusProduct(k,i)*plusProduct(i,l)*minusCurrent(j,i)+
		     minusProduct(l,k)*plusProduct(l,i)*minusCurrent(j,l))/
		    (kl+il+ik)-
		    (minusProduct(j,k)*plusProduct(l,j)*minusCurrent(j,i)+
		     minusProduct(j,k)*plusProduct(l,k)*minusCurrent(k,i))/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusProduct(l,i)*plusProduct(i,k)*minusCurrent(j,i)+
		     minusProduct(k,l)*plusProduct(k,i)*minusCurrent(j,k))/
		    (kl+il+ik)-
		    (minusProduct(j,l)*plusProduct(k,j)*minusCurrent(j,i)+
		     minusProduct(j,l)*plusProduct(k,l)*minusCurrent(l,i))/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}


const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftOneLoopCurrent(int q,    int qHel,
									int qbar, int qbarHel) {


  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  const LorentzVector<Complex>& tree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);
  if ( getCurrent(hash<1>(1,2,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(0.5*CF*(3.*log(abs(invariant(q,qbar))) - 8.)*tree);
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarLeftOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightOneLoopCurrent(int q,    int qHel,
									 int qbar, int qbarHel) { 


  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  const LorentzVector<Complex>& tree = qqbarRightCurrent(q,qHel,qbar,qbarHel);
  if ( getCurrent(hash<1>(2,2,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(0.5*CF*(3.*log(abs(invariant(q,qbar))) - 8.)*tree);
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarRightOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  return cachedCurrent();

}

// ln(s(a+i0))
inline Complex log(double s, double a) {
  return
    s < 0. ?
    Complex(log(abs(a)),-pi*theta(a)) :
    Complex(log(abs(a)),pi*theta(-a));
}

// ln(s(a+i0)/(b+i0))
inline Complex log(double s, double a, double b) {
  return
    s < 0. ?
    Complex(log(abs(a/b)),-pi*theta(a/b)*sign(b-a)) :
    Complex(log(abs(a/b)),pi*theta(-a/b)*sign(b-a));
}

// Li2(-(a+i0)/(b+i0))
inline Complex Li2(double a, double b) {
  if ( -a/b < 1. )
    return Complex(Herwig::Math::ReLi2(-a/b),0.0);
  return Complex(Herwig::Math::ReLi2(-a/b),-pi*log(-a/b)*sign(b-a));
}

Complex MatchboxCurrents::box6(int i, int j, int k) {

  const double sij = invariant(i,j);
  const double sik = invariant(i,k);
  const double sjk = invariant(j,k);

  return
    -( Li2(sik+sjk,sij) + Li2(sik+sij,sjk) + 0.5*csqr(log(1.,sij,sjk)) + sqr(pi)/6. )/8.;
      
}

void MatchboxCurrents::qqbargLoopCoefficients(int i, int j, int k) {

  // use a dummy cache entry to check if we need to get some work done
  static Complex dummy;

  if ( getAmplitude(hash<5>(1,2,i,0,j,0,k,0)) ) {
    dummy = 0.;
    cacheAmplitude(dummy);
    cachedAmplitude();
  } else {
    cachedAmplitude();
    return;
  }

  qqbargLoops.resize(13);

  // get the transcendentals
  const double ij = invariant(i,j);  
  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const double Q2 = ij + ik + jk;
  // checked for LEP that virtuals + I operator are mu2 independent
  //double xmu2 = 10*GeV2/sqr(amplitudeScale());
  const double xmu2 = 1.;

  const Complex Lijk = log(1.,-xmu2/Q2);

  const Complex Lij = log(1.,Q2,ij);
  const Complex Lik = log(1.,Q2,ik);
  const Complex Ljk = log(1.,Q2,jk);

  const Complex Box6ijk = box6(i,j,k);
  const Complex Box6ikj = box6(i,k,j);
  const Complex Box6jik = box6(j,i,k);

  // get the coefficients

  qqbargLoops[0] = 
    (2*CF*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (32*CA*Box6ijk*sqr(ij))/
    (sqr(ij + ik)*(ij + jk)) + 
    (64*CF*Box6ijk*sqr(ij))/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*sqr(ij))/
    (sqr(ij + ik)*(ij + jk)) + 
    (16*CF*Box6jik*sqr(ij))/
    (sqr(ij + ik)*(ij + jk)) + 
    (2*CA*Lij*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Lij*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lik*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Ljk*sqr(ij))/(sqr(ij + ik)*(ij + jk)) - 
    (16*CA*Box6ijk*pow(ij,3))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (32*CF*Box6ijk*pow(ij,3))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (CA*Lij*pow(ij,3))/
    (ik*sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*pow(ij,3))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (2*CF*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (16*CA*Box6ijk*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) + 
    (32*CF*Box6ijk*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (16*CA*Box6jik*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) + 
    (32*CF*Box6jik*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) + 
    (CA*Lij*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (2*CA*Lik*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Lik*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Ljk*ij*ik)/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*sqr(ik))/
    (sqr(ij + ik)*(ij + jk)) + 
    (16*CF*Box6jik*sqr(ik))/
    (sqr(ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ik))/(sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lik*sqr(ik))/(sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*pow(ij,3))/
    (sqr(ij + ik)*jk*(ij + jk)) + 
    (16*CF*Box6jik*pow(ij,3))/
    (sqr(ij + ik)*jk*(ij + jk)) - 
    (16*CA*Box6jik*sqr(ij)*ik)/
    (sqr(ij + ik)*jk*(ij + jk)) + 
    (32*CF*Box6jik*sqr(ij)*ik)/
    (sqr(ij + ik)*jk*(ij + jk)) - 
    (8*CA*Box6jik*ij*sqr(ik))/
    (sqr(ij + ik)*jk*(ij + jk)) + 
    (16*CF*Box6jik*ij*sqr(ik))/
    (sqr(ij + ik)*jk*(ij + jk)) + 
    (2*CF*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (40*CA*Box6ijk*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (80*CF*Box6ijk*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (2*CA*Lij*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Lij*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (CA*Lik*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Lik*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (12*CF*Ljk*ij*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*pow(ij,3)*jk)/
    (sqr(ik)*sqr(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*pow(ij,3)*jk)/
    (sqr(ik)*sqr(ij + ik)*(ij + jk)) - 
    (32*CA*Box6ijk*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (64*CF*Box6ijk*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (CA*Lij*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (CA*Ljk*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*sqr(ij)*jk)/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (2*CF*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (16*CA*Box6ijk*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (32*CF*Box6ijk*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (48*CA*Box6ikj*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (CA*Lij*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (2*CA*Lik*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CF*Lik*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CF*Ljk*ik*jk)/
    (sqr(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*sqr(ik)*jk)/
    (ij*sqr(ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ik)*jk)/
    (ij*sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Lik*sqr(ik)*jk)/
    (ij*sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(jk))/
    (sqr(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(jk))/
    (sqr(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*sqr(jk))/
    (sqr(ij + ik)*(ij + jk)) - 
    (8*CF*Ljk*sqr(jk))/(sqr(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(ij)*sqr(jk))/
    (sqr(ik)*sqr(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(ij)*sqr(jk))/
    (sqr(ik)*sqr(ij + ik)*(ij + jk)) - 
    (16*CA*Box6ijk*ij*sqr(jk))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (32*CF*Box6ijk*ij*sqr(jk))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (CA*Ljk*ij*sqr(jk))/
    (ik*sqr(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*ij*sqr(jk))/
    (ik*sqr(ij + ik)*(ij + jk)) + 
    (48*CA*Box6ikj*ik*sqr(jk))/
    (ij*sqr(ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik*sqr(jk))/
    (ij*sqr(ij + ik)*(ij + jk)) - 
    (4*CF*Ljk*ik*sqr(jk))/
    (ij*sqr(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (sqr(ij)*sqr(ij + ik)*(ij + jk));

  qqbargLoops[1] = 
    (-2*CF*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6ijk*sqr(ij))/
    ((ij + ik)*sqr(ij + jk)) - 
    (16*CF*Box6ijk*sqr(ij))/
    ((ij + ik)*sqr(ij + jk)) + 
    (32*CA*Box6jik*sqr(ij))/
    ((ij + ik)*sqr(ij + jk)) - 
    (64*CF*Box6jik*sqr(ij))/
    ((ij + ik)*sqr(ij + jk)) - 
    (2*CA*Lij*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Lij*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Lik*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (CA*Ljk*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (2*CF*Ljk*sqr(ij))/((ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6ijk*pow(ij,3))/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (16*CF*Box6ijk*pow(ij,3))/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (2*CF*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) - 
    (24*CA*Box6ikj*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) + 
    (40*CA*Box6jik*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) - 
    (80*CF*Box6jik*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) - 
    (2*CA*Lij*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Lij*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) + 
    (12*CF*Lik*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) + 
    (CA*Ljk*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Ljk*ij*ik)/
    ((ij + ik)*sqr(ij + jk)) - 
    (24*CA*Box6ikj*sqr(ik))/
    ((ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6jik*sqr(ik))/
    ((ij + ik)*sqr(ij + jk)) - 
    (16*CF*Box6jik*sqr(ik))/
    ((ij + ik)*sqr(ij + jk)) + 
    (8*CF*Lik*sqr(ik))/((ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6jik*pow(ij,3)*ik)/
    ((ij + ik)*sqr(jk)*sqr(ij + jk)) - 
    (16*CF*Box6jik*pow(ij,3)*ik)/
    ((ij + ik)*sqr(jk)*sqr(ij + jk)) + 
    (8*CA*Box6jik*sqr(ij)*sqr(ik))/
    ((ij + ik)*sqr(jk)*sqr(ij + jk)) - 
    (16*CF*Box6jik*sqr(ij)*sqr(ik))/
    ((ij + ik)*sqr(jk)*sqr(ij + jk)) + 
    (16*CA*Box6jik*pow(ij,3))/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (32*CF*Box6jik*pow(ij,3))/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (CA*Lij*pow(ij,3))/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (2*CF*Lij*pow(ij,3))/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (32*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (64*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (CA*Lij*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (2*CF*Lij*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (CA*Lik*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (2*CF*Lik*sqr(ij)*ik)/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (16*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (32*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (CA*Lik*ij*sqr(ik))/
    ((ij + ik)*jk*sqr(ij + jk)) + 
    (2*CF*Lik*ij*sqr(ik))/
    ((ij + ik)*jk*sqr(ij + jk)) - 
    (2*CF*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (16*CA*Box6ijk*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (32*CF*Box6ijk*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (16*CA*Box6jik*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (32*CF*Box6jik*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (CA*Lij*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (2*CF*Lij*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Lik*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (2*CA*Ljk*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (4*CF*Ljk*ij*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (16*CA*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (32*CF*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (2*CF*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (48*CA*Box6ikj*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (16*CA*Box6jik*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (32*CF*Box6jik*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (CA*Lij*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (2*CF*Lij*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (CA*Lik*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (8*CF*Lik*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (2*CA*Ljk*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) + 
    (8*CF*Ljk*ik*jk)/
    ((ij + ik)*sqr(ij + jk)) - 
    (48*CA*Box6ikj*sqr(ik)*jk)/
    (ij*(ij + ik)*sqr(ij + jk)) + 
    (CA*Lik*sqr(ik)*jk)/
    (ij*(ij + ik)*sqr(ij + jk)) + 
    (4*CF*Lik*sqr(ik)*jk)/
    (ij*(ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6ijk*sqr(jk))/
    ((ij + ik)*sqr(ij + jk)) - 
    (16*CF*Box6ijk*sqr(jk))/
    ((ij + ik)*sqr(ij + jk)) + 
    (CA*Ljk*sqr(jk))/((ij + ik)*sqr(ij + jk)) + 
    (2*CF*Ljk*sqr(jk))/((ij + ik)*sqr(ij + jk)) + 
    (8*CA*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (16*CF*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*sqr(ij + jk)) - 
    (24*CA*Box6ikj*ik*sqr(jk))/
    (ij*(ij + ik)*sqr(ij + jk)) + 
    (CA*Ljk*ik*sqr(jk))/
    (ij*(ij + ik)*sqr(ij + jk)) + 
    (4*CF*Ljk*ik*sqr(jk))/
    (ij*(ij + ik)*sqr(ij + jk)) - 
    (24*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (sqr(ij)*(ij + ik)*sqr(ij + jk));

  qqbargLoops[2] = -3*CF*Lijk +
    (-4*CA*Box6jik*pow(ij,3))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CF*Box6jik*pow(ij,3))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lij*pow(ij,3))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lij*pow(ij,3))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (16*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lij*sqr(ij)*ik)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lij*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*sqr(ij)*ik)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lik*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (4*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*ij*sqr(ik))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lik*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (4*CA*Box6jik*pow(ij,3)*ik)/
    ((ij + ik)*jk*(ij + jk)*
     (ik + jk)) + (8*CF*Box6jik*pow(ij,3)*ik)/
    ((ij + ik)*jk*(ij + jk)*
     (ik + jk)) - (4*CA*Box6jik*sqr(ij)*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)*
     (ik + jk)) + (8*CF*Box6jik*sqr(ij)*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)*
     (ik + jk)) + (CA*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (12*CA*Box6ijk*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (24*CF*Box6ijk*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (4*CA*Box6jik*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CF*Box6jik*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*sqr(ij)*jk)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lik*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (CA*Ljk*sqr(ij)*jk)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (CF*Ljk*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (4*CA*Box6ijk*pow(ij,3)*jk)/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (8*CF*Box6ijk*pow(ij,3)*jk)/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (CA*Lij*pow(ij,3)*jk)/
    (2.*ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (CF*Lij*pow(ij,3)*jk)/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (2*CA*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (18*CF*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (16*CA*Box6ijk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (32*CF*Box6ijk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (28*CA*Box6ikj*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (4*CA*Box6jik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CF*Box6jik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lij*ij*ik*jk)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CF*Lij*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (CF*Lik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (CA*Ljk*ij*ik*jk)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (3*CF*Ljk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (20*CA*Box6ikj*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*sqr(ik)*jk)/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (CA*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (12*CA*Box6ijk*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (24*CF*Box6ijk*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (20*CA*Box6ikj*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (CA*Lij*ij*sqr(jk))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (CF*Lij*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*ij*sqr(jk))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CA*Ljk*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (4*CF*Ljk*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (4*CA*Box6ijk*pow(ij,3)*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ij + jk)*
     (ik + jk)) - (8*CF*Box6ijk*pow(ij,3)*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ij + jk)*
     (ik + jk)) + (8*CA*Box6ijk*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (16*CF*Box6ijk*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (CA*Lij*sqr(ij)*sqr(jk))/
    (2.*ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (CF*Lij*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (CA*Ljk*sqr(ij)*sqr(jk))/
    (2.*ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (CF*Ljk*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (CA*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (32*CA*Box6ikj*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*Lik*ik*sqr(jk))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    - (CA*Ljk*ik*sqr(jk))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (3*CF*Ljk*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (12*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)*
     (ik + jk)) - (12*CA*Box6ikj*pow(jk,3))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (CA*Ljk*pow(jk,3))/
    (2.*(ij + ik)*(ij + jk)*(ik + jk))
    + (3*CF*Ljk*pow(jk,3))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (4*CA*Box6ijk*sqr(ij)*pow(jk,3))/
    (sqr(ik)*(ij + ik)*(ij + jk)*
     (ik + jk)) - (8*CF*Box6ijk*sqr(ij)*pow(jk,3))/
    (sqr(ik)*(ij + ik)*(ij + jk)*
     (ik + jk)) + (4*CA*Box6ijk*ij*pow(jk,3))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (8*CF*Box6ijk*ij*pow(jk,3))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (CA*Ljk*ij*pow(jk,3))/
    (2.*ik*(ij + ik)*(ij + jk)*
     (ik + jk)) + (CF*Ljk*ij*pow(jk,3))/
    (ik*(ij + ik)*(ij + jk)*
     (ik + jk)) - (12*CA*Box6ikj*ik*pow(jk,3))/
    (ij*(ij + ik)*(ij + jk)*
     (ik + jk));

  qqbargLoops[3] = 3*CF*Lijk +
    (8*CF*sqr(ij))/((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(ij))/((ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(ij))/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ikj*sqr(ij))/((ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*sqr(ij))/((ij + ik)*(ij + jk)) + 
    (16*CF*Box6jik*sqr(ij))/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lij*sqr(ij))/(2.*(ij + ik)*(ij + jk)) - 
    (CF*Lij*sqr(ij))/((ij + ik)*(ij + jk)) + 
    (8*CF*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ikj*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (12*CA*Box6jik*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (24*CF*Box6jik*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lij*ij*ik)/
    (2.*(ij + ik)*(ij + jk)) - 
    (CF*Lij*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lik*ij*ik)/
    (2.*(ij + ik)*(ij + jk)) - 
    (CF*Lik*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (4*CA*Box6jik*sqr(ik))/((ij + ik)*(ij + jk)) + 
    (8*CF*Box6jik*sqr(ik))/((ij + ik)*(ij + jk)) + 
    (CA*Lik*sqr(ik))/(2.*(ij + ik)*(ij + jk)) - 
    (CF*Lik*sqr(ik))/((ij + ik)*(ij + jk)) - 
    (4*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (4*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CF*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (12*CA*Box6ijk*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (24*CF*Box6ijk*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ikj*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (16*CF*Box6jik*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lij*ij*jk)/
    (2.*(ij + ik)*(ij + jk)) - 
    (CF*Lij*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Ljk*ij*jk)/
    (2.*(ij + ik)*(ij + jk)) - 
    (CF*Ljk*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (4*CA*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (8*CF*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (8*CF*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (4*CA*Box6ikj*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6jik*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (16*CF*Box6jik*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lij*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) - 
    (CF*Lij*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Lik*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) + 
    (2*CF*Lik*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (CA*Ljk*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) + 
    (2*CF*Ljk*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (12*CA*Box6ikj*sqr(ik)*jk)/
    (ij*(ij + ik)*(ij + jk)) + 
    (CA*Lik*sqr(ik)*jk)/
    (2.*ij*(ij + ik)*(ij + jk)) + 
    (2*CF*Lik*sqr(ik)*jk)/
    (ij*(ij + ik)*(ij + jk)) - 
    (4*CA*Box6ijk*sqr(jk))/((ij + ik)*(ij + jk)) + 
    (8*CF*Box6ijk*sqr(jk))/((ij + ik)*(ij + jk)) + 
    (CA*Ljk*sqr(jk))/(2.*(ij + ik)*(ij + jk)) - 
    (CF*Ljk*sqr(jk))/((ij + ik)*(ij + jk)) - 
    (4*CA*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) + 
    (8*CF*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) - 
    (12*CA*Box6ikj*ik*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) + 
    (CA*Ljk*ik*sqr(jk))/
    (2.*ij*(ij + ik)*(ij + jk)) + 
    (2*CF*Ljk*ik*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) - 
    (12*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (sqr(ij)*(ij + ik)*(ij + jk));

  qqbargLoops[4] = -3*CF*Lijk +
    (-8*CF*sqr(ij))/((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ijk*sqr(ij))/((ij + ik)*(ij + jk)) - 
    (16*CF*Box6ijk*sqr(ij))/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ikj*sqr(ij))/((ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*sqr(ij))/((ij + ik)*(ij + jk)) - 
    (16*CF*Box6jik*sqr(ij))/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lij*sqr(ij))/(2.*(ij + ik)*(ij + jk)) + 
    (CF*Lij*sqr(ij))/((ij + ik)*(ij + jk)) - 
    (8*CF*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ijk*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (16*CF*Box6ijk*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ikj*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (12*CA*Box6jik*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (24*CF*Box6jik*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lij*ij*ik)/
    (2.*(ij + ik)*(ij + jk)) + 
    (CF*Lij*ij*ik)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lik*ij*ik)/
    (2.*(ij + ik)*(ij + jk)) + 
    (CF*Lik*ij*ik)/
    ((ij + ik)*(ij + jk)) + 
    (4*CA*Box6jik*sqr(ik))/((ij + ik)*(ij + jk)) - 
    (8*CF*Box6jik*sqr(ik))/((ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ik))/(2.*(ij + ik)*(ij + jk)) + 
    (CF*Lik*sqr(ik))/((ij + ik)*(ij + jk)) + 
    (4*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (8*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ij + jk)) + 
    (4*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) - 
    (8*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) - 
    (8*CF*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (12*CA*Box6ijk*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (24*CF*Box6ijk*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ikj*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (16*CF*Box6jik*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lij*ij*jk)/
    (2.*(ij + ik)*(ij + jk)) + 
    (CF*Lij*ij*jk)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Ljk*ij*jk)/
    (2.*(ij + ik)*(ij + jk)) + 
    (CF*Ljk*ij*jk)/
    ((ij + ik)*(ij + jk)) + 
    (4*CA*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ij + jk)) - 
    (8*CF*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ij + jk)) - 
    (8*CF*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6ijk*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (16*CF*Box6ijk*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (4*CA*Box6ikj*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (16*CF*Box6jik*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lij*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) + 
    (CF*Lij*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Lik*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) - 
    (2*CF*Lik*ik*jk)/
    ((ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik*jk)/
    (2.*(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*ik*jk)/
    ((ij + ik)*(ij + jk)) + 
    (12*CA*Box6ikj*sqr(ik)*jk)/
    (ij*(ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ik)*jk)/
    (2.*ij*(ij + ik)*(ij + jk)) - 
    (2*CF*Lik*sqr(ik)*jk)/
    (ij*(ij + ik)*(ij + jk)) + 
    (4*CA*Box6ijk*sqr(jk))/((ij + ik)*(ij + jk)) - 
    (8*CF*Box6ijk*sqr(jk))/((ij + ik)*(ij + jk)) - 
    (CA*Ljk*sqr(jk))/(2.*(ij + ik)*(ij + jk)) + 
    (CF*Ljk*sqr(jk))/((ij + ik)*(ij + jk)) + 
    (4*CA*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) - 
    (8*CF*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) + 
    (12*CA*Box6ikj*ik*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik*sqr(jk))/
    (2.*ij*(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*ik*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) + 
    (12*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (sqr(ij)*(ij + ik)*(ij + jk));

  qqbargLoops[5] = 3*CF*Lijk +
    (-4*CA*Box6jik*sqr(ij))/((ij + ik)*(ik + jk)) + 
    (8*CF*Box6jik*sqr(ij))/((ij + ik)*(ik + jk)) + 
    (CA*Lij*sqr(ij))/(2.*(ij + ik)*(ik + jk)) - 
    (CF*Lij*sqr(ij))/((ij + ik)*(ik + jk)) - 
    (CA*ij*ik)/((ij + ik)*(ik + jk)) + 
    (9*CF*ij*ik)/
    ((ij + ik)*(ik + jk)) - 
    (8*CA*Box6ijk*ij*ik)/
    ((ij + ik)*(ik + jk)) + 
    (16*CF*Box6ijk*ij*ik)/
    ((ij + ik)*(ik + jk)) + 
    (8*CA*Box6ikj*ij*ik)/
    ((ij + ik)*(ik + jk)) - 
    (4*CA*Box6jik*ij*ik)/
    ((ij + ik)*(ik + jk)) + 
    (8*CF*Box6jik*ij*ik)/
    ((ij + ik)*(ik + jk)) + 
    (CA*Lij*ij*ik)/
    (2.*(ij + ik)*(ik + jk)) - 
    (CF*Lij*ij*ik)/
    ((ij + ik)*(ik + jk)) + 
    (CA*Lik*ij*ik)/
    (2.*(ij + ik)*(ik + jk)) - 
    (CF*Lik*ij*ik)/
    ((ij + ik)*(ik + jk)) - 
    (CA*sqr(ik))/((ij + ik)*(ik + jk)) + 
    (9*CF*sqr(ik))/((ij + ik)*(ik + jk)) - 
    (8*CA*Box6ijk*sqr(ik))/((ij + ik)*(ik + jk)) + 
    (16*CF*Box6ijk*sqr(ik))/
    ((ij + ik)*(ik + jk)) + 
    (8*CA*Box6ikj*sqr(ik))/((ij + ik)*(ik + jk)) + 
    (CA*Lik*sqr(ik))/(2.*(ij + ik)*(ik + jk)) - 
    (CF*Lik*sqr(ik))/((ij + ik)*(ik + jk)) - 
    (4*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ik + jk)) + 
    (8*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*jk*(ik + jk)) - 
    (4*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ik + jk)) + 
    (8*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*jk*(ik + jk)) - 
    (CA*ij*jk)/((ij + ik)*(ik + jk)) + 
    (9*CF*ij*jk)/
    ((ij + ik)*(ik + jk)) - 
    (4*CA*Box6ijk*ij*jk)/
    ((ij + ik)*(ik + jk)) + 
    (8*CF*Box6ijk*ij*jk)/
    ((ij + ik)*(ik + jk)) + 
    (8*CA*Box6ikj*ij*jk)/
    ((ij + ik)*(ik + jk)) - 
    (CA*Lij*ij*jk)/
    (2.*(ij + ik)*(ik + jk)) + 
    (CF*Lij*ij*jk)/
    ((ij + ik)*(ik + jk)) + 
    (CA*Lik*ij*jk)/
    (2.*(ij + ik)*(ik + jk)) - 
    (CF*Lik*ij*jk)/
    ((ij + ik)*(ik + jk)) - 
    (CA*Ljk*ij*jk)/
    (2.*(ij + ik)*(ik + jk)) + 
    (CF*Ljk*ij*jk)/
    ((ij + ik)*(ik + jk)) + 
    (4*CA*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ik + jk)) - 
    (8*CF*Box6ijk*sqr(ij)*jk)/
    (ik*(ij + ik)*(ik + jk)) - 
    (CA*Lij*sqr(ij)*jk)/
    (2.*ik*(ij + ik)*(ik + jk)) + 
    (CF*Lij*sqr(ij)*jk)/
    (ik*(ij + ik)*(ik + jk)) - 
    (CA*ik*jk)/((ij + ik)*(ik + jk)) + 
    (9*CF*ik*jk)/
    ((ij + ik)*(ik + jk)) - 
    (8*CA*Box6ijk*ik*jk)/
    ((ij + ik)*(ik + jk)) + 
    (16*CF*Box6ijk*ik*jk)/
    ((ij + ik)*(ik + jk)) + 
    (20*CA*Box6ikj*ik*jk)/
    ((ij + ik)*(ik + jk)) + 
    (CA*Lik*ik*jk)/
    (2.*(ij + ik)*(ik + jk)) - 
    (CF*Lik*ik*jk)/
    ((ij + ik)*(ik + jk)) - 
    (CA*Ljk*ik*jk)/
    (2.*(ij + ik)*(ik + jk)) - 
    (2*CF*Ljk*ik*jk)/
    ((ij + ik)*(ik + jk)) + 
    (12*CA*Box6ikj*sqr(ik)*jk)/
    (ij*(ij + ik)*(ik + jk)) + 
    (12*CA*Box6ikj*sqr(jk))/
    ((ij + ik)*(ik + jk)) - 
    (CA*Ljk*sqr(jk))/(2.*(ij + ik)*(ik + jk)) - 
    (2*CF*Ljk*sqr(jk))/((ij + ik)*(ik + jk)) + 
    (4*CA*Box6ijk*sqr(ij)*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ik + jk)) - 
    (8*CF*Box6ijk*sqr(ij)*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ik + jk)) + 
    (4*CA*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ik + jk)) - 
    (8*CF*Box6ijk*ij*sqr(jk))/
    (ik*(ij + ik)*(ik + jk)) - 
    (CA*Ljk*ij*sqr(jk))/
    (2.*ik*(ij + ik)*(ik + jk)) + 
    (CF*Ljk*ij*sqr(jk))/
    (ik*(ij + ik)*(ik + jk)) + 
    (12*CA*Box6ikj*ik*sqr(jk))/
    (ij*(ij + ik)*(ik + jk));

  qqbargLoops[6] = 
    (-2*CF*ij)/sqr(ij + ik) + 
    (32*CA*Box6ijk*ij)/sqr(ij + ik) - 
    (64*CF*Box6ijk*ij)/sqr(ij + ik) - 
    (4*CA*Lij*ij)/sqr(ij + ik) + 
    (8*CF*Lij*ij)/sqr(ij + ik) + 
    (4*CF*Ljk*ij)/sqr(ij + ik) + 
    (16*CA*Box6ijk*sqr(ij))/(ik*sqr(ij + ik)) - 
    (32*CF*Box6ijk*sqr(ij))/(ik*sqr(ij + ik)) - 
    (2*CA*Lij*sqr(ij))/(ik*sqr(ij + ik)) + 
    (4*CF*Lij*sqr(ij))/(ik*sqr(ij + ik)) - 
    (2*CF*ik)/sqr(ij + ik) + 
    (16*CA*Box6ijk*ik)/sqr(ij + ik) - 
    (32*CF*Box6ijk*ik)/sqr(ij + ik) - 
    (2*CA*Lij*ik)/sqr(ij + ik) + 
    (4*CF*Lij*ik)/sqr(ij + ik) + 
    (4*CF*Ljk*ik)/sqr(ij + ik) + 
    (16*CA*Box6ijk*jk)/sqr(ij + ik) - 
    (32*CF*Box6ijk*jk)/sqr(ij + ik) - 
    (2*CA*Ljk*jk)/sqr(ij + ik) + 
    (6*CF*Ljk*jk)/sqr(ij + ik) + 
    (16*CA*Box6ijk*sqr(ij)*jk)/
    (sqr(ik)*sqr(ij + ik)) - 
    (32*CF*Box6ijk*sqr(ij)*jk)/
    (sqr(ik)*sqr(ij + ik)) + 
    (32*CA*Box6ijk*ij*jk)/(ik*sqr(ij + ik)) - 
    (64*CF*Box6ijk*ij*jk)/(ik*sqr(ij + ik)) - 
    (2*CA*Ljk*ij*jk)/(ik*sqr(ij + ik)) + 
    (4*CF*Ljk*ij*jk)/(ik*sqr(ij + ik));

  qqbargLoops[7] = 
    (8*CA*Box6jik*ij)/((ij + ik)*(ij + jk)) - 
    (16*CF*Box6jik*ij)/((ij + ik)*(ij + jk)) + 
    (CA*Lij*ij)/((ij + ik)*(ij + jk)) - 
    (2*CF*Lij*ij)/((ij + ik)*(ij + jk)) + 
    (CA*Lik*ij)/((ij + ik)*(ij + jk)) + 
    (2*CF*Lik*ij)/((ij + ik)*(ij + jk)) + 
    (CA*Lij*sqr(ij))/
    (ik*(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*sqr(ij))/
    (ik*(ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*ik)/((ij + ik)*(ij + jk)) - 
    (16*CF*Box6jik*ik)/((ij + ik)*(ij + jk)) + 
    (CA*Lik*ik)/((ij + ik)*(ij + jk)) + 
    (2*CF*Lik*ik)/((ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*sqr(ij))/
    ((ij + ik)*jk*(ij + jk)) - 
    (16*CF*Box6jik*sqr(ij))/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CA*Box6jik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (16*CF*Box6jik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (24*CA*Box6ikj*jk)/((ij + ik)*(ij + jk)) + 
    (CA*Lij*jk)/((ij + ik)*(ij + jk)) - 
    (2*CF*Lij*jk)/((ij + ik)*(ij + jk)) + 
    (CA*Lik*jk)/((ij + ik)*(ij + jk)) + 
    (4*CF*Lik*jk)/((ij + ik)*(ij + jk)) + 
    (CA*Ljk*jk)/((ij + ik)*(ij + jk)) + 
    (4*CF*Ljk*jk)/((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(ij)*jk)/
    (sqr(ik)*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(ij)*jk)/
    (sqr(ik)*(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (CA*Lij*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) - 
    (2*CF*Lij*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (CA*Ljk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) - 
    (24*CA*Box6ikj*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) + 
    (CA*Lik*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) + 
    (4*CF*Lik*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) - 
    (24*CA*Box6ikj*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) + 
    (CA*Ljk*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) + 
    (4*CF*Ljk*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*ij*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ij*sqr(jk))/
    (sqr(ik)*(ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) + 
    (CA*Ljk*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*sqr(jk))/
    (ik*(ij + ik)*(ij + jk)) - 
    (24*CA*Box6ikj*ik*sqr(jk))/
    (sqr(ij)*(ij + ik)*(ij + jk));

  qqbargLoops[8] = 
    (-8*CA*Box6ijk*ij)/((ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ij)/((ij + ik)*(ij + jk)) - 
    (CA*Lij*ij)/((ij + ik)*(ij + jk)) + 
    (2*CF*Lij*ij)/((ij + ik)*(ij + jk)) - 
    (CA*Ljk*ij)/((ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*ij)/((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*sqr(ij))/
    (ik*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*sqr(ij))/
    (ik*(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*ik)/((ij + ik)*(ij + jk)) - 
    (CA*Lij*ik)/((ij + ik)*(ij + jk)) + 
    (2*CF*Lij*ik)/((ij + ik)*(ij + jk)) - 
    (CA*Lik*ik)/((ij + ik)*(ij + jk)) - 
    (4*CF*Lik*ik)/((ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik)/((ij + ik)*(ij + jk)) - 
    (4*CF*Ljk*ik)/((ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*sqr(ik))/
    (ij*(ij + ik)*(ij + jk)) - 
    (CA*Lik*sqr(ik))/
    (ij*(ij + ik)*(ij + jk)) - 
    (4*CF*Lik*sqr(ik))/
    (ij*(ij + ik)*(ij + jk)) + 
    (8*CA*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*sqr(jk)*(ij + jk)) - 
    (16*CF*Box6jik*sqr(ij)*ik)/
    ((ij + ik)*sqr(jk)*(ij + jk)) + 
    (8*CA*Box6jik*ij*sqr(ik))/
    ((ij + ik)*sqr(jk)*(ij + jk)) - 
    (16*CF*Box6jik*ij*sqr(ik))/
    ((ij + ik)*sqr(jk)*(ij + jk)) - 
    (CA*Lij*sqr(ij))/
    ((ij + ik)*jk*(ij + jk)) + 
    (2*CF*Lij*sqr(ij))/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CA*Box6jik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (16*CF*Box6jik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (CA*Lij*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) + 
    (2*CF*Lij*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) - 
    (CA*Lik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) + 
    (2*CF*Lik*ij*ik)/
    ((ij + ik)*jk*(ij + jk)) + 
    (8*CA*Box6jik*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) - 
    (16*CF*Box6jik*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) - 
    (CA*Lik*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) + 
    (2*CF*Lik*sqr(ik))/
    ((ij + ik)*jk*(ij + jk)) - 
    (8*CA*Box6ijk*jk)/((ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*jk)/((ij + ik)*(ij + jk)) - 
    (CA*Ljk*jk)/((ij + ik)*(ij + jk)) - 
    (2*CF*Ljk*jk)/((ij + ik)*(ij + jk)) - 
    (8*CA*Box6ijk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (16*CF*Box6ijk*ij*jk)/
    (ik*(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) - 
    (CA*Ljk*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) - 
    (4*CF*Ljk*ik*jk)/
    (ij*(ij + ik)*(ij + jk)) + 
    (24*CA*Box6ikj*sqr(ik)*jk)/
    (sqr(ij)*(ij + ik)*(ij + jk));

  qqbargLoops[9] = 
    (2*CF*ij)/sqr(ij + jk) - 
    (32*CA*Box6jik*ij)/sqr(ij + jk) + 
    (64*CF*Box6jik*ij)/sqr(ij + jk) + 
    (4*CA*Lij*ij)/sqr(ij + jk) - 
    (8*CF*Lij*ij)/sqr(ij + jk) - 
    (4*CF*Lik*ij)/sqr(ij + jk) - 
    (16*CA*Box6jik*ik)/sqr(ij + jk) + 
    (32*CF*Box6jik*ik)/sqr(ij + jk) + 
    (2*CA*Lik*ik)/sqr(ij + jk) - 
    (6*CF*Lik*ik)/sqr(ij + jk) - 
    (16*CA*Box6jik*sqr(ij)*ik)/
    (sqr(jk)*sqr(ij + jk)) + 
    (32*CF*Box6jik*sqr(ij)*ik)/
    (sqr(jk)*sqr(ij + jk)) - 
    (16*CA*Box6jik*sqr(ij))/(jk*sqr(ij + jk)) + 
    (32*CF*Box6jik*sqr(ij))/(jk*sqr(ij + jk)) + 
    (2*CA*Lij*sqr(ij))/(jk*sqr(ij + jk)) - 
    (4*CF*Lij*sqr(ij))/(jk*sqr(ij + jk)) - 
    (32*CA*Box6jik*ij*ik)/(jk*sqr(ij + jk)) + 
    (64*CF*Box6jik*ij*ik)/(jk*sqr(ij + jk)) + 
    (2*CA*Lik*ij*ik)/(jk*sqr(ij + jk)) - 
    (4*CF*Lik*ij*ik)/(jk*sqr(ij + jk)) + 
    (2*CF*jk)/sqr(ij + jk) - 
    (16*CA*Box6jik*jk)/sqr(ij + jk) + 
    (32*CF*Box6jik*jk)/sqr(ij + jk) + 
    (2*CA*Lij*jk)/sqr(ij + jk) - 
    (4*CF*Lij*jk)/sqr(ij + jk) - 
    (4*CF*Lik*jk)/sqr(ij + jk);

  qqbargLoops[10] = 
    (-8*CA*Box6ijk*sqr(ij)*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (16*CF*Box6ijk*sqr(ij)*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CA*Lij*sqr(ij)*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (4*CF*Lij*sqr(ij)*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (CA*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CF*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (8*CA*Box6ijk*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (16*CF*Box6ijk*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (3*CA*Lij*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (6*CF*Lij*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (CA*Ljk*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (2*CF*Ljk*ij*ik*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (CA*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CF*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (CA*Lij*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (2*CF*Lij*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) + 
    (CA*Ljk*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (CF*Ljk*sqr(ik)*jk)/
    ((ij + ik)*sqr(ik + jk)) - 
    (CA*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CF*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (16*CA*Box6ijk*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (32*CF*Box6ijk*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CA*Lij*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (4*CF*Lij*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CA*Ljk*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (4*CF*Ljk*ij*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (16*CA*Box6ijk*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*sqr(ik + jk)) + 
    (32*CF*Box6ijk*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*sqr(ik + jk)) + 
    (CA*Lij*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*sqr(ik + jk)) - 
    (2*CF*Lij*sqr(ij)*sqr(jk))/
    (ik*(ij + ik)*sqr(ik + jk)) - 
    (CA*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CF*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (CA*Lij*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (2*CF*Lij*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (2*CA*Ljk*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) - 
    (2*CF*Ljk*ik*sqr(jk))/
    ((ij + ik)*sqr(ik + jk)) + 
    (CA*Ljk*pow(jk,3))/
    ((ij + ik)*sqr(ik + jk)) - 
    (CF*Ljk*pow(jk,3))/
    ((ij + ik)*sqr(ik + jk)) - 
    (8*CA*Box6ijk*sqr(ij)*pow(jk,3))/
    (sqr(ik)*(ij + ik)*sqr(ik + jk)) + 
    (16*CF*Box6ijk*sqr(ij)*pow(jk,3))/
    (sqr(ik)*(ij + ik)*sqr(ik + jk)) - 
    (8*CA*Box6ijk*ij*pow(jk,3))/
    (ik*(ij + ik)*sqr(ik + jk)) + 
    (16*CF*Box6ijk*ij*pow(jk,3))/
    (ik*(ij + ik)*sqr(ik + jk)) + 
    (CA*Ljk*ij*pow(jk,3))/
    (ik*(ij + ik)*sqr(ik + jk)) - 
    (2*CF*Ljk*ij*pow(jk,3))/
    (ik*(ij + ik)*sqr(ik + jk));

  qqbargLoops[11] = 
    (16*CA*Box6jik*sqr(ij)*ik)/
    ((ij + jk)*sqr(ik + jk)) - 
    (32*CF*Box6jik*sqr(ij)*ik)/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lij*sqr(ij)*ik)/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lij*sqr(ij)*ik)/
    ((ij + jk)*sqr(ik + jk)) + 
    (8*CA*Box6jik*ij*sqr(ik))/
    ((ij + jk)*sqr(ik + jk)) - 
    (16*CF*Box6jik*ij*sqr(ik))/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lik*ij*sqr(ik))/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lik*ij*sqr(ik))/
    ((ij + jk)*sqr(ik + jk)) + 
    (8*CA*Box6jik*sqr(ij)*sqr(ik))/
    (jk*(ij + jk)*sqr(ik + jk)) - 
    (16*CF*Box6jik*sqr(ij)*sqr(ik))/
    (jk*(ij + jk)*sqr(ik + jk)) + 
    (8*CA*Box6jik*sqr(ij)*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (16*CF*Box6jik*sqr(ij)*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CA*Lij*sqr(ij)*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (4*CF*Lij*sqr(ij)*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (CA*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CF*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (16*CA*Box6jik*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (32*CF*Box6jik*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CA*Lij*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (4*CF*Lij*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CA*Lik*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (4*CF*Lik*ij*ik*jk)/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lik*sqr(ik)*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (CF*Lik*sqr(ik)*jk)/
    ((ij + jk)*sqr(ik + jk)) + 
    (CA*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CF*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (8*CA*Box6jik*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (16*CF*Box6jik*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (3*CA*Lij*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (6*CF*Lij*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lik*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lik*ij*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (CA*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CF*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lij*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lij*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) - 
    (2*CA*Lik*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lik*ik*sqr(jk))/
    ((ij + jk)*sqr(ik + jk)) + 
    (CA*pow(jk,3))/((ij + jk)*sqr(ik + jk)) - 
    (2*CF*pow(jk,3))/((ij + jk)*sqr(ik + jk)) - 
    (CA*Lij*pow(jk,3))/
    ((ij + jk)*sqr(ik + jk)) + 
    (2*CF*Lij*pow(jk,3))/
    ((ij + jk)*sqr(ik + jk)) - 
    (CA*Lik*pow(jk,3))/
    ((ij + jk)*sqr(ik + jk)) + 
    (CF*Lik*pow(jk,3))/((ij + jk)*sqr(ik + jk));

  qqbargLoops[12] = -3*CF*Lijk +
    (CA*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*sqr(ij)*ik)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*ij*sqr(ik))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6jik*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6jik*sqr(ij)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (2*CA*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (18*CF*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (40*CA*Box6ikj*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6jik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6jik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (3*CF*Lik*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (3*CF*Ljk*ij*ik*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6ijk*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6ijk*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (32*CA*Box6ikj*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (3*CF*Lik*sqr(ik)*jk)/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (8*CA*Box6ikj*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6jik*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6jik*ij*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (CA*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (9*CF*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (32*CA*Box6ikj*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (8*CA*Box6jik*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (16*CF*Box6jik*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) + 
    (3*CF*Ljk*ik*sqr(jk))/
    ((ij + ik)*(ij + jk)*(ik + jk)) - 
    (24*CA*Box6ikj*sqr(ik)*sqr(jk))/
    (ij*(ij + ik)*(ij + jk)*
     (ik + jk));

  /* // idendities implied by gauge invariance and current conservation; checked analytically and numerically
  Complex c1 = qqbargLoops[0] + qqbargLoops[6] + qqbargLoops[7];
  Complex c2 = qqbargLoops[1] + qqbargLoops[8] + qqbargLoops[9];
  Complex c3 = qqbargLoops[3] + qqbargLoops[4];
  Complex c4 = qqbargLoops[2] + qqbargLoops[5] + qqbargLoops[10] + qqbargLoops[11];
  Complex c5 = 
    2.*qqbargLoops[3]/ik +
    2.*qqbargLoops[5]/jk +
    qqbargLoops[6]*(1.+ij/ik) +
    qqbargLoops[8]*(jk+ij)/ik +
    2.*qqbargLoops[10]*(1./ik+1./jk) +
    2.*qqbargLoops[12]*(1./ik+1./jk);
  Complex c6 = 
    2.*qqbargLoops[4]/jk +
    2.*qqbargLoops[5]/jk +
    qqbargLoops[7]*(ik+ij)/jk +
    qqbargLoops[9]*(1.+ij/jk) +
    2.*qqbargLoops[11]*(ik/sqr(jk)+1./jk);
  Complex c7 =
    0.5*qqbargLoops[0]*(ij+ik) +
    0.5*qqbargLoops[1]*(ij+jk) +
    qqbargLoops[2]*(1.+ik/jk) -
    qqbargLoops[12]*(1.+ik/jk);

  double x1 = c1 != 0. ? log(abs(real(c1*conj(c1)))) : 0.;
  double x2 = c2 != 0. ? log(abs(real(c2*conj(c2)))) : 0.;
  double x3 = c3 != 0. ? log(abs(real(c3*conj(c3)))) : 0.;
  double x4 = c4 != 0. ? log(abs(real(c4*conj(c4)))) : 0.;
  double x5 = c5 != 0. ? log(abs(real(c5*conj(c5)))) : 0.;
  double x6 = c6 != 0. ? log(abs(real(c6*conj(c6)))) : 0.;
  double x7 = c7 != 0. ? log(abs(real(c7*conj(c7)))) : 0.;

  cerr << x1 << " " << x2 << " " << x3 << " " << x4 << " "
       << x5 << " " << x6 << " " << x7 << "\n";
  */

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralLeftLoopCurrent(int i, int,
								      int j, int,
								      int k, int gHel,
								      int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      (sqrt(2)*c6*plusProduct(j,k)*minusCurrent(n,k)*minusProduct(i,k))/(jk*minusProduct(k,n)) + 
      (sqrt(2)*c1*plusProduct(j,k)*momentum(i)*minusProduct(i,n))/minusProduct(k,n) + 
      (sqrt(2)*c2*plusProduct(j,k)*momentum(j)*minusProduct(i,n))/minusProduct(k,n) + 
      (2*sqrt(2)*c3*plusProduct(j,k)*momentum(k)*minusProduct(i,n))/(jk*minusProduct(k,n)) + 
      (sqrt(2)*c4*plusProduct(i,k)*minusCurrent(i,j)*minusProduct(i,n))/(ik*minusProduct(k,n)) - 
      (sqrt(2)*c7*plusProduct(i,k)*plusProduct(j,k)*momentum(i)*minusProduct(i,k)*minusProduct(i,n))/(ik*minusProduct(k,n)) - 
      (sqrt(2)*c9*plusProduct(i,k)*plusProduct(j,k)*momentum(j)*minusProduct(i,k)*minusProduct(i,n))/(ik*minusProduct(k,n)) - 
      (2*sqrt(2)*c11*plusProduct(i,k)*plusProduct(j,k)*momentum(k)*minusProduct(i,k)*minusProduct(i,n))/(ik*jk*minusProduct(k,n)) + 
      (sqrt(2)*c5*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (sqrt(2)*c8*sqr(plusProduct(j,k))*momentum(i)*minusProduct(i,k)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (sqrt(2)*c10*sqr(plusProduct(j,k))*momentum(j)*minusProduct(i,k)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (2*sqrt(2)*c12*sqr(plusProduct(j,k))*momentum(k)*minusProduct(i,k)*minusProduct(j,n))/(sqr(jk)*minusProduct(k,n));
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c1*plusProduct(j,n)*momentum(i)*minusProduct(i,k))/plusProduct(k,n)) - 
      (sqrt(2)*c2*plusProduct(j,n)*momentum(j)*minusProduct(i,k))/plusProduct(k,n) - 
      (2*sqrt(2)*c3*plusProduct(j,n)*momentum(k)*minusProduct(i,k))/(jk*plusProduct(k,n)) - 
      (sqrt(2)*c4*plusProduct(i,n)*minusCurrent(i,j)*minusProduct(i,k))/(ik*plusProduct(k,n)) + 
      (sqrt(2)*c13*minusCurrent(k,j)*minusProduct(i,k))/ik + (sqrt(2)*c13*minusCurrent(k,j)*minusProduct(i,k))/jk - 
      (sqrt(2)*c6*plusProduct(j,k)*minusCurrent(k,n)*minusProduct(i,k))/(jk*plusProduct(k,n)) + 
      (sqrt(2)*c7*plusProduct(i,n)*plusProduct(j,k)*momentum(i)*sqr(minusProduct(i,k)))/(ik*plusProduct(k,n)) + 
      (sqrt(2)*c9*plusProduct(i,n)*plusProduct(j,k)*momentum(j)*sqr(minusProduct(i,k)))/(ik*plusProduct(k,n)) + 
      (2*sqrt(2)*c11*plusProduct(i,n)*plusProduct(j,k)*momentum(k)*sqr(minusProduct(i,k)))/(ik*jk*plusProduct(k,n)) - 
      (sqrt(2)*c5*plusProduct(j,n)*minusCurrent(i,j)*minusProduct(j,k))/(jk*plusProduct(k,n)) + 
      (sqrt(2)*c8*plusProduct(j,k)*plusProduct(j,n)*momentum(i)*minusProduct(i,k)*minusProduct(j,k))/(jk*plusProduct(k,n)) + 
      (sqrt(2)*c10*plusProduct(j,k)*plusProduct(j,n)*momentum(j)*minusProduct(i,k)*minusProduct(j,k))/(jk*plusProduct(k,n)) + 
      (2*sqrt(2)*c12*plusProduct(j,k)*plusProduct(j,n)*momentum(k)*minusProduct(i,k)*minusProduct(j,k))/(sqr(jk)*plusProduct(k,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedLeftLoopCurrent(int i, int,
								    int j, int,
								    int k, int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c6*plusProduct(j,k)*minusCurrent(i,k))/jk) 
      - (sqrt(2)*c8*sqr(plusProduct(j,k))*momentum(i)*minusProduct(i,j))/jk - 
      (sqrt(2)*c10*sqr(plusProduct(j,k))*momentum(j)*minusProduct(i,j))/jk - 
      (2*sqrt(2)*c12*sqr(plusProduct(j,k))*momentum(k)*minusProduct(i,j))/sqr(jk) + 
      (sqrt(2)*c5*plusProduct(j,k)*minusCurrent(i,j)*minusProduct(i,j))/(jk*minusProduct(i,k));
  }

  if ( gHel == -1 ) {
    return
      (sqrt(2)*c4*plusProduct(i,j)*minusCurrent(i,j)*minusProduct(i,k))/(ik*plusProduct(j,k)) + 
      (sqrt(2)*c13*minusCurrent(k,j)*minusProduct(i,k))/ik + (sqrt(2)*c13*minusCurrent(k,j)*minusProduct(i,k))/jk + 
      (sqrt(2)*c6*minusCurrent(k,j)*minusProduct(i,k))/jk - (sqrt(2)*c7*plusProduct(i,j)*momentum(i)*
								       sqr(minusProduct(i,k)))/ik - 
      (sqrt(2)*c9*plusProduct(i,j)*momentum(j)*sqr(minusProduct(i,k)))/ik - 
      (2*sqrt(2)*c11*plusProduct(i,j)*momentum(k)*sqr(minusProduct(i,k)))/(ik*jk);
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralRightLoopCurrent(int i,    int,
								       int j,    int,
								       int k,    int gHel,
								       int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c13*plusProduct(i,k)*minusCurrent(j,k))/ik) - (sqrt(2)*c13*plusProduct(i,k)*minusCurrent(j,k))/jk + 
      (sqrt(2)*c4*plusProduct(i,k)*minusCurrent(j,i)*minusProduct(i,n))/(ik*minusProduct(k,n)) + 
      (sqrt(2)*c6*plusProduct(i,k)*minusCurrent(n,k)*minusProduct(j,k))/(jk*minusProduct(k,n)) - 
      (sqrt(2)*c7*sqr(plusProduct(i,k))*momentum(i)*minusProduct(i,n)*minusProduct(j,k))/(ik*minusProduct(k,n)) - 
      (sqrt(2)*c9*sqr(plusProduct(i,k))*momentum(j)*minusProduct(i,n)*minusProduct(j,k))/(ik*minusProduct(k,n)) - 
      (2*sqrt(2)*c11*sqr(plusProduct(i,k))*momentum(k)*minusProduct(i,n)*minusProduct(j,k))/(ik*jk*minusProduct(k,n)) + 
      (sqrt(2)*c1*plusProduct(i,k)*momentum(i)*minusProduct(j,n))/minusProduct(k,n) + 
      (sqrt(2)*c2*plusProduct(i,k)*momentum(j)*minusProduct(j,n))/minusProduct(k,n) + 
      (2*sqrt(2)*c3*plusProduct(i,k)*momentum(k)*minusProduct(j,n))/(jk*minusProduct(k,n)) + 
      (sqrt(2)*c5*plusProduct(j,k)*minusCurrent(j,i)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (sqrt(2)*c8*plusProduct(i,k)*plusProduct(j,k)*momentum(i)*minusProduct(j,k)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (sqrt(2)*c10*plusProduct(i,k)*plusProduct(j,k)*momentum(j)*minusProduct(j,k)*minusProduct(j,n))/(jk*minusProduct(k,n)) - 
      (2*sqrt(2)*c12*plusProduct(i,k)*plusProduct(j,k)*momentum(k)*minusProduct(j,k)*minusProduct(j,n))/(sqr(jk)*minusProduct(k,n));
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c4*plusProduct(i,n)*minusCurrent(j,i)*minusProduct(i,k))/(ik*plusProduct(k,n))) - 
      (sqrt(2)*c1*plusProduct(i,n)*momentum(i)*minusProduct(j,k))/plusProduct(k,n) - 
      (sqrt(2)*c2*plusProduct(i,n)*momentum(j)*minusProduct(j,k))/plusProduct(k,n) - 
      (2*sqrt(2)*c3*plusProduct(i,n)*momentum(k)*minusProduct(j,k))/(jk*plusProduct(k,n)) - 
      (sqrt(2)*c5*plusProduct(j,n)*minusCurrent(j,i)*minusProduct(j,k))/(jk*plusProduct(k,n)) - 
      (sqrt(2)*c6*plusProduct(i,k)*minusCurrent(k,n)*minusProduct(j,k))/(jk*plusProduct(k,n)) + 
      (sqrt(2)*c7*plusProduct(i,k)*plusProduct(i,n)*momentum(i)*minusProduct(i,k)*minusProduct(j,k))/(ik*plusProduct(k,n)) + 
      (sqrt(2)*c9*plusProduct(i,k)*plusProduct(i,n)*momentum(j)*minusProduct(i,k)*minusProduct(j,k))/(ik*plusProduct(k,n)) + 
      (2*sqrt(2)*c11*plusProduct(i,k)*plusProduct(i,n)*momentum(k)*minusProduct(i,k)*minusProduct(j,k))/
      (ik*jk*plusProduct(k,n)) + 
      (sqrt(2)*c8*plusProduct(i,k)*plusProduct(j,n)*momentum(i)*sqr(minusProduct(j,k)))/(jk*plusProduct(k,n)) + 
      (sqrt(2)*c10*plusProduct(i,k)*plusProduct(j,n)*momentum(j)*sqr(minusProduct(j,k)))/(jk*plusProduct(k,n)) + 
      (2*sqrt(2)*c12*plusProduct(i,k)*plusProduct(j,n)*momentum(k)*sqr(minusProduct(j,k)))/(sqr(jk)*plusProduct(k,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedRightLoopCurrent(int i, int,
								     int j, int,
								     int k, int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c13*plusProduct(i,k)*minusCurrent(j,k))/ik) - 
      (sqrt(2)*c13*plusProduct(i,k)*minusCurrent(j,k))/jk - 
      (sqrt(2)*c6*plusProduct(i,k)*minusCurrent(j,k))/jk + 
      (sqrt(2)*c7*sqr(plusProduct(i,k))*momentum(i)*minusProduct(i,j))/ik + 
      (sqrt(2)*c9*sqr(plusProduct(i,k))*momentum(j)*minusProduct(i,j))/ik + 
      (2*sqrt(2)*c11*sqr(plusProduct(i,k))*momentum(k)*minusProduct(i,j))/(ik*jk) - 
      (sqrt(2)*c4*plusProduct(i,k)*minusCurrent(j,i)*minusProduct(i,j))/(ik*minusProduct(j,k));
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c5*plusProduct(i,j)*minusCurrent(j,i)*minusProduct(j,k))/(jk*plusProduct(i,k))) + 
      (sqrt(2)*c6*minusCurrent(k,i)*minusProduct(j,k))/jk + 
      (sqrt(2)*c8*plusProduct(i,j)*momentum(i)*sqr(minusProduct(j,k)))/jk + 
      (sqrt(2)*c10*plusProduct(i,j)*momentum(j)*sqr(minusProduct(j,k)))/jk + 
      (2*sqrt(2)*c12*plusProduct(i,j)*momentum(k)*sqr(minusProduct(j,k)))/sqr(jk);
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftOneLoopCurrent(int q,    int qHel,
									 int qbar, int qbarHel,
									 int g1,   int g1Hel) {
  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<2>(1,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5)*qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5)*qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5)*qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5)*qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5)*qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t()*conj(x1.t())) + real(x1.x()*conj(x1.x())) + real(x1.y()*conj(x1.y())) + real(x1.z()*conj(x1.z()));
    double c2 = 
      real(x2.t()*conj(x2.t())) + real(x2.x()*conj(x2.x())) + real(x2.y()*conj(x2.y())) + real(x2.z()*conj(x2.z()));
    double c3 = 
      real(x3.t()*conj(x3.t())) + real(x3.x()*conj(x3.x())) + real(x3.y()*conj(x3.y())) + real(x3.z()*conj(x3.z()));
    double c4 = 
      real(x4.t()*conj(x4.t())) + real(x4.x()*conj(x4.x())) + real(x4.y()*conj(x4.y())) + real(x4.z()*conj(x4.z()));
    ostream& ncheck = checkStream("qqbargLeftLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5)*qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightOneLoopCurrent(int q,    int qHel,
									  int qbar, int qbarHel,
									  int g1,   int g1Hel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<2>(2,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5)*qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5)*qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5)*qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5)*qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5)*qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t()*conj(x1.t())) + real(x1.x()*conj(x1.x())) + real(x1.y()*conj(x1.y())) + real(x1.z()*conj(x1.z()));
    double c2 = 
      real(x2.t()*conj(x2.t())) + real(x2.x()*conj(x2.x())) + real(x2.y()*conj(x2.y())) + real(x2.z()*conj(x2.z()));
    double c3 = 
      real(x3.t()*conj(x3.t())) + real(x3.x()*conj(x3.x())) + real(x3.y()*conj(x3.y())) + real(x3.z()*conj(x3.z()));
    double c4 = 
      real(x4.t()*conj(x4.t())) + real(x4.x()*conj(x4.x())) + real(x4.y()*conj(x4.y())) + real(x4.z()*conj(x4.z()));
    ostream& ncheck = checkStream("qqbargRightLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5)*qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

#ifdef CHECK_MatchboxCurrents

map<string,ofstream*>& MatchboxCurrents::checkStreams() {
  static map<string,ofstream*> theMap;
  return theMap;
}

ostream& MatchboxCurrents::checkStream(const string& id) {
  map<string,ofstream*>::iterator ret = checkStreams().find(id);
  if ( ret == checkStreams().end() ) {
    checkStreams()[id] = new ofstream(id.c_str());
    ret = checkStreams().find(id);
  }
  return *(ret->second);
}

void MatchboxCurrents::checkCurrent(const string& id,
				    const LorentzVector<Complex>& current,
				    const LorentzVector<double>& q) {
  
  Complex c = current.dot(q);
  double ac = abs(real(conj(c)*c));

  if ( isnan(ac) || isinf(ac) ) {
    cerr << "ooops ... nan encountered in current conservation\n" << flush;
    return;
  }

  checkStream(id) << (ac > 0. ? log10(ac) : 0.) << "\n";

}

#endif // CHECK_MatchboxCurrents
