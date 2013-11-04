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
      cacheCurrent(Complex(0.,1.)*plusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,2.)*mass(lbar)/minusProduct(l,lbar))*momentum(l));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,-2.)*mass(l)/plusProduct(l,lbar))*momentum(lbar));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,1.)*mass(l)*mass(lbar)/invariant(l,lbar))*plusCurrent(lbar,l));
  }
  return cachedCurrent();
    
}

const LorentzVector<Complex>& MatchboxCurrents::llbarRightCurrent(int l,    int lHel,
								  int lbar, int lbarHel) {
    
  if ( getCurrent(hash<0>(2,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,1.)*mass(l)*mass(lbar)/invariant(l,lbar))*plusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,-2.)*mass(l)/minusProduct(l,lbar))*momentum(lbar));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,2.)*mass(lbar)/plusProduct(l,lbar))*momentum(l));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent(Complex(0.,1.)*plusCurrent(lbar,l));
  }
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftCurrent(int q,    int qHel,
								 int qbar, int qbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<1>(1,1,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(Complex(0.,1.)*plusCurrent(q,qbar));
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
    cacheCurrent(Complex(0.,1.)*plusCurrent(qbar,q));
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

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*plusCurrent(q,qbar)
		    +(1./plusProduct(g,qbar))*plusCurrent(q,g)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*plusCurrent(q,qbar)
		    +(1./minusProduct(q,g))*plusCurrent(g,qbar)));
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

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*plusCurrent(qbar,q)
		    +(1./plusProduct(q,g))*plusCurrent(qbar,g)));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif
    return cachedCurrent();
  }

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*plusCurrent(qbar,q)
		    +(1./minusProduct(g,qbar))*plusCurrent(g,q)));
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
  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(i,n))/
      (invariant(i,k)*invariant(j,l)*plusProduct(k,n)) + 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(k,j)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(l,j)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,k)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,l)*plusCurrent(i,j)*sqr(plusProduct(i,n)))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(k,j)*sqr(plusProduct(i,n)))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,n)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,j)*sqr(plusProduct(j,n)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(l,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*minusProduct(j,k)*minusProduct(j,n)*plusCurrent(i,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(k,j)*plusProduct(i,l)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(l,j)*plusProduct(i,l)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(i,n)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(l,j)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(k,j)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,j)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(k,n)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,n)*plusCurrent(i,l)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(i,k)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(k,n)*plusCurrent(i,k)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(i,l)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*minusProduct(i,n)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(i,k))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,l)*plusProduct(i,k))/(invariant(i,k)*invariant(j,l)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(k,j)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(j,n)*plusCurrent(i,l)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,l)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(k,j)*plusProduct(i,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(i,k)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(l,n)*plusCurrent(i,l)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n));
  }

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*sqr(minusProduct(i,n))*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(l,n)) + 
      (Complex(0,2)*sqr(minusProduct(j,n))*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(k,n)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(k,n)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*plusCurrent(k,j)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(l,n)*plusCurrent(l,j)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*sqr(minusProduct(j,n))*plusCurrent(i,l)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(i,k)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(i,l)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedLeftCurrent(int i, int,
								 int j, int,
								 int k, int g1Hel,
								 int l, int g2Hel) {

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,j))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(i,j))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,k)*plusProduct(i,j))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(i,j)*sqr(plusProduct(i,j)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(i,k)*plusProduct(i,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(i,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(i,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-1)*sqr(minusProduct(i,k))*plusCurrent(i,j)*plusProduct(i,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(l,j)*plusProduct(i,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,k)*plusCurrent(i,j)*sqr(plusProduct(i,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(k,j)*sqr(plusProduct(i,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,j)*sqr(plusProduct(i,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,k)*plusCurrent(i,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,j)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,k)*plusCurrent(i,l)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(i,k)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(i,k)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(i,l)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,1)*sqr(minusProduct(i,l))*plusCurrent(i,j)*plusProduct(i,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,l)*plusProduct(i,k))/(invariant(i,k)*invariant(j,l)) + 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,j)*plusProduct(i,j)*plusProduct(i,k))/
      (invariant(i,k)*invariant(j,l)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,l)*plusCurrent(i,j)*sqr(plusProduct(i,k)))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*plusCurrent(i,j)*sqr(plusProduct(i,k)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,j)*sqr(plusProduct(i,k)))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,j)*sqr(plusProduct(i,k)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(l,j)*sqr(plusProduct(i,k)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,j)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,l)*plusCurrent(i,l)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(i,k)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,k)*plusCurrent(i,l)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(i,l)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l));
  }

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*minusProduct(i,j)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(i,l)) + 
      (Complex(0,2)*sqr(minusProduct(i,j))*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,2)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*plusCurrent(k,j)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,l)*plusCurrent(l,j)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(i,j)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(i,j)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*sqr(minusProduct(i,j))*plusCurrent(i,l)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*minusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(i,k)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(i,l)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralRightCurrent(int i, int,
								    int j, int,
								    int k, int g1Hel,
								    int l, int g2Hel,
								    int n) {

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(j,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,l)*plusCurrent(j,l)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(k,i)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,l)*plusCurrent(j,k)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(j,k)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*plusProduct(l,n)) + 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(l,i)*plusProduct(j,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,l)*plusCurrent(j,i)*sqr(plusProduct(i,n)))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,n)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(j,i)*sqr(plusProduct(j,n)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(l,i)*sqr(plusProduct(j,n)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,k)*plusProduct(k,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(k,n));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,n)*plusCurrent(j,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(l,k)*plusProduct(j,l))/(invariant(i,k)*invariant(j,l)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,n)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(i,n)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(l,i)*plusProduct(i,n)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(i,n)*plusCurrent(j,k)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(l,i)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(k,n)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(k,n)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(j,k)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,n)*plusCurrent(j,k)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,n)*plusCurrent(j,l)*plusProduct(k,l)*plusProduct(l,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)*plusProduct(k,n));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2)*minusProduct(i,l)*minusProduct(i,n)*plusCurrent(j,l)*plusProduct(i,k))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(i,n)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(i,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,n))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(j,n))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(i,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(k,i)*plusProduct(j,n)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(i,n)*plusCurrent(j,k)*plusProduct(i,k)*plusProduct(k,n))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,n)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(j,k)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(j,l)*minusProduct(l,n)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(k,l)*minusProduct(l,n)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(l,n)*plusCurrent(j,k)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(k,l)*plusCurrent(j,l)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(l,n)*plusCurrent(j,l)*plusProduct(k,l)*plusProduct(k,n))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*plusProduct(l,n));
  }

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*sqr(minusProduct(i,n))*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,n)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(l,i)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(k,n)) + 
      (Complex(0,2)*sqr(minusProduct(j,n))*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n)) + 
      (Complex(0,2)*sqr(minusProduct(i,n))*plusCurrent(j,k)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)*minusProduct(l,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*minusProduct(k,n)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(l,n)) - 
      (Complex(0,2)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(l,n)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(l,n)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,n)) + 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(j,k)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,n)) - 
      (Complex(0,2)*minusProduct(i,n)*plusCurrent(j,l)*sqr(plusProduct(k,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(l,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedRightCurrent(int i, int,
								  int j, int,
								  int k, int g1Hel,
								  int l, int g2Hel) {

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(j,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(k,i)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(j,k)*plusProduct(i,j))/
      (invariant(i,k)*invariant(j,l)*plusProduct(i,l)) + 
      (Complex(0,2)*sqr(minusProduct(k,l))*plusCurrent(l,i)*plusProduct(i,j))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(j,k)*minusProduct(j,l)*plusCurrent(j,i)*sqr(plusProduct(i,j)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(j,l)*minusProduct(k,l)*plusCurrent(l,i)*sqr(plusProduct(i,j)))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,k)*plusProduct(i,k))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(j,k)*plusProduct(i,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(j,l)*plusProduct(i,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,k));
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2)*sqr(minusProduct(i,k))*plusCurrent(j,i)*plusProduct(i,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,1)*sqr(minusProduct(i,k))*plusCurrent(j,i)*plusProduct(i,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,k)*plusCurrent(j,i)*sqr(plusProduct(i,l)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,l)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(i,k)*plusCurrent(j,k)*plusProduct(j,l))/
      (invariant(i,k)*invariant(j,l)*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(l,k)*plusProduct(j,l))/(invariant(i,k)*invariant(j,l)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,j)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,j)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(k,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,l)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,2)*sqr(minusProduct(i,k))*plusCurrent(j,k)*plusProduct(k,l))/
      (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,2)*sqr(minusProduct(i,k))*plusCurrent(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,k)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k)) + 
      (Complex(0,1)*minusProduct(i,k)*plusCurrent(j,k)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,k)) - 
      (Complex(0,1)*sqr(minusProduct(i,k))*plusCurrent(j,l)*plusProduct(i,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,l)*plusProduct(i,k));
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,1)*sqr(minusProduct(i,l))*plusCurrent(j,i)*plusProduct(i,k))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,1)*minusProduct(i,l)*plusCurrent(j,i)*sqr(plusProduct(i,k)))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(j,l)*plusCurrent(k,i)*plusProduct(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(l,i)*plusProduct(i,j)*plusProduct(j,k))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,j)*minusProduct(k,l)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(k,l)*plusCurrent(k,i)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(l,i)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,1)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(k,i)*plusProduct(i,k)*plusProduct(j,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,2)*sqr(minusProduct(i,l))*plusCurrent(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(j,l)*plusCurrent(j,i)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,l)*minusProduct(k,l)*plusCurrent(k,i)*plusProduct(i,j)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) - 
      (Complex(0,1)*sqr(minusProduct(i,l))*plusCurrent(j,k)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(i,k)*plusProduct(i,l)) + 
      (Complex(0,1)*minusProduct(i,l)*plusCurrent(j,l)*plusProduct(i,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(i,l));
  }

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2)*sqr(minusProduct(i,j))*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) + 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(j,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(j,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,2)*minusProduct(i,k)*plusCurrent(k,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,l)) - 
      (Complex(0,2)*plusCurrent(l,i)*plusProduct(j,k)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,j)*plusCurrent(j,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
      (Complex(0,2)*minusProduct(i,l)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k)) - 
      (Complex(0,2)*minusProduct(i,l)*plusCurrent(l,i)*plusProduct(j,l)*plusProduct(k,l))/
      (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(i,k));
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

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(k,i)*minusProduct(i,l)*plusCurrent(i,j)+
		     plusProduct(i,k)*minusProduct(l,k)*plusCurrent(k,j))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,k)*minusProduct(l,j)*plusCurrent(i,j)+
		     plusProduct(l,k)*minusProduct(l,j)*plusCurrent(i,l))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(l,i)*minusProduct(i,k)*plusCurrent(i,j)+
		     plusProduct(i,l)*minusProduct(k,l)*plusCurrent(l,j))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,l)*minusProduct(k,j)*plusCurrent(i,j)+
		     plusProduct(k,l)*minusProduct(k,j)*plusCurrent(i,k))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
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

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(k,i)*minusProduct(i,l)*plusCurrent(j,i)+
		     plusProduct(l,k)*minusProduct(l,i)*plusCurrent(j,l))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,k)*minusProduct(l,j)*plusCurrent(j,i)+
		     plusProduct(j,k)*minusProduct(l,k)*plusCurrent(k,i))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(l,i)*minusProduct(i,k)*plusCurrent(j,i)+
		     plusProduct(k,l)*minusProduct(k,i)*plusCurrent(j,k))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,l)*minusProduct(k,j)*plusCurrent(j,i)+
		     plusProduct(j,l)*minusProduct(k,l)*plusCurrent(l,i))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
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

  double sij = invariant(i,j);
  double sik = invariant(i,k);
  double sjk = invariant(j,k);

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

  double Q2 = invariant(i,j) + invariant(i,k) + invariant(j,k);
  // checked for LEP that virtuals + I operator are mu2 independent
  //double xmu2 = 10*GeV2/sqr(amplitudeScale());
  double xmu2 = 1.;

  Complex Lijk = log(1.,-xmu2/Q2);

  Complex Lij = log(1.,Q2,invariant(i,j));
  Complex Lik = log(1.,Q2,invariant(i,k));
  Complex Ljk = log(1.,Q2,invariant(j,k));

  Complex Box6ijk = box6(i,j,k);
  Complex Box6ikj = box6(i,k,j);
  Complex Box6jik = box6(j,i,k);

  // get the coefficients

  qqbargLoops[0] = 
    (2*CF*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (32*CA*Box6ijk*sqr(invariant(i,j)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (64*CF*Box6ijk*sqr(invariant(i,j)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*sqr(invariant(i,j)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*sqr(invariant(i,j)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Lij*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lij*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lik*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Ljk*sqr(invariant(i,j)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6ijk*pow(invariant(i,j),3))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6ijk*pow(invariant(i,j),3))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*pow(invariant(i,j),3))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*pow(invariant(i,j),3))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6ijk*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6ijk*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CA*Lik*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Ljk*invariant(i,j)*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*sqr(invariant(i,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*sqr(invariant(i,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lik*sqr(invariant(i,k)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*pow(invariant(i,j),3))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*pow(invariant(i,j),3))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (40*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (80*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Lij*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lij*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CF*Ljk*invariant(i,j)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*pow(invariant(i,j),3)*invariant(j,k))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*pow(invariant(i,j),3)*invariant(j,k))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (32*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (64*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6ijk*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6ijk*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (48*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CA*Lik*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Lik*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Ljk*invariant(i,k)*invariant(j,k))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(j,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(j,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*sqr(invariant(j,k)))/
    (sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Ljk*sqr(invariant(j,k)))/(sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (48*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (sqr(invariant(i,j))*sqr(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k)));

  qqbargLoops[1] = 
    (-2*CF*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (32*CA*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (64*CF*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (2*CA*Lij*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lij*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lik*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Ljk*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*pow(invariant(i,j),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*pow(invariant(i,j),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (2*CF*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (40*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (80*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (2*CA*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (12*CF*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Ljk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Lik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*pow(invariant(i,j),3)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*pow(invariant(i,j),3)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6jik*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6jik*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (32*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (64*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (2*CF*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Ljk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Ljk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (2*CF*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (48*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (16*CA*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CF*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Lik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Ljk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Ljk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (48*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Ljk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (sqr(invariant(i,j))*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,j) + invariant(j,k)));

  qqbargLoops[2] = -3*CF*Lijk +
    (-4*CA*Box6jik*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*pow(invariant(i,j),3))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lij*pow(invariant(i,j),3))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j))*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lij*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,j))*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*pow(invariant(i,j),3)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (8*CF*Box6jik*pow(invariant(i,j),3)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (4*CA*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (8*CF*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CA*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (12*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (24*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,j))*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*sqr(invariant(i,j))*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (CF*Ljk*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*pow(invariant(i,j),3)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (8*CF*Box6ijk*pow(invariant(i,j),3)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (CA*Lij*pow(invariant(i,j),3)*invariant(j,k))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CF*Lij*pow(invariant(i,j),3)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (2*CA*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (18*CF*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (16*CA*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (32*CF*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (28*CA*Box6ikj*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CF*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (3*CF*Ljk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (20*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (CA*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (12*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (24*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (20*CA*Box6ikj*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (CF*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*sqr(invariant(j,k)))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CA*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CF*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*pow(invariant(i,j),3)*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (8*CF*Box6ijk*pow(invariant(i,j),3)*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (8*CA*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (16*CF*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (CA*Lij*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CF*Lij*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (CA*Ljk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CF*Ljk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CA*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (32*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k)*sqr(invariant(j,k)))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    - (CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (3*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (12*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (12*CA*Box6ikj*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*pow(invariant(j,k),3))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k)))
    + (3*CF*Ljk*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*sqr(invariant(i,j))*pow(invariant(j,k),3))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (8*CF*Box6ijk*sqr(invariant(i,j))*pow(invariant(j,k),3))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (4*CA*Box6ijk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (8*CF*Box6ijk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (CA*Ljk*invariant(i,j)*pow(invariant(j,k),3))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) + (CF*Ljk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k))) - (12*CA*Box6ikj*invariant(i,k)*pow(invariant(j,k),3))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k)));

  qqbargLoops[3] = 3*CF*Lijk +
    (8*CF*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ikj*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lij*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ikj*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6jik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6jik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ikj*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lij*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Ljk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Lij*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Ljk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (2.*invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6ijk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6ijk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(j,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CF*Ljk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (2.*invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (12*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (sqr(invariant(i,j))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k)));

  qqbargLoops[4] = -3*CF*Lijk +
    (-8*CF*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ikj*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*sqr(invariant(i,j)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lij*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ikj*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (12*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (24*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6jik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6jik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (12*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (24*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ikj*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lij*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Ljk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Lij*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (12*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (2.*invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6ijk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6ijk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*sqr(invariant(j,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CF*Ljk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (12*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (2.*invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (12*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (sqr(invariant(i,j))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k)));

  qqbargLoops[5] = 3*CF*Lijk +
    (-4*CA*Box6jik*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lij*sqr(invariant(i,j)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,j)*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (9*CF*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ikj*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*invariant(i,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (9*CF*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ikj*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*sqr(invariant(i,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lik*sqr(invariant(i,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,j)*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (9*CF*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (4*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ikj*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CF*Lij*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lik*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,j)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CF*Ljk*invariant(i,j)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*sqr(invariant(i,j))*invariant(j,k))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CF*Lij*sqr(invariant(i,j))*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,k)*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (9*CF*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (20*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CF*Lik*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*invariant(j,k))/
    (2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (12*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (12*CA*Box6ikj*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*sqr(invariant(j,k)))/(2.*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Ljk*sqr(invariant(j,k)))/((invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CF*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (4*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    (2.*invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (CF*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k))) + 
    (12*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,k) + invariant(j,k)));

  qqbargLoops[6] = 
    (-2*CF*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) + 
    (32*CA*Box6ijk*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) - 
    (64*CF*Box6ijk*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) - 
    (4*CA*Lij*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) + 
    (8*CF*Lij*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) + 
    (4*CF*Ljk*invariant(i,j))/sqr(invariant(i,j) + invariant(i,k)) + 
    (16*CA*Box6ijk*sqr(invariant(i,j)))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) - 
    (32*CF*Box6ijk*sqr(invariant(i,j)))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) - 
    (2*CA*Lij*sqr(invariant(i,j)))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) + 
    (4*CF*Lij*sqr(invariant(i,j)))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) - 
    (2*CF*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (16*CA*Box6ijk*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) - 
    (32*CF*Box6ijk*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) - 
    (2*CA*Lij*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (4*CF*Lij*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (4*CF*Ljk*invariant(i,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (16*CA*Box6ijk*invariant(j,k))/sqr(invariant(i,j) + invariant(i,k)) - 
    (32*CF*Box6ijk*invariant(j,k))/sqr(invariant(i,j) + invariant(i,k)) - 
    (2*CA*Ljk*invariant(j,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (6*CF*Ljk*invariant(j,k))/sqr(invariant(i,j) + invariant(i,k)) + 
    (16*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))) - 
    (32*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (sqr(invariant(i,k))*sqr(invariant(i,j) + invariant(i,k))) + 
    (32*CA*Box6ijk*invariant(i,j)*invariant(j,k))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) - 
    (64*CF*Box6ijk*invariant(i,j)*invariant(j,k))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) - 
    (2*CA*Ljk*invariant(i,j)*invariant(j,k))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k))) + 
    (4*CF*Ljk*invariant(i,j)*invariant(j,k))/(invariant(i,k)*sqr(invariant(i,j) + invariant(i,k)));

  qqbargLoops[7] = 
    (8*CA*Box6jik*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*sqr(invariant(i,j)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lik*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Ljk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lij*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Lij*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Lik*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Lik*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (4*CF*Ljk*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (24*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    (sqr(invariant(i,j))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k)));

  qqbargLoops[8] = 
    (-8*CA*Box6ijk*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,j))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,j)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Ljk*invariant(i,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*sqr(invariant(i,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*sqr(invariant(i,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(j,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*sqr(invariant(i,j)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,j)*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) + 
    (2*CF*Lik*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*invariant(j,k)*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(j,k))/((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*invariant(j,k))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (CA*Ljk*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Ljk*invariant(i,k)*invariant(j,k))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))) + 
    (24*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    (sqr(invariant(i,j))*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k)));

  qqbargLoops[9] = 
    (2*CF*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) - 
    (32*CA*Box6jik*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) + 
    (64*CF*Box6jik*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) + 
    (4*CA*Lij*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) - 
    (8*CF*Lij*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) - 
    (4*CF*Lik*invariant(i,j))/sqr(invariant(i,j) + invariant(j,k)) - 
    (16*CA*Box6jik*invariant(i,k))/sqr(invariant(i,j) + invariant(j,k)) + 
    (32*CF*Box6jik*invariant(i,k))/sqr(invariant(i,j) + invariant(j,k)) + 
    (2*CA*Lik*invariant(i,k))/sqr(invariant(i,j) + invariant(j,k)) - 
    (6*CF*Lik*invariant(i,k))/sqr(invariant(i,j) + invariant(j,k)) - 
    (16*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    (sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    (sqr(invariant(j,k))*sqr(invariant(i,j) + invariant(j,k))) - 
    (16*CA*Box6jik*sqr(invariant(i,j)))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (32*CF*Box6jik*sqr(invariant(i,j)))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Lij*sqr(invariant(i,j)))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lij*sqr(invariant(i,j)))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (32*CA*Box6jik*invariant(i,j)*invariant(i,k))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (64*CF*Box6jik*invariant(i,j)*invariant(i,k))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CA*Lik*invariant(i,j)*invariant(i,k))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) - 
    (4*CF*Lik*invariant(i,j)*invariant(i,k))/(invariant(j,k)*sqr(invariant(i,j) + invariant(j,k))) + 
    (2*CF*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k)) - 
    (16*CA*Box6jik*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k)) + 
    (32*CF*Box6jik*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k)) + 
    (2*CA*Lij*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k)) - 
    (4*CF*Lij*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k)) - 
    (4*CF*Lik*invariant(j,k))/sqr(invariant(i,j) + invariant(j,k));

  qqbargLoops[10] = 
    (-8*CA*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CA*Lij*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (4*CF*Lij*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (3*CA*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (6*CF*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Lij*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Ljk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CF*Ljk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CA*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (32*CF*Box6ijk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CA*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (4*CF*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CA*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (4*CF*Ljk*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CA*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (32*CF*Box6ijk*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Lij*sqr(invariant(i,j))*sqr(invariant(j,k)))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Lij*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Lij*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CA*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Ljk*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CF*Ljk*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*sqr(invariant(i,j))*pow(invariant(j,k),3))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*sqr(invariant(i,j))*pow(invariant(j,k),3))/
    (sqr(invariant(i,k))*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ijk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (16*CF*Box6ijk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*Ljk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*Ljk*invariant(i,j)*pow(invariant(j,k),3))/
    (invariant(i,k)*(invariant(i,j) + invariant(i,k))*sqr(invariant(i,k) + invariant(j,k)));

  qqbargLoops[11] = 
    (16*CA*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (32*CF*Box6jik*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lij*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    (invariant(j,k)*(invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j))*sqr(invariant(i,k)))/
    (invariant(j,k)*(invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CA*Lij*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (4*CF*Lij*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (16*CA*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (32*CF*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CA*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (4*CF*Lij*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CA*Lik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (4*CF*Lik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lik*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (3*CA*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (6*CF*Lij*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lij*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CA*Lik*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lik*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CA*pow(invariant(j,k),3))/((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (2*CF*pow(invariant(j,k),3))/((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lij*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (2*CF*Lij*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) - 
    (CA*Lik*pow(invariant(j,k),3))/
    ((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k))) + 
    (CF*Lik*pow(invariant(j,k),3))/((invariant(i,j) + invariant(j,k))*sqr(invariant(i,k) + invariant(j,k)));

  qqbargLoops[12] = -3*CF*Lijk +
    (CA*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*sqr(invariant(i,j))*invariant(i,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*invariant(i,j)*sqr(invariant(i,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*sqr(invariant(i,j))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (2*CA*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (18*CF*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (40*CA*Box6ikj*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (3*CF*Lik*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (3*CF*Ljk*invariant(i,j)*invariant(i,k)*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6ijk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6ijk*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (32*CA*Box6ikj*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (3*CF*Lik*sqr(invariant(i,k))*invariant(j,k))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (8*CA*Box6ikj*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,j)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (CA*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (9*CF*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (32*CA*Box6ikj*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (8*CA*Box6jik*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (16*CF*Box6jik*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) + 
    (3*CF*Ljk*invariant(i,k)*sqr(invariant(j,k)))/
    ((invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*(invariant(i,k) + invariant(j,k))) - 
    (24*CA*Box6ikj*sqr(invariant(i,k))*sqr(invariant(j,k)))/
    (invariant(i,j)*(invariant(i,j) + invariant(i,k))*(invariant(i,j) + invariant(j,k))*
     (invariant(i,k) + invariant(j,k)));

  /* // idendities implied by gauge invariance and current conservation; checked analytically and numerically
  Complex c1 = qqbargLoops[0] + qqbargLoops[6] + qqbargLoops[7];
  Complex c2 = qqbargLoops[1] + qqbargLoops[8] + qqbargLoops[9];
  Complex c3 = qqbargLoops[3] + qqbargLoops[4];
  Complex c4 = qqbargLoops[2] + qqbargLoops[5] + qqbargLoops[10] + qqbargLoops[11];
  Complex c5 = 
    2.*qqbargLoops[3]/invariant(i,k) +
    2.*qqbargLoops[5]/invariant(j,k) +
    qqbargLoops[6]*(1.+invariant(i,j)/invariant(i,k)) +
    qqbargLoops[8]*(invariant(j,k)+invariant(i,j))/invariant(i,k) +
    2.*qqbargLoops[10]*(1./invariant(i,k)+1./invariant(j,k)) +
    2.*qqbargLoops[12]*(1./invariant(i,k)+1./invariant(j,k));
  Complex c6 = 
    2.*qqbargLoops[4]/invariant(j,k) +
    2.*qqbargLoops[5]/invariant(j,k) +
    qqbargLoops[7]*(invariant(i,k)+invariant(i,j))/invariant(j,k) +
    qqbargLoops[9]*(1.+invariant(i,j)/invariant(j,k)) +
    2.*qqbargLoops[11]*(invariant(i,k)/sqr(invariant(j,k))+1./invariant(j,k));
  Complex c7 =
    0.5*qqbargLoops[0]*(invariant(i,j)+invariant(i,k)) +
    0.5*qqbargLoops[1]*(invariant(i,j)+invariant(j,k)) +
    qqbargLoops[2]*(1.+invariant(i,k)/invariant(j,k)) -
    qqbargLoops[12]*(1.+invariant(i,k)/invariant(j,k));

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

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == -1 ) {
    return
      (sqrt(2)*c6*minusProduct(j,k)*plusCurrent(n,k)*plusProduct(i,k))/(invariant(j,k)*plusProduct(k,n)) + 
      (sqrt(2)*c1*minusProduct(j,k)*momentum(i)*plusProduct(i,n))/plusProduct(k,n) + 
      (sqrt(2)*c2*minusProduct(j,k)*momentum(j)*plusProduct(i,n))/plusProduct(k,n) + 
      (2*sqrt(2)*c3*minusProduct(j,k)*momentum(k)*plusProduct(i,n))/(invariant(j,k)*plusProduct(k,n)) + 
      (sqrt(2)*c4*minusProduct(i,k)*plusCurrent(i,j)*plusProduct(i,n))/(invariant(i,k)*plusProduct(k,n)) - 
      (sqrt(2)*c7*minusProduct(i,k)*minusProduct(j,k)*momentum(i)*plusProduct(i,k)*plusProduct(i,n))/(invariant(i,k)*plusProduct(k,n)) - 
      (sqrt(2)*c9*minusProduct(i,k)*minusProduct(j,k)*momentum(j)*plusProduct(i,k)*plusProduct(i,n))/(invariant(i,k)*plusProduct(k,n)) - 
      (2*sqrt(2)*c11*minusProduct(i,k)*minusProduct(j,k)*momentum(k)*plusProduct(i,k)*plusProduct(i,n))/(invariant(i,k)*invariant(j,k)*plusProduct(k,n)) + 
      (sqrt(2)*c5*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (sqrt(2)*c8*sqr(minusProduct(j,k))*momentum(i)*plusProduct(i,k)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (sqrt(2)*c10*sqr(minusProduct(j,k))*momentum(j)*plusProduct(i,k)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (2*sqrt(2)*c12*sqr(minusProduct(j,k))*momentum(k)*plusProduct(i,k)*plusProduct(j,n))/(sqr(invariant(j,k))*plusProduct(k,n));
  }

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c1*minusProduct(j,n)*momentum(i)*plusProduct(i,k))/minusProduct(k,n)) - 
      (sqrt(2)*c2*minusProduct(j,n)*momentum(j)*plusProduct(i,k))/minusProduct(k,n) - 
      (2*sqrt(2)*c3*minusProduct(j,n)*momentum(k)*plusProduct(i,k))/(invariant(j,k)*minusProduct(k,n)) - 
      (sqrt(2)*c4*minusProduct(i,n)*plusCurrent(i,j)*plusProduct(i,k))/(invariant(i,k)*minusProduct(k,n)) + 
      (sqrt(2)*c13*plusCurrent(k,j)*plusProduct(i,k))/invariant(i,k) + (sqrt(2)*c13*plusCurrent(k,j)*plusProduct(i,k))/invariant(j,k) - 
      (sqrt(2)*c6*minusProduct(j,k)*plusCurrent(k,n)*plusProduct(i,k))/(invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c7*minusProduct(i,n)*minusProduct(j,k)*momentum(i)*sqr(plusProduct(i,k)))/(invariant(i,k)*minusProduct(k,n)) + 
      (sqrt(2)*c9*minusProduct(i,n)*minusProduct(j,k)*momentum(j)*sqr(plusProduct(i,k)))/(invariant(i,k)*minusProduct(k,n)) + 
      (2*sqrt(2)*c11*minusProduct(i,n)*minusProduct(j,k)*momentum(k)*sqr(plusProduct(i,k)))/(invariant(i,k)*invariant(j,k)*minusProduct(k,n)) - 
      (sqrt(2)*c5*minusProduct(j,n)*plusCurrent(i,j)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c8*minusProduct(j,k)*minusProduct(j,n)*momentum(i)*plusProduct(i,k)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c10*minusProduct(j,k)*minusProduct(j,n)*momentum(j)*plusProduct(i,k)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) + 
      (2*sqrt(2)*c12*minusProduct(j,k)*minusProduct(j,n)*momentum(k)*plusProduct(i,k)*plusProduct(j,k))/(sqr(invariant(j,k))*minusProduct(k,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedLeftLoopCurrent(int i, int,
								    int j, int,
								    int k, int gHel) {

  qqbargLoopCoefficients(i,j,k);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c6*minusProduct(j,k)*plusCurrent(i,k))/invariant(j,k)) 
      - (sqrt(2)*c8*sqr(minusProduct(j,k))*momentum(i)*plusProduct(i,j))/invariant(j,k) - 
      (sqrt(2)*c10*sqr(minusProduct(j,k))*momentum(j)*plusProduct(i,j))/invariant(j,k) - 
      (2*sqrt(2)*c12*sqr(minusProduct(j,k))*momentum(k)*plusProduct(i,j))/sqr(invariant(j,k)) + 
      (sqrt(2)*c5*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,j))/(invariant(j,k)*plusProduct(i,k));
  }

  if ( gHel == 1 ) {
    return
      (sqrt(2)*c4*minusProduct(i,j)*plusCurrent(i,j)*plusProduct(i,k))/(invariant(i,k)*minusProduct(j,k)) + 
      (sqrt(2)*c13*plusCurrent(k,j)*plusProduct(i,k))/invariant(i,k) + (sqrt(2)*c13*plusCurrent(k,j)*plusProduct(i,k))/invariant(j,k) + 
      (sqrt(2)*c6*plusCurrent(k,j)*plusProduct(i,k))/invariant(j,k) - (sqrt(2)*c7*minusProduct(i,j)*momentum(i)*
								       sqr(plusProduct(i,k)))/invariant(i,k) - 
      (sqrt(2)*c9*minusProduct(i,j)*momentum(j)*sqr(plusProduct(i,k)))/invariant(i,k) - 
      (2*sqrt(2)*c11*minusProduct(i,j)*momentum(k)*sqr(plusProduct(i,k)))/(invariant(i,k)*invariant(j,k));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralRightLoopCurrent(int i,    int,
								       int j,    int,
								       int k,    int gHel,
								       int n) {

  qqbargLoopCoefficients(i,j,k);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c13*minusProduct(i,k)*plusCurrent(j,k))/invariant(i,k)) - (sqrt(2)*c13*minusProduct(i,k)*plusCurrent(j,k))/invariant(j,k) + 
      (sqrt(2)*c4*minusProduct(i,k)*plusCurrent(j,i)*plusProduct(i,n))/(invariant(i,k)*plusProduct(k,n)) + 
      (sqrt(2)*c6*minusProduct(i,k)*plusCurrent(n,k)*plusProduct(j,k))/(invariant(j,k)*plusProduct(k,n)) - 
      (sqrt(2)*c7*sqr(minusProduct(i,k))*momentum(i)*plusProduct(i,n)*plusProduct(j,k))/(invariant(i,k)*plusProduct(k,n)) - 
      (sqrt(2)*c9*sqr(minusProduct(i,k))*momentum(j)*plusProduct(i,n)*plusProduct(j,k))/(invariant(i,k)*plusProduct(k,n)) - 
      (2*sqrt(2)*c11*sqr(minusProduct(i,k))*momentum(k)*plusProduct(i,n)*plusProduct(j,k))/(invariant(i,k)*invariant(j,k)*plusProduct(k,n)) + 
      (sqrt(2)*c1*minusProduct(i,k)*momentum(i)*plusProduct(j,n))/plusProduct(k,n) + 
      (sqrt(2)*c2*minusProduct(i,k)*momentum(j)*plusProduct(j,n))/plusProduct(k,n) + 
      (2*sqrt(2)*c3*minusProduct(i,k)*momentum(k)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) + 
      (sqrt(2)*c5*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (sqrt(2)*c8*minusProduct(i,k)*minusProduct(j,k)*momentum(i)*plusProduct(j,k)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (sqrt(2)*c10*minusProduct(i,k)*minusProduct(j,k)*momentum(j)*plusProduct(j,k)*plusProduct(j,n))/(invariant(j,k)*plusProduct(k,n)) - 
      (2*sqrt(2)*c12*minusProduct(i,k)*minusProduct(j,k)*momentum(k)*plusProduct(j,k)*plusProduct(j,n))/(sqr(invariant(j,k))*plusProduct(k,n));
  }

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c4*minusProduct(i,n)*plusCurrent(j,i)*plusProduct(i,k))/(invariant(i,k)*minusProduct(k,n))) - 
      (sqrt(2)*c1*minusProduct(i,n)*momentum(i)*plusProduct(j,k))/minusProduct(k,n) - 
      (sqrt(2)*c2*minusProduct(i,n)*momentum(j)*plusProduct(j,k))/minusProduct(k,n) - 
      (2*sqrt(2)*c3*minusProduct(i,n)*momentum(k)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) - 
      (sqrt(2)*c5*minusProduct(j,n)*plusCurrent(j,i)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) - 
      (sqrt(2)*c6*minusProduct(i,k)*plusCurrent(k,n)*plusProduct(j,k))/(invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c7*minusProduct(i,k)*minusProduct(i,n)*momentum(i)*plusProduct(i,k)*plusProduct(j,k))/(invariant(i,k)*minusProduct(k,n)) + 
      (sqrt(2)*c9*minusProduct(i,k)*minusProduct(i,n)*momentum(j)*plusProduct(i,k)*plusProduct(j,k))/(invariant(i,k)*minusProduct(k,n)) + 
      (2*sqrt(2)*c11*minusProduct(i,k)*minusProduct(i,n)*momentum(k)*plusProduct(i,k)*plusProduct(j,k))/
      (invariant(i,k)*invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c8*minusProduct(i,k)*minusProduct(j,n)*momentum(i)*sqr(plusProduct(j,k)))/(invariant(j,k)*minusProduct(k,n)) + 
      (sqrt(2)*c10*minusProduct(i,k)*minusProduct(j,n)*momentum(j)*sqr(plusProduct(j,k)))/(invariant(j,k)*minusProduct(k,n)) + 
      (2*sqrt(2)*c12*minusProduct(i,k)*minusProduct(j,n)*momentum(k)*sqr(plusProduct(j,k)))/(sqr(invariant(j,k))*minusProduct(k,n));
  }

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedRightLoopCurrent(int i, int,
								     int j, int,
								     int k, int gHel) {

  qqbargLoopCoefficients(i,j,k);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == -1 ) {
    return
      -((sqrt(2)*c13*minusProduct(i,k)*plusCurrent(j,k))/invariant(i,k)) - 
      (sqrt(2)*c13*minusProduct(i,k)*plusCurrent(j,k))/invariant(j,k) - 
      (sqrt(2)*c6*minusProduct(i,k)*plusCurrent(j,k))/invariant(j,k) + 
      (sqrt(2)*c7*sqr(minusProduct(i,k))*momentum(i)*plusProduct(i,j))/invariant(i,k) + 
      (sqrt(2)*c9*sqr(minusProduct(i,k))*momentum(j)*plusProduct(i,j))/invariant(i,k) + 
      (2*sqrt(2)*c11*sqr(minusProduct(i,k))*momentum(k)*plusProduct(i,j))/(invariant(i,k)*invariant(j,k)) - 
      (sqrt(2)*c4*minusProduct(i,k)*plusCurrent(j,i)*plusProduct(i,j))/(invariant(i,k)*plusProduct(j,k));
  }

  if ( gHel == 1 ) {
    return
      -((sqrt(2)*c5*minusProduct(i,j)*plusCurrent(j,i)*plusProduct(j,k))/(invariant(j,k)*minusProduct(i,k))) + 
      (sqrt(2)*c6*plusCurrent(k,i)*plusProduct(j,k))/invariant(j,k) + 
      (sqrt(2)*c8*minusProduct(i,j)*momentum(i)*sqr(plusProduct(j,k)))/invariant(j,k) + 
      (sqrt(2)*c10*minusProduct(i,j)*momentum(j)*sqr(plusProduct(j,k)))/invariant(j,k) + 
      (2*sqrt(2)*c12*minusProduct(i,j)*momentum(k)*sqr(plusProduct(j,k)))/sqr(invariant(j,k));
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
