// - * - C++ - * -
//
// MatchboxCurrents.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxCurrents.h"
#include "Herwig/Utilities/Maths.h"

using namespace Herwig;
using namespace Herwig::Math;

using Constants::pi;

namespace {

  const static LorentzVector<Complex> czero(0.,0.,0.,0.);

  inline Complex csqr(const Complex & a) {
    return a * a;
  }

  inline double theta(const double x) {
    if ( x >= 0. )
      return 1.;
    return 0.;
  }
  inline double sign(const double x) {
    if ( x >= 0. )
      return 1.;
    return -1.;
  }

  // quick'n'dirty fix to template troubles

  Complex operator * (const Complex& a, const double b) {
    return Complex(a.real() * b,a.imag() * b);
  }

  Complex operator * (const double b, const Complex& a) {
    return Complex(a.real() * b,a.imag() * b);
  }

  Complex operator+(const Complex& a, const double b) {
    return Complex(a.real()+b,a.imag());
  }

  Complex operator+(const double b, const Complex& a) {
    return Complex(a.real()+b,a.imag());
  }

  Complex operator-(const Complex& a, const double b) {
    return Complex(a.real()-b,a.imag());
  }

  Complex operator-(const double b, const Complex& a) {
    return Complex(b-a.real(),-a.imag());
  }


  // end fix, needs to be looked at in ThePEG/Config/

}


void MatchboxCurrents::setupLeptons(const int l,    const Lorentz5Momentum& pl,
				    const int lbar, const Lorentz5Momentum& plbar) {

  const Energy4 Delta = (sqr(pl*plbar) - (pl*pl)*(plbar*plbar));
  const Energy2 prod = pl*plbar;

  // Variable to contain the sign of pl*plbar
  double sgn;
  if (prod < ZERO ) {sgn = -1;}
  else if (prod > ZERO) {sgn = 1;}
  else {sgn = 0;}

  InvEnergy2 fact = 0.5/(sgn*sqrt(Delta));

  Lorentz5Momentum lmassless    = ( double(fact*(sgn*sqrt(Delta) + prod))*pl    - double(fact*(   pl*pl))*plbar );
  Lorentz5Momentum lbarmassless = ( double(fact*(sgn*sqrt(Delta) + prod))*plbar - double(fact*(plbar*plbar))*pl );

  lmassless.setMass(ZERO); lmassless.rescaleEnergy();
  lbarmassless.setMass(ZERO); lbarmassless.rescaleEnergy();

  if ( pl.t() < ZERO )
    lmassless.setT(-lmassless.t());

  if ( plbar.t() < ZERO )
    lbarmassless.setT(-lbarmassless.t());

  momentum(l,lmassless,true,pl.mass());
  momentum(lbar,lbarmassless,true,plbar.mass());
}


void MatchboxCurrents::setupQuarks(const int q,    const Lorentz5Momentum& pq,
				   const int qbar, const Lorentz5Momentum& pqbar) {

  const Energy4 Delta = (sqr(pq*pqbar) - (pq*pq)*(pqbar*pqbar));
  const Energy2 prod = pq*pqbar;

  // Variable to contain the sign of pq*pqbar
  double sgn;
  if (prod < ZERO) {sgn = -1;}
  else if (prod > ZERO) {sgn = 1;}
  else {sgn = 0;}
  
  InvEnergy2 fact =  0.5/(sgn*sqrt(Delta));

  Lorentz5Momentum qmassless    = ( double(fact*(sgn*sqrt(Delta) + prod))*pq    - double(fact*(pq*pq))*pqbar );
  Lorentz5Momentum qbarmassless = ( double(fact*(sgn*sqrt(Delta) + prod))*pqbar - double(fact*(pqbar*pqbar))*pq );

  qmassless.setMass(ZERO); qmassless.rescaleEnergy();
  qbarmassless.setMass(ZERO); qbarmassless.rescaleEnergy();

  if ( pq.t() < ZERO )
    qmassless.setT(-qmassless.t());

  if ( pqbar.t() < ZERO )
    qbarmassless.setT(-qbarmassless.t());

  momentum(q,qmassless,true,pq.mass());
  momentum(qbar,qbarmassless,true,pqbar.mass());
}


const LorentzVector<Complex>& MatchboxCurrents::llbarLeftCurrent(const int l,    const int lHel,
								 const int lbar, const int lbarHel) {

  if ( getCurrent(hash<0>(1,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,2.) * mass(lbar)/plusProduct(l,lbar)) * momentum(l));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,-2.) * mass(l)/minusProduct(l,lbar)) * momentum(lbar));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,1.) * mass(l) * mass(lbar)/invariant(l,lbar)) * minusCurrent(lbar,l));
  }

  return cachedCurrent();
}

const LorentzVector<Complex>& MatchboxCurrents::llbarRightCurrent(const int l,    const int lHel,
								  const int lbar, const int lbarHel) {
    
  if ( getCurrent(hash<0>(2,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,1.) * mass(l) * mass(lbar)/invariant(l,lbar)) * minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,-2.) * mass(l)/plusProduct(l,lbar)) * momentum(lbar));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,2.) * mass(lbar)/minusProduct(l,lbar)) * momentum(l));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(lbar,l));
  }

  return cachedCurrent();
}


const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftCurrent(const int q,    const int qHel,
								 const int qbar, const int qbarHel) {

  if ( getCurrent(hash<1>(1,1,q,qHel,qbar,qbarHel)) ) {
    if ( qHel == 1 && qbarHel == 1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(q,qbar));
    if ( qHel == 1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,2.) * mass(qbar)/plusProduct(q,qbar)) * momentum(q));
    if ( qHel == -1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,-2.) * mass(q)/minusProduct(q,qbar)) * momentum(qbar));
    if ( qHel == -1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,1.) * mass(q) * mass(qbar)/invariant(q,qbar)) * minusCurrent(qbar,q));
  }

#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  
  return cachedCurrent();
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightCurrent(const int q,    const int qHel,
								  const int qbar, const int qbarHel) {
    
  if ( getCurrent(hash<1>(2,1,q,qHel,qbar,qbarHel)) ) {
    if ( qHel == 1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,1.) * mass(q) * mass(qbar)/invariant(q,qbar)) * minusCurrent(q,qbar));
    if ( qHel == 1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,-2.) * mass(q)/plusProduct(q,qbar)) * momentum(qbar));
    if ( qHel == -1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,2.) * mass(qbar)/minusProduct(q,qbar)) * momentum(q));
    if ( qHel == -1 && qbarHel == -1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(qbar,q));
  }
 
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  
  return cachedCurrent();
}


const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftCurrent(const int q,    const int qHel,
								  const int qbar, const int qbarHel,
								  const int g,    const int gHel) {

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the negative helicity gluon 
      // Note that the gluon is outgoing so the polarisation vector of the hel=+1 gluon is conjugated to give the hel=-1 vector
      const Complex cminus = sqrt(2.0) / minusProduct(g,q);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*minusCurrent(q, qbar) - (minusProduct(g,q)*plusProduct(g,qbar)/den_j)*minusCurrent(q,g) ) ); 
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(-mass(qbar)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*2*momentum(q) + (invariant(q,g)/den_j)*minusCurrent(q,g) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(q)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (plusProduct(g,qbar)*minusProduct(qbar,q)/den_j))*2*momentum(qbar) - (minusProduct(g,q)*plusProduct(g,qbar)/den_j)*minusCurrent(qbar,g) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (minusProduct(q,qbar)*plusProduct(qbar,g)/den_j))*minusCurrent(qbar,q) + (invariant(q,g)/den_j)*minusCurrent(qbar,g) ) );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) { 

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cplus = sqrt(2.0) / plusProduct(q,g);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*minusCurrent(q, qbar) - (invariant(q,g)/den_i)*minusCurrent(g,qbar) ) ); 
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(-mass(qbar)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*2*momentum(q) - (invariant(q,g)/den_i)*minusCurrent(g,q) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(q)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*2*momentum(qbar) + (minusProduct(qbar,g)*plusProduct(q,g)/den_i)*minusCurrent(g,qbar) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*minusCurrent(qbar, q) + (minusProduct(qbar,g)*plusProduct(q,g)/den_i)*minusCurrent(g,q) ) );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  return czero;
}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightCurrent(const int q,    const int qHel,
								   const int qbar, const int qbarHel,
								   const int g,    const int gHel) {

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);
     
      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cminus = sqrt(2.0) / minusProduct(g,q);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*plusCurrent(qbar, q) + (plusProduct(qbar,g)*minusProduct(q,g)/den_i)*plusCurrent(g,q) ) );
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(q)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*2*momentum(qbar) + (plusProduct(qbar,g)*minusProduct(q,g)/den_i)*plusCurrent(g,qbar) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(-mass(qbar)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*2*momentum(q) - (invariant(q,g)/den_i)*plusCurrent(g,q) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*plusCurrent(q, qbar) - (invariant(q,g)/den_i)*plusCurrent(g,qbar) ) ); 
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cplus = sqrt(2.0) / plusProduct(q,g);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (plusProduct(q,qbar)*minusProduct(qbar,g)/den_j))*plusCurrent(qbar, q) + (invariant(q,g)/den_j)*plusCurrent(qbar,g) ) );
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(q)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (minusProduct(g,qbar)*plusProduct(qbar,q)/den_j))*2*momentum(qbar) - (plusProduct(g,q)*minusProduct(g,qbar)/den_j)*plusCurrent(qbar,g) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(-mass(qbar)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*2*momentum(q) + (invariant(q,g)/den_j)*plusCurrent(q,g) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*plusCurrent(q, qbar) - (plusProduct(g,q)*minusProduct(g,qbar)/den_j)*plusCurrent(q,g) ) ); 
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  return czero;
}


LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralLeftCurrent(const int i, const int,
								   const int j, const int,
								   const int k, const int g1Hel,
								   const int l, const int g2Hel,
								   const int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_kn = plusProduct(k,n);
  const Complex plusP_ln = plusProduct(l,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_kn = minusProduct(k,n);
  const Complex minusP_ln = minusProduct(l,n);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kl = minusCurrent(k,l);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_jl * plusP_kl * minusC_ik)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ik)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_il)/
      (kl * (jk + jl + kl)) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ij * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_il * minusP_in)/
      (ik * jl * minusP_kn) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_kj * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_il * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ij * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_lj * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ij * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ik * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_ij * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_kj * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_ij * minusP_in * minusP_jn)/
      (ik * jl * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ij * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ik * minusP_kn)/
      (kl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_ln)/
      (kl * (jk + jl + kl) * minusP_kn);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_jk * plusP_jn * minusC_ik * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_ij * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_kj * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * minusC_lj * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_ij * minusP_in * minusP_jl)/
      (ik * jl * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_ij * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_ij * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_kl * plusP_kn * minusC_lj * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ij * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_ij * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_lj * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ij * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_ij * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_kj * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_kl * minusC_lj * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_ij * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ij * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ij * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_il * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ik * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ik * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_il * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_in * plusP_jl * minusC_il * minusP_ik)/
      (ik * jl * plusP_kn) + 
      (Complex(0,2) * plusP_jl * minusC_kl * minusP_ik)/(ik * jl) - 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_il * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_ij * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_kj * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_jl * minusC_ij * minusP_ik * minusP_jn)/
      (ik * jl * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * minusC_kj * minusP_ik * minusP_jn)/
      (ik * jl * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_ij * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_ij * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_kl * plusP_ln * minusC_kj * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ij * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_il * minusP_jn * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * minusC_ij * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ij * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_ij * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_kl * minusC_kj * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_kl * minusC_kj * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_lj * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_ij * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_kj * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ij * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_ij * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ij * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ik * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_il * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_il * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_in) * minusC_ij * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_kj * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_lj * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_jn * minusC_ij * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_kj * minusP_ik * minusP_jl)/
      (ik * jl * plusP_ln) + 
      (Complex(0,2) * sqr(plusP_jn) * minusC_ij * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_kn * minusC_kj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_kn * minusC_kj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) + 
      (Complex(0,2) * minusC_kj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ln * minusC_lj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ij * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_ij * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * sqr(plusP_jn) * minusC_il * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_ik * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jn * minusC_il * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_ln);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedLeftCurrent(const int i, const int,
								 const int j, const int,
								 const int k, const int g1Hel,
								 const int l, const int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);

  const Complex plusP_kl = plusProduct(k,l);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_kl = minusProduct(k,l);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kl = minusCurrent(k,l);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_jl * plusP_kl * minusC_ik)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ik)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_il)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_il * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ij * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ik * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ij * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ik * minusP_ik)/
      (kl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_il)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_il)/
      (kl * (jk + jl + kl) * minusP_ik);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-1) * sqr(plusP_ik) * minusC_ij * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_lj * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * minusC_ij * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_kj * sqr(minusP_il))/
      (kl * (ik + il + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_kl * minusC_lj * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_ij * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ik * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ij * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ij * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ij * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_il * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ij * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ik * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ik * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_il * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1) * sqr(plusP_il) * minusC_ij * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_kj * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_jl * minusC_kl * minusP_ik)/(ik * jl) + 
      (Complex(0,2) * plusP_jl * minusC_kj * minusP_ij * minusP_ik)/
      (ik * jl * minusP_il) - 
      (Complex(0,2) * plusP_il * minusC_ij * sqr(minusP_ik))/
      (ik * (ik + il + kl) * minusP_il) - 
      (Complex(0,1) * plusP_il * minusC_ij * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_kl * minusC_kj * sqr(minusP_ik))/
      (ik * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_kl * minusC_kj * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_lj * sqr(minusP_ik))/
      (kl * (ik + il + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_ij * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_ij * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ij * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ij * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_il * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ij * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_il * minusP_ij * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ik * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_il * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_il * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_ij * minusC_kj * minusP_ik * minusP_jl)/
      (ik * jl * plusP_il) + 
      (Complex(0,2) * sqr(plusP_ij) * minusC_ij * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_kj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_kj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * minusC_kj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_il * minusC_lj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ij * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ij * minusC_ij * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * sqr(plusP_ij) * minusC_il * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) - 
      (Complex(0,2) * plusP_ij * minusC_ik * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_ij * minusC_il * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_il);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralRightCurrent(const int i, const int,
								    const int j, const int,
								    const int k, const int g1Hel,
								    const int l, const int g2Hel,
								    const int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_kn = plusProduct(k,n);
  const Complex plusP_ln = plusProduct(l,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_kn = minusProduct(k,n);
  const Complex minusP_ln = minusProduct(l,n);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);
  const LorentzVector<Complex> & minusC_lk = minusCurrent(l,k);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jk)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ji * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_jl * minusP_in)/
      (ik * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ji * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_ki * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ji * minusP_in)/
      (ik * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ji * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_jk * minusP_in)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ji * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_jk * minusP_jn)/
      (ik * jl * minusP_ln) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_li * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_ji * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_ji * minusP_in * minusP_jn)/
      (ik * jl * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ji * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_li * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_kn)/
      (kl * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jl * minusP_ln)/
      (kl * (ik + il + kl) * minusP_kn);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_ik * plusP_kn * minusC_ji * minusP_il)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_jk * minusP_jl)/
      (ik * jl * plusP_ln) + 
      (Complex(0,2) * plusP_ik * minusC_lk * minusP_jl)/(ik * jl) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_jk * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_jk * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_ji * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_ji * minusP_in * minusP_jl)/
      (ik * jl * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_ik * minusC_li * minusP_in * minusP_jl)/
      (ik * jl * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_ji * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_li * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_ji * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_jk * minusP_in * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ji * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_kl * plusP_kn * minusC_li * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_ji * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ji * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_ji * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_ji * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_li * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * minusC_ji * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ji * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ji * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_ki * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_jk * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_jk * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_jl * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_il * plusP_in * minusC_jl * minusP_ik)/
      (ik * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_jl * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_ji * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_jl * minusC_ji * minusP_ik * minusP_jn)/
      (ik * jl * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_ji * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * minusC_ki * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_li * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_ji * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ji * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_kl * plusP_ln * minusC_ki * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ji * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_ji * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_jk * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_ji * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ji * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_ji * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_kl * minusC_ki * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_li * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ji * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_ki * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_jk * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_jl * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_jl * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_in) * minusC_ji * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_jn * minusC_ji * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_in * minusC_li * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn) + 
      (Complex(0,2) * sqr(plusP_jn) * minusC_ji * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_jn * minusC_ki * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_jn * minusC_li * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_in * minusC_ji * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * sqr(plusP_in) * minusC_jk * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ji * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_kn * minusC_ki * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * minusC_li * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_ln * minusC_li * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_ln * minusC_li * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_in * minusC_jk * sqr(minusP_kl))/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_in * minusC_jl * sqr(minusP_kl))/
      (kl * (ik + il + kl) * plusP_ln);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedRightCurrent(const int i, const int,
								  const int j, const int,
								  const int k, const int g1Hel,
								  const int l, const int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);

  const Complex plusP_kl = plusProduct(k,l);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_kl = minusProduct(k,l);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);
  const LorentzVector<Complex> & minusC_lk = minusCurrent(l,k);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jk)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (kl * (ik + il + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ji * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_ki * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ji * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_jk * minusP_ij)/
      (ik * jl * minusP_il) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_li * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ji * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_li * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_ik)/
      (ik * (ik + il + kl) * minusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_ik)/
      (kl * (ik + il + kl) * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jl * minusP_il)/
      (kl * (ik + il + kl) * minusP_ik);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * sqr(plusP_ik) * minusC_ji * minusP_il)/
      (ik * (ik + il + kl) * plusP_il) - 
      (Complex(0,1) * sqr(plusP_ik) * minusC_ji * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * minusC_ji * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_ji * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_li * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_ik * minusC_jk * minusP_jl)/
      (ik * jl * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_lk * minusP_jl)/(ik * jl) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ji * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_li * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_jk * minusC_ji * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ji * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ji * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_ki * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * minusP_ik) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * sqr(plusP_ik) * minusC_jk * minusP_kl)/
      (ik * (ik + il + kl) * plusP_il) - 
      (Complex(0,2) * sqr(plusP_ik) * minusC_jk * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ji * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_li * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * minusC_jk * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * minusP_ik) - 
      (Complex(0,1) * sqr(plusP_ik) * minusC_jl * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il * minusP_ik);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1) * sqr(plusP_il) * minusC_ji * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) - 
      (Complex(0,1) * plusP_il * minusC_ji * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_ji * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jl * minusC_ki * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_li * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_ji * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ji * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_kl * minusC_ki * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_li * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ji * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_ki * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * sqr(plusP_il) * minusC_jl * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ji * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ki * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * sqr(plusP_il) * minusC_jk * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * minusC_jl * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * minusP_il);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_ij) * minusC_ji * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) + 
      (Complex(0,2) * plusP_ij * minusC_ki * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ij * minusC_li * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ik * minusC_ki * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * minusC_li * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_il * minusC_li * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_il * minusC_li * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik);
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggLeftCurrent(const int q,    const int qHel,
								   const int qbar, const int qbarHel,
								   const int g1,   const int g1Hel,
								   const int g2,   const int g2Hel) {
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
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
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

const LorentzVector<Complex>& MatchboxCurrents::qqbarggRightCurrent(const int q,    const int qHel,
								    const int qbar, const int qbarHel,
								    const int g1,   const int g1Hel,
								    const int g2,   const int g2Hel) {

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
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
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

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarLeftCurrent(const int q,    const int qHel,
								      const int qbar, const int qbarHel,
								      const int k,    const int kHel,
								      const int kbar, const int kbarHel) {

  if ( qHel != 1 || qbarHel != 1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  const int i = q; const int j = qbar; const int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_kj = plusProduct(k,j);
  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_lj = plusProduct(l,j);
  const Complex plusP_lk = plusProduct(l,k);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_ki = minusProduct(k,i);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_li = minusProduct(l,i);
  const Complex minusP_lk = minusProduct(l,k);

  
  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_ki * plusP_il * minusC_ij+
		     minusP_ik * plusP_lk * minusC_kj)/
		    (kl+il+ik)-
		    (minusP_jk * plusP_lj * minusC_ij+
		     minusP_lk * plusP_lj * minusC_il)/
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
		   ((minusP_li * plusP_ik * minusC_ij+
		     minusP_il * plusP_kl * minusC_lj)/
		    (kl+il+ik)-
		    (minusP_jl * plusP_kj * minusC_ij+
		     minusP_kl * plusP_kj * minusC_ik)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarRightCurrent(const int q,    const int qHel,
								       const int qbar, const int qbarHel,
								       const int k,    const int kHel,
								       const int kbar, const int kbarHel) {

  if ( qHel != -1 || qbarHel != -1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  const int i = q; const int j = qbar; const int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_ki = plusProduct(k,i);
  const Complex plusP_kj = plusProduct(k,j);
  const Complex plusP_kl = plusProduct(k,l);

  const Complex plusP_li = plusProduct(l,i);
  const Complex plusP_lj = plusProduct(l,j);
  const Complex plusP_lk = plusProduct(l,k);

  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_ki = minusProduct(k,i);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_li = minusProduct(l,i);
  const Complex minusP_lk = minusProduct(l,k);

  
  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_ki * plusP_il * minusC_ji+
		     minusP_lk * plusP_li * minusC_jl)/
		    (kl+il+ik)-
		    (minusP_jk * plusP_lj * minusC_ji+
		     minusP_jk * plusP_lk * minusC_ki)/
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
		   ((minusP_li * plusP_ik * minusC_ji+
		     minusP_kl * plusP_ki * minusC_jk)/
		    (kl+il+ik)-
		    (minusP_jl * plusP_kj * minusC_ji+
		     minusP_jl * plusP_kl * minusC_li)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}


// Definition of sqrt to enable calculation of the sqrt of a negative double
inline Complex sqrt1 (double a) {

  if (a > 0.) { return Complex(sqrt(a), 0.) ;}
  else if (a < 0.) { return Complex(0., sqrt(abs(a))) ;}
  else { return Complex(0., 0.); }
}

// Definition of sqrt to enable calculation of the sqrt of Complex arguments
inline Complex sqrt1 (Complex a) {
  
  const double real_part = sqrt(abs(a))*cos(0.5*arg(a));
  const double imag_part = sqrt(abs(a))*sin(0.5*arg(a)); 
  return Complex(real_part, imag_part) ;
}

// Definition of log to enable continuation of the log of a negative double
inline Complex log1 (double a) {

  if (a < 0.) { return Complex(log(abs(a)), Constants::pi) ;}
  else { return Complex(log(a), 0.) ;}
}

// Definition of log to enable continuation of the log of a Complex argument with a negative real part
inline Complex log1 (Complex a) {

  return Complex(log(abs(a)), arg(a)) ;
}


const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftOneLoopCurrent(const int q,    const int qHel,
									const int qbar, const int qbarHel) {

  // Note this cannot currently handle the case of one massive quark and one massless quark
  assert( (mass(q) == 0 && mass(qbar) == 0) || (mass(q) != 0 && mass(qbar) != 0) );

  // Massless quarks
  if ( mass(q) == 0 && mass(qbar) == 0 ) {

    if ( qHel != 1 || qbarHel != 1 )
      return czero;

    const LorentzVector<Complex>& tree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(1,2,q,qHel,qbar,qbarHel)) ) {
      cacheCurrent( 0.5*CF*( -8. - 3.*log1(-1./invariant(q,qbar)) - sqr(log1(-1./invariant(q,qbar))) ) * tree );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarLeftOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }

  // Massive quarks
  else {
    const LorentzVector<Complex>& momQ = momentum(q) + (sqr(mass(q))/invariant(q,qbar))*momentum(qbar);
    const LorentzVector<Complex>& momQbar = momentum(qbar) + (sqr(mass(qbar))/invariant(q,qbar))*momentum(q);

    const Complex s = (momQ+momQbar).dot(momQ+momQbar);
    const Complex inv12 = s - sqr(mass(q)) - sqr(mass(qbar));

    // Coefficient of the left-handed born-level current
    const Complex coeffLeftTree = -1.0*log1(1./sqr(mass(q)))-1.0*log1(1./sqr(mass(qbar)))-4.0 + 0.5*((2.*log1(sqr(mass(q))/sqr(mass(qbar)))*(0.5*inv12+sqr(mass(q))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))-(2.*inv12*Li2(0.5-(0.25*inv12)/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*sqr(mass(qbar)))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*Li2((0.5*(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(1.*inv12*log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*log1((2*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))*log1((-0.5*inv12-1.*sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(1.*inv12*log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(4*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(-0.25*(3+2*log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*sqr(inv12)+sqr(mass(q))*sqr(mass(qbar))-0.5*inv12*(1.+log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*(sqr(mass(q))+sqr(mass(qbar)))))/((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))));

    // Coefficient of the right-handed born-level current
    const Complex coeffRightTree = (2*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*mass(q)*mass(qbar))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)));      
       
    const LorentzVector<Complex>& leftTree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);
    const LorentzVector<Complex>& rightTree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(1,2,q,qHel,qbar,qbarHel)) ) {

      if ( qHel == 1 && qbarHel == 1 ) {
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree) );
      }

      if ( qHel == 1 && qbarHel == -1 ) {
	// Coefficients of the left and right handed products of massive spinors
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffRightProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * mass(q)*mass(qbar)/plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == 1 ){
	// Coefficients of the left and right handed products of massive spinors
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffRightProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * mass(q)*mass(qbar)/minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == -1 ){
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarLeftOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightOneLoopCurrent(const int q,    const int qHel,
									 const int qbar, const int qbarHel) {

  // Note this cannot currently handle the case of one massive quark and one massless quark
  assert( (mass(q) == 0 && mass(qbar) == 0) || (mass(q) != 0 && mass(qbar) != 0) );

  // Massless quarks
  if ( mass(q) == 0 && mass(qbar) ==0 ) {

    if ( qHel != -1 || qbarHel != -1 )
      return czero;

    const LorentzVector<Complex>& tree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(2,2,q,qHel,qbar,qbarHel)) ) {
      cacheCurrent( 0.5*CF*( -8. - 3.*log1(-1./invariant(q,qbar)) - sqr(log1(-1./invariant(q,qbar))) ) * tree );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarRightOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }

  // Massive quarks
  else {
    const LorentzVector<Complex>& momQ = momentum(q) + (sqr(mass(q))/invariant(q,qbar))*momentum(qbar);
    const LorentzVector<Complex>& momQbar = momentum(qbar) + (sqr(mass(qbar))/invariant(q,qbar))*momentum(q);

    const Complex s = (momQ+momQbar).dot(momQ+momQbar);
    const Complex inv12 = s - sqr(mass(q)) - sqr(mass(qbar));

    // Coefficient of the right-handed born-level current
    const Complex coeffRightTree = -1.0*log1(1./sqr(mass(q)))-1.0*log1(1./sqr(mass(qbar)))-4.0 + 0.5*((2.*log1(sqr(mass(q))/sqr(mass(qbar)))*(0.5*inv12+sqr(mass(q))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))-(2.*inv12*Li2(0.5-(0.25*inv12)/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*sqr(mass(qbar)))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*Li2((0.5*(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(1.*inv12*log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*log1((2*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))*log1((-0.5*inv12-1.*sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(1.*inv12*log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(4*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(-0.25*(3+2*log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*sqr(inv12)+sqr(mass(q))*sqr(mass(qbar))-0.5*inv12*(1.+log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*(sqr(mass(q))+sqr(mass(qbar)))))/((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))));

    // Coefficient of the left-handed born-level current
    const Complex coeffLeftTree = (2*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*mass(q)*mass(qbar))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)));      
       
    const LorentzVector<Complex>& leftTree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);
    const LorentzVector<Complex>& rightTree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(2,2,q,qHel,qbar,qbarHel)) ) {

      if ( qHel == 1 && qbarHel == 1 ) {
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

      if ( qHel == 1 && qbarHel == -1 ) {
	// Coefficients of the right and left handed products of massive spinors
	const LorentzVector<Complex>& coeffRightProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * mass(q)*mass(qbar)/plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == 1 ){
	// Coefficients of the right and left handed products of massive spinors
	const LorentzVector<Complex>& coeffRightProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * mass(q)*mass(qbar)/minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == -1 ){
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

    }
  
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarRightOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }
}


// ln(s(a+i0))
inline Complex log(double s, double a) {
  return
    s < 0. ?
	Complex(log(abs(a)),-pi * theta(a)) :
    Complex(log(abs(a)),pi * theta(-a));
}

// ln(s(a+i0)/(b+i0))
inline Complex log(double s, double a, double b) {
  return
    s < 0. ?
	Complex(log(abs(a/b)),-pi * theta(a/b) * sign(b-a)) :
    Complex(log(abs(a/b)),pi * theta(-a/b) * sign(b-a));
}


// Li2(-(a+i0)/(b+i0))
inline Complex Li2(double a, double b) {
  if ( -a/b < 1. )
    return Complex(Herwig::Math::ReLi2(-a/b),0.0);
  return Complex(Herwig::Math::ReLi2(-a/b),-pi * log(-a/b) * sign(b-a));
}

Complex MatchboxCurrents::box6(const int i, const int j, const int k) {

  const double sij = invariant(i,j);
  const double sik = invariant(i,k);
  const double sjk = invariant(j,k);

  return
    -( Li2(sik+sjk,sij) + Li2(sik+sij,sjk) + 0.5 * csqr(log(1.,sij,sjk)) + sqr(pi)/6. )/8.;
      
}

void MatchboxCurrents::qqbargLoopCoefficients(const int i, const int j, const int k) {

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
  const double ij2 = sqr(ij);
  const double ij3 = ij2 * ij;

  const double ik = invariant(i,k);
  const double ik2 = sqr(ik);
  //const double ik3 = ik2 * ik;

  const double jk = invariant(j,k);
  const double jk2 = sqr(jk);
  const double jk3 = jk2 * jk;

  const double ij_ik = ij + ik;
  const double ij_ik_2 = sqr(ij_ik);

  const double ij_jk = ij + jk;
  const double ij_jk_2 = sqr(ij_jk);

  const double ik_jk = ik + jk;
  const double ik_jk_2 = sqr(ik_jk);


  const double Q2 = ij + ik + jk;
  // checked for LEP that virtuals + I operator are mu2 independent
  //double xmu2 = 10 * GeV2/sqr(amplitudeScale());
  const double xmu2 = 1.;

  const Complex Lijk = log(1.,-xmu2/Q2);

  const Complex Lij = log(1.,Q2,ij);
  const Complex Lik = log(1.,Q2,ik);
  const Complex Ljk = log(1.,Q2,jk);

  const Complex Box6ijk = box6(i,j,k);
  const Complex Box6ikj = box6(i,k,j);
  const Complex Box6jik = box6(j,i,k);

  // get the coefficients
  qqbargLoops[0] = (
		    (2 * CF * ij2) 			- 
		    (32 * CA * Box6ijk * ij2)		+ 
		    (64 * CF * Box6ijk * ij2)		- 
		    (8 * CA * Box6jik * ij2)		+ 
		    (16 * CF * Box6jik * ij2)		+ 
		    (2 * CA * Lij * ij2)		- 
		    (4 * CF * Lij * ij2)		- 
		    (CA * Lik * ij2)			- 
		    (2 * CF * Lik * ij2)		- 
		    (4 * CF * Ljk * ij2)		- 
		    (16 * CA * Box6ijk * ij3)		/ ik + 
		    (32 * CF * Box6ijk * ij3)		/ ik + 
		    (CA * Lij * ij3)			/ ik - 
		    (2 * CF * Lij * ij3)		/ ik + 
		    (2 * CF * ij * ik)			- 
		    (16 * CA * Box6ijk * ij * ik)	+ 
		    (32 * CF * Box6ijk * ij * ik)	- 
		    (16 * CA * Box6jik * ij * ik)	+ 
		    (32 * CF * Box6jik * ij * ik)	+ 
		    (CA * Lij * ij * ik)		- 
		    (2 * CF * Lij * ij * ik)		- 
		    (2 * CA * Lik * ij * ik)		- 
		    (4 * CF * Lik * ij * ik)		- 
		    (4 * CF * Ljk * ij * ik)		- 
		    (8 * CA * Box6jik * ik2)		+ 
		    (16 * CF * Box6jik * ik2)		- 
		    (CA * Lik * ik2)			- 
		    (2 * CF * Lik * ik2)		- 
		    (8 * CA * Box6jik * ij3)		/ jk + 
		    (16 * CF * Box6jik * ij3)		/ jk - 
		    (16 * CA * Box6jik * ij2 * ik)	/ jk + 
		    (32 * CF * Box6jik * ij2 * ik)	/ jk - 
		    (8 * CA * Box6jik * ij * ik2)	/ jk + 
		    (16 * CF * Box6jik * ij * ik2)	/ jk + 
		    (2 * CF * ij * jk)			- 
		    (40 * CA * Box6ijk * ij * jk)	+ 
		    (80 * CF * Box6ijk * ij * jk)	+ 
		    (24 * CA * Box6ikj * ij * jk)	+ 
		    (2 * CA * Lij * ij * jk)		- 
		    (4 * CF * Lij * ij * jk)		- 
		    (CA * Lik * ij * jk)		- 
		    (4 * CF * Lik * ij * jk)		- 
		    (12 * CF * Ljk * ij * jk)		- 
		    (8 * CA * Box6ijk * ij3 * jk)	/ ik2 + 
		    (16 * CF * Box6ijk * ij3 * jk)	/ ik2 - 
		    (32 * CA * Box6ijk * ij2 * jk)	/ ik + 
		    (64 * CF * Box6ijk * ij2 * jk)	/ ik + 
		    (CA * Lij * ij2 * jk)		/ ik - 
		    (2 * CF * Lij * ij2 * jk)		/ ik + 
		    (CA * Ljk * ij2 * jk)		/ ik - 
		    (2 * CF * Ljk * ij2 * jk)		/ ik + 
		    (2 * CF * ik * jk)			- 
		    (16 * CA * Box6ijk * ik * jk)	+ 
		    (32 * CF * Box6ijk * ik * jk)	+ 
		    (48 * CA * Box6ikj * ik * jk)	+ 
		    (CA * Lij * ik * jk)		- 
		    (2 * CF * Lij * ik * jk)		- 
		    (2 * CA * Lik * ik * jk)		- 
		    (8 * CF * Lik * ik * jk)		- 
		    (CA * Ljk * ik * jk)		- 
		    (8 * CF * Ljk * ik * jk)		+ 
		    (24 * CA * Box6ikj * ik2 * jk)	/ ij - 
		    (CA * Lik * ik2 * jk)		/ ij - 
		    (4 * CF * Lik * ik2 * jk)		/ ij - 
		    (8 * CA * Box6ijk * jk2)		+ 
		    (16 * CF * Box6ijk * jk2)		+ 
		    (24 * CA * Box6ikj * jk2)		- 
		    (8 * CF * Ljk * jk2)		- 
		    (8 * CA * Box6ijk * ij2 * jk2)	/ ik2 + 
		    (16 * CF * Box6ijk * ij2 * jk2)	/ ik2 - 
		    (16 * CA * Box6ijk * ij * jk2)	/ ik + 
		    (32 * CF * Box6ijk * ij * jk2)	/ ik + 
		    (CA * Ljk * ij * jk2)		/ ik - 
		    (2 * CF * Ljk * ij * jk2)		/ ik + 
		    (48 * CA * Box6ikj * ik * jk2)	/ ij - 
		    (CA * Ljk * ik * jk2)		/ ij - 
		    (4 * CF * Ljk * ik * jk2)		/ ij + 
		    (24 * CA * Box6ikj * ik2 * jk2)	/ ij2
		    ) / (ij_ik_2 * ij_jk);

  qqbargLoops[1] = (
		    (-2 * CF * ij2)			+ 
		    (8 * CA * Box6ijk * ij2)		- 
		    (16 * CF * Box6ijk * ij2)		+ 
		    (32 * CA * Box6jik * ij2)		- 
		    (64 * CF * Box6jik * ij2)		- 
		    (2 * CA * Lij * ij2)		+ 
		    (4 * CF * Lij * ij2)		+ 
		    (4 * CF * Lik * ij2)		+ 
		    (CA * Ljk * ij2)			+ 
		    (2 * CF * Ljk * ij2)		+ 
		    (8 * CA * Box6ijk * ij3)		/ ik - 
		    (16 * CF * Box6ijk * ij3)		/ ik - 
		    (2 * CF * ij * ik)			- 
		    (24 * CA * Box6ikj * ij * ik)	+ 
		    (40 * CA * Box6jik * ij * ik)	- 
		    (80 * CF * Box6jik * ij * ik)	- 
		    (2 * CA * Lij * ij * ik)		+ 
		    (4 * CF * Lij * ij * ik)		+ 
		    (12 * CF * Lik * ij * ik)		+ 
		    (CA * Ljk * ij * ik)		+ 
		    (4 * CF * Ljk * ij * ik)		- 
		    (24 * CA * Box6ikj * ik2)		+ 
		    (8 * CA * Box6jik * ik2)		- 
		    (16 * CF * Box6jik * ik2)		+ 
		    (8 * CF * Lik * ik2)		+ 
		    (8 * CA * Box6jik * ij3 * ik)	/ jk2 - 
		    (16 * CF * Box6jik * ij3 * ik)	/ jk2 + 
		    (8 * CA * Box6jik * ij2 * ik2)	/ jk2 - 
		    (16 * CF * Box6jik * ij2 * ik2)	/ jk2 + 
		    (16 * CA * Box6jik * ij3)		/ jk - 
		    (32 * CF * Box6jik * ij3)		/ jk - 
		    (CA * Lij * ij3)			/ jk + 
		    (2 * CF * Lij * ij3)		/ jk + 
		    (32 * CA * Box6jik * ij2 * ik)	/ jk - 
		    (64 * CF * Box6jik * ij2 * ik)	/ jk - 
		    (CA * Lij * ij2 * ik)		/ jk + 
		    (2 * CF * Lij * ij2 * ik)		/ jk - 
		    (CA * Lik * ij2 * ik)		/ jk + 
		    (2 * CF * Lik * ij2 * ik)		/ jk + 
		    (16 * CA * Box6jik * ij * ik2)	/ jk - 
		    (32 * CF * Box6jik * ij * ik2)	/ jk - 
		    (CA * Lik * ij * ik2)		/ jk + 
		    (2 * CF * Lik * ij * ik2)		/ jk - 
		    (2 * CF * ij * jk)			+ 
		    (16 * CA * Box6ijk * ij * jk)	- 
		    (32 * CF * Box6ijk * ij * jk)	+ 
		    (16 * CA * Box6jik * ij * jk)	- 
		    (32 * CF * Box6jik * ij * jk)	- 
		    (CA * Lij * ij * jk)		+ 
		    (2 * CF * Lij * ij * jk)		+ 
		    (4 * CF * Lik * ij * jk)		+ 
		    (2 * CA * Ljk * ij * jk)		+ 
		    (4 * CF * Ljk * ij * jk)		+ 
		    (16 * CA * Box6ijk * ij2 * jk)	/ ik - 
		    (32 * CF * Box6ijk * ij2 * jk)	/ ik - 
		    (2 * CF * ik * jk)			- 
		    (48 * CA * Box6ikj * ik * jk)	+ 
		    (16 * CA * Box6jik * ik * jk)	- 
		    (32 * CF * Box6jik * ik * jk)	- 
		    (CA * Lij * ik * jk)		+ 
		    (2 * CF * Lij * ik * jk)		+ 
		    (CA * Lik * ik * jk)		+ 
		    (8 * CF * Lik * ik * jk)		+ 
		    (2 * CA * Ljk * ik * jk)		+ 
		    (8 * CF * Ljk * ik * jk)		- 
		    (48 * CA * Box6ikj * ik2 * jk)	/ ij + 
		    (CA * Lik * ik2 * jk)		/ ij + 
		    (4 * CF * Lik * ik2 * jk)		/ ij + 
		    (8 * CA * Box6ijk * jk2)		- 
		    (16 * CF * Box6ijk * jk2)		+ 
		    (CA * Ljk * jk2)			+ 
		    (2 * CF * Ljk * jk2)		+ 
		    (8 * CA * Box6ijk * ij * jk2)	/ ik - 
		    (16 * CF * Box6ijk * ij * jk2)	/ ik - 
		    (24 * CA * Box6ikj * ik * jk2)	/ ij + 
		    (CA * Ljk * ik * jk2)		/ ij + 
		    (4 * CF * Ljk * ik * jk2)		/ ij - 
		    (24 * CA * Box6ikj * ik2 * jk2)	/ ij2
		    ) / (ij_ik * ij_jk_2);

  qqbargLoops[2] = -3 * CF * Lijk + (
				     (-4 * CA * Box6jik * ij3)		+ 
				     (8 * CF * Box6jik * ij3)		+ 
				     (CA * Lij * ij3)			/ 2. -
				     (CF * Lij * ij3)			+ 
				     (CA * ij2 * ik)			- 
				     (9 * CF * ij2 * ik)			+ 
				     (8 * CA * Box6ijk * ij2 * ik)	- 
				     (16 * CF * Box6ijk * ij2 * ik)	- 
				     (8 * CA * Box6ikj * ij2 * ik)	- 
				     (8 * CA * Box6jik * ij2 * ik)	+ 
				     (16 * CF * Box6jik * ij2 * ik)	+ 
				     (CA * Lij * ij2 * ik)		/ 2. -
				     (CF * Lij * ij2 * ik)		+ 
				     (CA * Lik * ij2 * ik)		/ 2. -
				     (CF * Lik * ij2 * ik)		+ 
				     (CA * ij * ik2)			- 
				     (9 * CF * ij * ik2)			+ 
				     (8 * CA * Box6ijk * ij * ik2)	- 
				     (16 * CF * Box6ijk * ij * ik2)	- 
				     (8 * CA * Box6ikj * ij * ik2)	- 
				     (4 * CA * Box6jik * ij * ik2)	+ 
				     (8 * CF * Box6jik * ij * ik2)	+ 
				     (CA * Lik * ij * ik2)		/ 2. -
				     (CF * Lik * ij * ik2)		- 
				     (4 * CA * Box6jik * ij3 * ik)	/ jk + 
				     (8 * CF * Box6jik * ij3 * ik)	/ jk - 
				     (4 * CA * Box6jik * ij2 * ik2)	/ jk + 
				     (8 * CF * Box6jik * ij2 * ik2)	/ jk + 
				     (CA * ij2 * jk)			- 
				     (9 * CF * ij2 * jk)			+ 
				     (12 * CA * Box6ijk * ij2 * jk)	- 
				     (24 * CF * Box6ijk * ij2 * jk)	- 
				     (8 * CA * Box6ikj * ij2 * jk)	- 
				     (4 * CA * Box6jik * ij2 * jk)	+ 
				     (8 * CF * Box6jik * ij2 * jk)	+ 
				     (CA * Lik * ij2 * jk)		/ 2. -
				     (CF * Lik * ij2 * jk)		- 
				     (CA * Ljk * ij2 * jk)		/ 2. +
				     (CF * Ljk * ij2 * jk)		+ 
				     (4 * CA * Box6ijk * ij3 * jk)	/ ik - 
				     (8 * CF * Box6ijk * ij3 * jk)	/ ik - 
				     (CA * Lij * ij3 * jk)		/ (2. * ik) + 
				     (CF * Lij * ij3 * jk)		/ ik + 
				     (2 * CA * ij * ik * jk)		- 
				     (18 * CF * ij * ik * jk)		+ 
				     (16 * CA * Box6ijk * ij * ik * jk)	- 
				     (32 * CF * Box6ijk * ij * ik * jk)	- 
				     (28 * CA * Box6ikj * ij * ik * jk)	- 
				     (4 * CA * Box6jik * ij * ik * jk)	+ 
				     (8 * CF * Box6jik * ij * ik * jk)	+ 
				     (CA * Lij * ij * ik * jk)		/ 2. -
				     (CF * Lij * ij * ik * jk)		+ 
				     (CA * Lik * ij * ik * jk)		- 
				     (CF * Lik * ij * ik * jk)		- 
				     (CA * Ljk * ij * ik * jk)		/ 2. +
				     (3 * CF * Ljk * ij * ik * jk)	+ 
				     (CA * ik2 * jk)			- 
				     (9 * CF * ik2 * jk)			+ 
				     (8 * CA * Box6ijk * ik2 * jk)	- 
				     (16 * CF * Box6ijk * ik2 * jk)	- 
				     (20 * CA * Box6ikj * ik2 * jk)	+ 
				     (CA * Lik * ik2 * jk)		/ 2. +
				     (CA * ij * jk2)			- 
				     (9 * CF * ij * jk2)			+ 
				     (12 * CA * Box6ijk * ij * jk2)	- 
				     (24 * CF * Box6ijk * ij * jk2)	- 
				     (20 * CA * Box6ikj * ij * jk2)	- 
				     (CA * Lij * ij * jk2)		/ 2. + 
				     (CF * Lij * ij * jk2)		+ 
				     (CA * Lik * ij * jk2)		/ 2. - 
				     (CA * Ljk * ij * jk2)		+ 
				     (4 * CF * Ljk * ij * jk2)		+ 
				     (4 * CA * Box6ijk * ij3 * jk2)	/ ik2 - 
				     (8 * CF * Box6ijk * ij3 * jk2)	/ ik2 + 
				     (8 * CA * Box6ijk * ij2 * jk2)	/ ik - 
				     (16 * CF * Box6ijk * ij2 * jk2)	/ ik - 
				     (CA * Lij * ij2 * jk2)		/ (2. * ik) + 
				     (CF * Lij * ij2 * jk2)		/ ik - 
				     (CA * Ljk * ij2 * jk2)		/ (2. * ik) + 
				     (CF * Ljk * ij2 * jk2)		/ ik + 
				     (CA * ik * jk2)			- 
				     (9 * CF * ik * jk2)			+ 
				     (8 * CA * Box6ijk * ik * jk2)	- 
				     (16 * CF * Box6ijk * ik * jk2)	- 
				     (32 * CA * Box6ikj * ik * jk2)	+ 
				     (CA * Lik * ik * jk2)		/ 2. - 
				     (CA * Ljk * ik * jk2)		/ 2. + 
				     (3 * CF * Ljk * ik * jk2)		- 
				     (12 * CA * Box6ikj * ik2 * jk2)	/ ij - 
				     (12 * CA * Box6ikj * jk3)		- 
				     (CA * Ljk * jk3)			/ 2. +
				     (3 * CF * Ljk * jk3)		+ 
				     (4 * CA * Box6ijk * ij2 * jk3)	/ ik2 - 
				     (8 * CF * Box6ijk * ij2 * jk3)	/ ik2 + 
				     (4 * CA * Box6ijk * ij * jk3)	/ ik - 
				     (8 * CF * Box6ijk * ij * jk3)	/ ik - 
				     (CA * Ljk * ij * jk3)		/ (2. * ik) + 
				     (CF * Ljk * ij * jk3)		/ ik - 
				     (12 * CA * Box6ikj * ik * jk3)	/ ij
				     ) / (ij_ik * ij_jk * ik_jk);

  qqbargLoops[3] = 3 * CF * Lijk + (
				    (8 * CF * ij2)			- 
				    (8 * CA * Box6ijk * ij2)		+ 
				    (16 * CF * Box6ijk * ij2)		+ 
				    (8 * CA * Box6ikj * ij2)		- 
				    (8 * CA * Box6jik * ij2)		+ 
				    (16 * CF * Box6jik * ij2)		+ 
				    (CA * Lij * ij2)			/ 2. - 
				    (CF * Lij * ij2)			+ 
				    (8 * CF * ij * ik)			- 
				    (8 * CA * Box6ijk * ij * ik)	+ 
				    (16 * CF * Box6ijk * ij * ik)	+ 
				    (8 * CA * Box6ikj * ij * ik)	- 
				    (12 * CA * Box6jik * ij * ik)	+ 
				    (24 * CF * Box6jik * ij * ik)	+ 
				    (CA * Lij * ij * ik)		/ 2. - 
				    (CF * Lij * ij * ik)		+ 
				    (CA * Lik * ij * ik)		/ 2. - 
				    (CF * Lik * ij * ik)		- 
				    (4 * CA * Box6jik * ik2)		+ 
				    (8 * CF * Box6jik * ik2)		+ 
				    (CA * Lik * ik2)			/ 2. - 
				    (CF * Lik * ik2)			- 
				    (4 * CA * Box6jik * ij2 * ik)	/ jk + 
				    (8 * CF * Box6jik * ij2 * ik)	/ jk - 
				    (4 * CA * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * ij * jk)			- 
				    (12 * CA * Box6ijk * ij * jk)	+ 
				    (24 * CF * Box6ijk * ij * jk)	+ 
				    (8 * CA * Box6ikj * ij * jk)	- 
				    (8 * CA * Box6jik * ij * jk)	+ 
				    (16 * CF * Box6jik * ij * jk)	+ 
				    (CA * Lij * ij * jk)		/ 2. - 
				    (CF * Lij * ij * jk)		+ 
				    (CA * Ljk * ij * jk)		/ 2. - 
				    (CF * Ljk * ij * jk)		- 
				    (4 * CA * Box6ijk * ij2 * jk)	/ ik + 
				    (8 * CF * Box6ijk * ij2 * jk)	/ ik + 
				    (8 * CF * ik * jk)			- 
				    (8 * CA * Box6ijk * ik * jk)	+ 
				    (16 * CF * Box6ijk * ik * jk)	- 
				    (4 * CA * Box6ikj * ik * jk)	- 
				    (8 * CA * Box6jik * ik * jk)	+ 
				    (16 * CF * Box6jik * ik * jk)	+ 
				    (CA * Lij * ik * jk)		/ 2. - 
				    (CF * Lij * ik * jk)		+ 
				    (CA * Lik * ik * jk)		/ 2. + 
				    (2 * CF * Lik * ik * jk)		+ 
				    (CA * Ljk * ik * jk)		/ 2. + 
				    (2 * CF * Ljk * ik * jk)		- 
				    (12 * CA * Box6ikj * ik2 * jk)	/ ij + 
				    (CA * Lik * ik2 * jk)		/ (2. * ij) + 
				    (2 * CF * Lik * ik2 * jk)		/ ij - 
				    (4 * CA * Box6ijk * jk2)		+ 
				    (8 * CF * Box6ijk * jk2)		+ 
				    (CA * Ljk * jk2)			/ 2. - 
				    (CF * Ljk * jk2)			- 
				    (4 * CA * Box6ijk * ij * jk2)	/ ik + 
				    (8 * CF * Box6ijk * ij * jk2)	/ ik - 
				    (12 * CA * Box6ikj * ik * jk2)	/ ij + 
				    (CA * Ljk * ik * jk2)		/ (2. * ij) + 
				    (2 * CF * Ljk * ik * jk2)		/ ij - 
				    (12 * CA * Box6ikj * ik2 * jk2)	/ ij2
				    ) / (ij_ik * ij_jk);

  qqbargLoops[4] = -3 * CF * Lijk + (
				     (-8 * CF * ij2)			+ 
				     (8 * CA * Box6ijk * ij2)		- 
				     (16 * CF * Box6ijk * ij2)		- 
				     (8 * CA * Box6ikj * ij2)		+ 
				     (8 * CA * Box6jik * ij2)		- 
				     (16 * CF * Box6jik * ij2)		- 
				     (CA * Lij * ij2)			/ 2. + 
				     (CF * Lij * ij2)			- 
				     (8 * CF * ij * ik)			+ 
				     (8 * CA * Box6ijk * ij * ik)	- 
				     (16 * CF * Box6ijk * ij * ik)	- 
				     (8 * CA * Box6ikj * ij * ik)	+ 
				     (12 * CA * Box6jik * ij * ik)	- 
				     (24 * CF * Box6jik * ij * ik)	- 
				     (CA * Lij * ij * ik)		/ 2. + 
				     (CF * Lij * ij * ik)		- 
				     (CA * Lik * ij * ik)		/ 2. + 
				     (CF * Lik * ij * ik)		+ 
				     (4 * CA * Box6jik * ik2)		- 
				     (8 * CF * Box6jik * ik2)		- 
				     (CA * Lik * ik2)			/ 2. + 
				     (CF * Lik * ik2)			+ 
				     (4 * CA * Box6jik * ij2 * ik)	/ jk - 
				     (8 * CF * Box6jik * ij2 * ik)	/ jk + 
				     (4 * CA * Box6jik * ij * ik2)	/ jk - 
				     (8 * CF * Box6jik * ij * ik2)	/ jk - 
				     (8 * CF * ij * jk)			+ 
				     (12 * CA * Box6ijk * ij * jk)	- 
				     (24 * CF * Box6ijk * ij * jk)	- 
				     (8 * CA * Box6ikj * ij * jk)	+ 
				     (8 * CA * Box6jik * ij * jk)	- 
				     (16 * CF * Box6jik * ij * jk)	- 
				     (CA * Lij * ij * jk)		/ 2. + 
				     (CF * Lij * ij * jk)		- 
				     (CA * Ljk * ij * jk)		/ 2. + 
				     (CF * Ljk * ij * jk)		+ 
				     (4 * CA * Box6ijk * ij2 * jk)	/ ik - 
				     (8 * CF * Box6ijk * ij2 * jk)	/ ik - 
				     (8 * CF * ik * jk)			+ 
				     (8 * CA * Box6ijk * ik * jk)	- 
				     (16 * CF * Box6ijk * ik * jk)	+ 
				     (4 * CA * Box6ikj * ik * jk)	+ 
				     (8 * CA * Box6jik * ik * jk)	- 
				     (16 * CF * Box6jik * ik * jk)	- 
				     (CA * Lij * ik * jk)		/ 2. + 
				     (CF * Lij * ik * jk)		- 
				     (CA * Lik * ik * jk)		/ 2. - 
				     (2 * CF * Lik * ik * jk)		- 
				     (CA * Ljk * ik * jk)		/ 2. - 
				     (2 * CF * Ljk * ik * jk)		+ 
				     (12 * CA * Box6ikj * ik2 * jk)	/ ij - 
				     (CA * Lik * ik2 * jk)		/ (2. * ij) - 
				     (2 * CF * Lik * ik2 * jk)		/ ij + 
				     (4 * CA * Box6ijk * jk2)		- 
				     (8 * CF * Box6ijk * jk2)		- 
				     (CA * Ljk * jk2)			/ 2. + 
				     (CF * Ljk * jk2)			+ 
				     (4 * CA * Box6ijk * ij * jk2)	/ ik - 
				     (8 * CF * Box6ijk * ij * jk2)	/ ik + 
				     (12 * CA * Box6ikj * ik * jk2)	/ ij - 
				     (CA * Ljk * ik * jk2)		/ (2. * ij) - 
				     (2 * CF * Ljk * ik * jk2)		/ ij + 
				     (12 * CA * Box6ikj * ik2 * jk2)	/ ij2
				     ) / (ij_ik * ij_jk);

  qqbargLoops[5] = 3 * CF * Lijk + (
				    (-4 * CA * Box6jik * ij2)		+ 
				    (8 * CF * Box6jik * ij2)		+ 
				    (CA * Lij * ij2)			/ 2. - 
				    (CF * Lij * ij2)			- 
				    (CA * ij * ik)			+ 
				    (9 * CF * ij * ik)			- 
				    (8 * CA * Box6ijk * ij * ik)	+ 
				    (16 * CF * Box6ijk * ij * ik)	+ 
				    (8 * CA * Box6ikj * ij * ik)	- 
				    (4 * CA * Box6jik * ij * ik)	+ 
				    (8 * CF * Box6jik * ij * ik)	+ 
				    (CA * Lij * ij * ik)		/ 2. - 
				    (CF * Lij * ij * ik)		+ 
				    (CA * Lik * ij * ik)		/ 2. - 
				    (CF * Lik * ij * ik)		- 
				    (CA * ik2)				+ 
				    (9 * CF * ik2)			- 
				    (8 * CA * Box6ijk * ik2)		+ 
				    (16 * CF * Box6ijk * ik2)		+ 
				    (8 * CA * Box6ikj * ik2)		+ 
				    (CA * Lik * ik2)			/ 2. - 
				    (CF * Lik * ik2)			- 
				    (4 * CA * Box6jik * ij2 * ik)	/ jk + 
				    (8 * CF * Box6jik * ij2 * ik)	/ jk - 
				    (4 * CA * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * Box6jik * ij * ik2)	/ jk - 
				    (CA * ij * jk)			+ 
				    (9 * CF * ij * jk)			- 
				    (4 * CA * Box6ijk * ij * jk)	+ 
				    (8 * CF * Box6ijk * ij * jk)	+ 
				    (8 * CA * Box6ikj * ij * jk)	- 
				    (CA * Lij * ij * jk)		/ 2. + 
				    (CF * Lij * ij * jk)		+ 
				    (CA * Lik * ij * jk)		/ 2. - 
				    (CF * Lik * ij * jk)		- 
				    (CA * Ljk * ij * jk)		/ 2. + 
				    (CF * Ljk * ij * jk)		+ 
				    (4 * CA * Box6ijk * ij2 * jk)	/ ik - 
				    (8 * CF * Box6ijk * ij2 * jk)	/ ik - 
				    (CA * Lij * ij2 * jk)		/ (2. * ik) + 
				    (CF * Lij * ij2 * jk)		/ ik - 
				    (CA * ik * jk)			+ 
				    (9 * CF * ik * jk)			- 
				    (8 * CA * Box6ijk * ik * jk)	+ 
				    (16 * CF * Box6ijk * ik * jk)	+ 
				    (20 * CA * Box6ikj * ik * jk)	+ 
				    (CA * Lik * ik * jk)		/ 2. - 
				    (CF * Lik * ik * jk)		- 
				    (CA * Ljk * ik * jk)		/ 2. - 
				    (2 * CF * Ljk * ik * jk)		+ 
				    (12 * CA * Box6ikj * ik2 * jk)	/ ij + 
				    (12 * CA * Box6ikj * jk2)		- 
				    (CA * Ljk * jk2)			/ 2. - 
				    (2 * CF * Ljk * jk2)		+ 
				    (4 * CA * Box6ijk * ij2 * jk2)	/ ik2 - 
				    (8 * CF * Box6ijk * ij2 * jk2)	/ ik2 + 
				    (4 * CA * Box6ijk * ij * jk2)	/ ik - 
				    (8 * CF * Box6ijk * ij * jk2)	/ ik - 
				    (CA * Ljk * ij * jk2)		/ (2. * ik) + 
				    (CF * Ljk * ij * jk2)		/ ik + 
				    (12 * CA * Box6ikj * ik * jk2)	/ ij
				    ) / (ij_ik * ik_jk);

  qqbargLoops[6] = (
		    (-2 * CF * ij)			+ 
		    (32 * CA * Box6ijk * ij)		- 
		    (64 * CF * Box6ijk * ij)		- 
		    (4 * CA * Lij * ij)			+ 
		    (8 * CF * Lij * ij)			+ 
		    (4 * CF * Ljk * ij)			+ 
		    (16 * CA * Box6ijk * ij2)		/ ik - 
		    (32 * CF * Box6ijk * ij2)		/ ik - 
		    (2 * CA * Lij * ij2)		/ ik + 
		    (4 * CF * Lij * ij2)		/ ik - 
		    (2 * CF * ik)			+ 
		    (16 * CA * Box6ijk * ik)		- 
		    (32 * CF * Box6ijk * ik)		- 
		    (2 * CA * Lij * ik)			+ 
		    (4 * CF * Lij * ik)			+ 
		    (4 * CF * Ljk * ik)			+ 
		    (16 * CA * Box6ijk * jk)		- 
		    (32 * CF * Box6ijk * jk)		- 
		    (2 * CA * Ljk * jk)			+ 
		    (6 * CF * Ljk * jk)			+ 
		    (16 * CA * Box6ijk * ij2 * jk)	/ ik2 - 
		    (32 * CF * Box6ijk * ij2 * jk)	/ ik2 + 
		    (32 * CA * Box6ijk * ij * jk)	/ ik - 
		    (64 * CF * Box6ijk * ij * jk)	/ ik - 
		    (2 * CA * Ljk * ij * jk)		/ ik + 
		    (4 * CF * Ljk * ij * jk)		/ ik
		    ) / ij_ik_2;

  qqbargLoops[7] = (
		    (8 * CA * Box6jik * ij)		- 
		    (16 * CF * Box6jik * ij)		+ 
		    (CA * Lij * ij)			- 
		    (2 * CF * Lij * ij)			+ 
		    (CA * Lik * ij)			+ 
		    (2 * CF * Lik * ij)			+ 
		    (CA * Lij * ij2)			/ ik - 
		    (2 * CF * Lij * ij2)		/ ik + 
		    (8 * CA * Box6jik * ik)		- 
		    (16 * CF * Box6jik * ik)		+ 
		    (CA * Lik * ik)			+ 
		    (2 * CF * Lik * ik)			+ 
		    (8 * CA * Box6jik * ij2)		/ jk - 
		    (16 * CF * Box6jik * ij2)		/ jk + 
		    (8 * CA * Box6jik * ij * ik)	/ jk - 
		    (16 * CF * Box6jik * ij * ik)	/ jk - 
		    (24 * CA * Box6ikj * jk)		+ 
		    (CA * Lij * jk)			- 
		    (2 * CF * Lij * jk)			+ 
		    (CA * Lik * jk)			+ 
		    (4 * CF * Lik * jk)			+ 
		    (CA * Ljk * jk)			+ 
		    (4 * CF * Ljk * jk)			- 
		    (8 * CA * Box6ijk * ij2 * jk)	/ ik2 + 
		    (16 * CF * Box6ijk * ij2 * jk)	/ ik2 - 
		    (8 * CA * Box6ijk * ij * jk)	/ ik + 
		    (16 * CF * Box6ijk * ij * jk)	/ ik + 
		    (CA * Lij * ij * jk)		/ ik - 
		    (2 * CF * Lij * ij * jk)		/ ik + 
		    (CA * Ljk * ij * jk)		/ ik - 
		    (2 * CF * Ljk * ij * jk)		/ ik - 
		    (24 * CA * Box6ikj * ik * jk)	/ ij + 
		    (CA * Lik * ik * jk)		/ ij + 
		    (4 * CF * Lik * ik * jk)		/ ij - 
		    (24 * CA * Box6ikj * jk2)		/ ij + 
		    (CA * Ljk * jk2)			/ ij + 
		    (4 * CF * Ljk * jk2)		/ ij - 
		    (8 * CA * Box6ijk * ij * jk2)	/ ik2 + 
		    (16 * CF * Box6ijk * ij * jk2)	/ ik2 - 
		    (8 * CA * Box6ijk * jk2)		/ ik + 
		    (16 * CF * Box6ijk * jk2)		/ ik + 
		    (CA * Ljk * jk2)			/ ik - 
		    (2 * CF * Ljk * jk2)		/ ik - 
		    (24 * CA * Box6ikj * ik * jk2)	/ ij2
		    ) / (ij_ik * ij_jk);

  qqbargLoops[8] = (
		    (-8 * CA * Box6ijk * ij)		+ 
		    (16 * CF * Box6ijk * ij)		- 
		    (CA * Lij * ij)			+ 
		    (2 * CF * Lij * ij)			- 
		    (CA * Ljk * ij)			- 
		    (2 * CF * Ljk * ij)			- 
		    (8 * CA * Box6ijk * ij2)		/ ik + 
		    (16 * CF * Box6ijk * ij2)		/ ik + 
		    (24 * CA * Box6ikj * ik)		- 
		    (CA * Lij * ik)			+ 
		    (2 * CF * Lij * ik)			- 
		    (CA * Lik * ik)			- 
		    (4 * CF * Lik * ik)			- 
		    (CA * Ljk * ik)			- 
		    (4 * CF * Ljk * ik)			+ 
		    (24 * CA * Box6ikj * ik2)		/ ij - 
		    (CA * Lik * ik2)			/ ij - 
		    (4 * CF * Lik * ik2)		/ ij + 
		    (8 * CA * Box6jik * ij2 * ik)	/ jk2 - 
		    (16 * CF * Box6jik * ij2 * ik)	/ jk2 + 
		    (8 * CA * Box6jik * ij * ik2)	/ jk2 - 
		    (16 * CF * Box6jik * ij * ik2)	/ jk2 - 
		    (CA * Lij * ij2)			/ jk + 
		    (2 * CF * Lij * ij2)		/ jk + 
		    (8 * CA * Box6jik * ij * ik)	/ jk - 
		    (16 * CF * Box6jik * ij * ik)	/ jk - 
		    (CA * Lij * ij * ik)		/ jk + 
		    (2 * CF * Lij * ij * ik)		/ jk - 
		    (CA * Lik * ij * ik)		/ jk + 
		    (2 * CF * Lik * ij * ik)		/ jk + 
		    (8 * CA * Box6jik * ik2)		/ jk - 
		    (16 * CF * Box6jik * ik2)		/ jk - 
		    (CA * Lik * ik2)			/ jk + 
		    (2 * CF * Lik * ik2)		/ jk - 
		    (8 * CA * Box6ijk * jk)		+ 
		    (16 * CF * Box6ijk * jk)		- 
		    (CA * Ljk * jk)			- 
		    (2 * CF * Ljk * jk)			- 
		    (8 * CA * Box6ijk * ij * jk)	/ ik + 
		    (16 * CF * Box6ijk * ij * jk)	/ ik + 
		    (24 * CA * Box6ikj * ik * jk)	/ ij - 
		    (CA * Ljk * ik * jk)		/ ij - 
		    (4 * CF * Ljk * ik * jk)		/ ij + 
		    (24 * CA * Box6ikj * ik2 * jk)	/ ij2
		    ) / (ij_ik * ij_jk);

  qqbargLoops[9] = (
		    (2 * CF * ij) 			- 
		    (32 * CA * Box6jik * ij) 		+ 
		    (64 * CF * Box6jik * ij) 		+ 
		    (4 * CA * Lij * ij) 		- 
		    (8 * CF * Lij * ij) 		- 
		    (4 * CF * Lik * ij) 		- 
		    (16 * CA * Box6jik * ik)		+ 
		    (32 * CF * Box6jik * ik)		+ 
		    (2 * CA * Lik * ik) 		- 
		    (6 * CF * Lik * ik) 		- 
		    (16 * CA * Box6jik * ij2 * ik)	/ jk2 + 
		    (32 * CF * Box6jik * ij2 * ik)	/ jk2 - 
		    (16 * CA * Box6jik * ij2)		/ jk + 
		    (32 * CF * Box6jik * ij2)		/ jk + 
		    (2 * CA * Lij * ij2)		/ jk - 
		    (4 * CF * Lij * ij2)		/ jk - 
		    (32 * CA * Box6jik * ij * ik)	/ jk + 
		    (64 * CF * Box6jik * ij * ik)	/ jk + 
		    (2 * CA * Lik * ij * ik)		/ jk - 
		    (4 * CF * Lik * ij * ik)		/ jk + 
		    (2 * CF * jk) - 
		    (16 * CA * Box6jik * jk) + 
		    (32 * CF * Box6jik * jk) + 
		    (2 * CA * Lij * jk) - 
		    (4 * CF * Lij * jk) - 
		    (4 * CF * Lik * jk)
		    ) / ij_jk_2;

  qqbargLoops[10] = (
		     (-8 * CA * Box6ijk * ij2 * jk)	+ 
		     (16 * CF * Box6ijk * ij2 * jk)	+ 
		     (2 * CA * Lij * ij2 * jk)		- 
		     (4 * CF * Lij * ij2 * jk)		- 
		     (CA * ij * ik * jk)			+ 
		     (2 * CF * ij * ik * jk)		- 
		     (8 * CA * Box6ijk * ij * ik * jk)	+ 
		     (16 * CF * Box6ijk * ij * ik * jk)	+ 
		     (3 * CA * Lij * ij * ik * jk)	- 
		     (6 * CF * Lij * ij * ik * jk)	+ 
		     (CA * Ljk * ij * ik * jk)		- 
		     (2 * CF * Ljk * ij * ik * jk)	- 
		     (CA * ik2 * jk)			+ 
		     (2 * CF * ik2 * jk)			+ 
		     (CA * Lij * ik2 * jk)		- 
		     (2 * CF * Lij * ik2 * jk)		+ 
		     (CA * Ljk * ik2 * jk)		- 
		     (CF * Ljk * ik2 * jk)		- 
		     (CA * ij * jk2)			+ 
		     (2 * CF * ij * jk2)			- 
		     (16 * CA * Box6ijk * ij * jk2)	+ 
		     (32 * CF * Box6ijk * ij * jk2)	+ 
		     (2 * CA * Lij * ij * jk2)		- 
		     (4 * CF * Lij * ij * jk2)		+ 
		     (2 * CA * Ljk * ij * jk2)		- 
		     (4 * CF * Ljk * ij * jk2)		- 
		     (16 * CA * Box6ijk * ij2 * jk2)	/ ik + 
		     (32 * CF * Box6ijk * ij2 * jk2)	/ ik + 
		     (CA * Lij * ij2 * jk2)		/ ik - 
		     (2 * CF * Lij * ij2 * jk2)		/ ik - 
		     (CA * ik * jk2)			+ 
		     (2 * CF * ik * jk2)			+ 
		     (CA * Lij * ik * jk2)		- 
		     (2 * CF * Lij * ik * jk2)		+ 
		     (2 * CA * Ljk * ik * jk2)		- 
		     (2 * CF * Ljk * ik * jk2)		+ 
		     (CA * Ljk * jk3)			- 
		     (CF * Ljk * jk3)			- 
		     (8 * CA * Box6ijk * ij2 * jk3)	/ ik2 + 
		     (16 * CF * Box6ijk * ij2 * jk3)	/ ik2 - 
		     (8 * CA * Box6ijk * ij * jk3)	/ ik + 
		     (16 * CF * Box6ijk * ij * jk3)	/ ik + 
		     (CA * Ljk * ij * jk3)		/ ik - 
		     (2 * CF * Ljk * ij * jk3)		/ ik
		     ) / (ij_ik * ik_jk_2);

  qqbargLoops[11] = (
		     (16 * CA * Box6jik * ij2 * ik)	- 
		     (32 * CF * Box6jik * ij2 * ik)	- 
		     (CA * Lij * ij2 * ik)		+ 
		     (2 * CF * Lij * ij2 * ik)		+ 
		     (8 * CA * Box6jik * ij * ik2)	- 
		     (16 * CF * Box6jik * ij * ik2)	- 
		     (CA * Lik * ij * ik2)		+ 
		     (2 * CF * Lik * ij * ik2)		+ 
		     (8 * CA * Box6jik * ij2 * ik2)	/ jk - 
		     (16 * CF * Box6jik * ij2 * ik2)	/ jk + 
		     (8 * CA * Box6jik * ij2 * jk)	- 
		     (16 * CF * Box6jik * ij2 * jk)	- 
		     (2 * CA * Lij * ij2 * jk)		+ 
		     (4 * CF * Lij * ij2 * jk)		+ 
		     (CA * ij * ik * jk)			- 
		     (2 * CF * ij * ik * jk)		+ 
		     (16 * CA * Box6jik * ij * ik * jk)	- 
		     (32 * CF * Box6jik * ij * ik * jk)	- 
		     (2 * CA * Lij * ij * ik * jk)	+ 
		     (4 * CF * Lij * ij * ik * jk)	- 
		     (2 * CA * Lik * ij * ik * jk)	+ 
		     (4 * CF * Lik * ij * ik * jk)	- 
		     (CA * Lik * ik2 * jk)		+ 
		     (CF * Lik * ik2 * jk)		+ 
		     (CA * ij * jk2)			- 
		     (2 * CF * ij * jk2)			+ 
		     (8 * CA * Box6jik * ij * jk2)	- 
		     (16 * CF * Box6jik * ij * jk2)	- 
		     (3 * CA * Lij * ij * jk2)		+ 
		     (6 * CF * Lij * ij * jk2)		- 
		     (CA * Lik * ij * jk2)		+ 
		     (2 * CF * Lik * ij * jk2)		+ 
		     (CA * ik * jk2)			- 
		     (2 * CF * ik * jk2)			- 
		     (CA * Lij * ik * jk2)		+ 
		     (2 * CF * Lij * ik * jk2)		- 
		     (2 * CA * Lik * ik * jk2)		+ 
		     (2 * CF * Lik * ik * jk2)		+ 
		     (CA * jk3)				- 
		     (2 * CF * jk3)			- 
		     (CA * Lij * jk3)			+ 
		     (2 * CF * Lij * jk3)		- 
		     (CA * Lik * jk3)		 	+ 
		     (CF * Lik * jk3)
		     ) / (ij_jk * ik_jk_2);

  qqbargLoops[12] = -3 * CF * Lijk + (
				      (CA * ij2 * ik)			- 
				      (9 * CF * ij2 * ik)			+ 
				      (8 * CA * Box6ijk * ij2 * ik)	- 
				      (16 * CF * Box6ijk * ij2 * ik)	- 
				      (8 * CA * Box6ikj * ij2 * ik)	+ 
				      (CA * ij * ik2)			- 
				      (9 * CF * ij * ik2)			+ 
				      (8 * CA * Box6ijk * ij * ik2)	- 
				      (16 * CF * Box6ijk * ij * ik2)	- 
				      (8 * CA * Box6ikj * ij * ik2)	+ 
				      (CA * ij2 * jk)			- 
				      (9 * CF * ij2 * jk)			- 
				      (8 * CA * Box6ikj * ij2 * jk)	+ 
				      (8 * CA * Box6jik * ij2 * jk)	- 
				      (16 * CF * Box6jik * ij2 * jk)	+ 
				      (2 * CA * ij * ik * jk)		- 
				      (18 * CF * ij * ik * jk)		+ 
				      (8 * CA * Box6ijk * ij * ik * jk)	- 
				      (16 * CF * Box6ijk * ij * ik * jk)	- 
				      (40 * CA * Box6ikj * ij * ik * jk)	+ 
				      (8 * CA * Box6jik * ij * ik * jk)	- 
				      (16 * CF * Box6jik * ij * ik * jk)	+ 
				      (3 * CF * Lik * ij * ik * jk)	+ 
				      (3 * CF * Ljk * ij * ik * jk)	+ 
				      (CA * ik2 * jk)			- 
				      (9 * CF * ik2 * jk)			+ 
				      (8 * CA * Box6ijk * ik2 * jk)	- 
				      (16 * CF * Box6ijk * ik2 * jk)	- 
				      (32 * CA * Box6ikj * ik2 * jk)	+ 
				      (3 * CF * Lik * ik2 * jk)		+ 
				      (CA * ij * jk2)			- 
				      (9 * CF * ij * jk2)			- 
				      (8 * CA * Box6ikj * ij * jk2)	+ 
				      (8 * CA * Box6jik * ij * jk2)	- 
				      (16 * CF * Box6jik * ij * jk2)	+ 
				      (CA * ik * jk2)			- 
				      (9 * CF * ik * jk2)			- 
				      (32 * CA * Box6ikj * ik * jk2)	+ 
				      (8 * CA * Box6jik * ik * jk2)	- 
				      (16 * CF * Box6jik * ik * jk2)	+ 
				      (3 * CF * Ljk * ik * jk2)		- 
				      (24 * CA * Box6ikj * ik2 * jk2)	/ ij
				      ) / (ij_ik * ij_jk * ik_jk);

  /* // idendities implied by gauge invariance and current conservation; checked analytically and numerically
     Complex c1 = qqbargLoops[0] + qqbargLoops[6] + qqbargLoops[7];
     Complex c2 = qqbargLoops[1] + qqbargLoops[8] + qqbargLoops[9];
     Complex c3 = qqbargLoops[3] + qqbargLoops[4];
     Complex c4 = qqbargLoops[2] + qqbargLoops[5] + qqbargLoops[10] + qqbargLoops[11];
     Complex c5 = 
     2. * qqbargLoops[3]/ik +
     2. * qqbargLoops[5]/jk +
     qqbargLoops[6] * (1.+ij/ik) +
     qqbargLoops[8] * (jk+ij)/ik +
     2. * qqbargLoops[10] * (1./ik+1./jk) +
     2. * qqbargLoops[12] * (1./ik+1./jk);
     Complex c6 = 
     2. * qqbargLoops[4]/jk +
     2. * qqbargLoops[5]/jk +
     qqbargLoops[7] * (ik+ij)/jk +
     qqbargLoops[9] * (1.+ij/jk) +
     2. * qqbargLoops[11] * (ik/jk2+1./jk);
     Complex c7 =
     0.5 * qqbargLoops[0] * (ij+ik) +
     0.5 * qqbargLoops[1] * (ij+jk) +
     qqbargLoops[2] * (1.+ik/jk) -
     qqbargLoops[12] * (1.+ik/jk);

     double x1 = c1 != 0. ? log(abs(real(c1 * conj(c1)))) : 0.;
     double x2 = c2 != 0. ? log(abs(real(c2 * conj(c2)))) : 0.;
     double x3 = c3 != 0. ? log(abs(real(c3 * conj(c3)))) : 0.;
     double x4 = c4 != 0. ? log(abs(real(c4 * conj(c4)))) : 0.;
     double x5 = c5 != 0. ? log(abs(real(c5 * conj(c5)))) : 0.;
     double x6 = c6 != 0. ? log(abs(real(c6 * conj(c6)))) : 0.;
     double x7 = c7 != 0. ? log(abs(real(c7 * conj(c7)))) : 0.;

     cerr << x1 << " " << x2 << " " << x3 << " " << x4 << " "
     << x5 << " " << x6 << " " << x7 << "\n";
  */

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralLeftLoopCurrent(const int i, const int,
								      const int j, const int,
								      const int k, const int gHel,
								      const int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kn = plusProduct(k,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kn = minusProduct(k,n);
  
  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);

  const LorentzVector<Complex> & minusC_nk = minusCurrent(n,k);
  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kn = minusCurrent(k,n);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      (sqrt(2) * c6 * plusP_jk * minusC_nk * minusP_ik)/(jk * minusP_kn) + 
      (sqrt(2) * c1 * plusP_jk * momentum(i) * minusP_in)/minusP_kn + 
      (sqrt(2) * c2 * plusP_jk * momentum(j) * minusP_in)/minusP_kn + 
      (2 * sqrt(2) * c3 * plusP_jk * momentum(k) * minusP_in)/(jk * minusP_kn) + 
      (sqrt(2) * c4 * plusP_ik * minusC_ij * minusP_in)/(ik * minusP_kn) - 
      (sqrt(2) * c7 * plusP_ik * plusP_jk * momentum(i) * minusP_ik * minusP_in)/(ik * minusP_kn) - 
      (sqrt(2) * c9 * plusP_ik * plusP_jk * momentum(j) * minusP_ik * minusP_in)/(ik * minusP_kn) - 
      (2 * sqrt(2) * c11 * plusP_ik * plusP_jk * momentum(k) * minusP_ik * minusP_in)/(ik * jk * minusP_kn) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ij * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c8 * sqr(plusP_jk) * momentum(i) * minusP_ik * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c10 * sqr(plusP_jk) * momentum(j) * minusP_ik * minusP_jn)/(jk * minusP_kn) - 
      (2 * sqrt(2) * c12 * sqr(plusP_jk) * momentum(k) * minusP_ik * minusP_jn)/(sqr(jk) * minusP_kn);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c1 * plusP_jn * momentum(i) * minusP_ik)/plusP_kn) - 
      (sqrt(2) * c2 * plusP_jn * momentum(j) * minusP_ik)/plusP_kn - 
      (2 * sqrt(2) * c3 * plusP_jn * momentum(k) * minusP_ik)/(jk * plusP_kn) - 
      (sqrt(2) * c4 * plusP_in * minusC_ij * minusP_ik)/(ik * plusP_kn) + 
      (sqrt(2) * c13 * minusC_kj * minusP_ik)/ik + (sqrt(2) * c13 * minusC_kj * minusP_ik)/jk - 
      (sqrt(2) * c6 * plusP_jk * minusC_kn * minusP_ik)/(jk * plusP_kn) + 
      (sqrt(2) * c7 * plusP_in * plusP_jk * momentum(i) * sqr(minusP_ik))/(ik * plusP_kn) + 
      (sqrt(2) * c9 * plusP_in * plusP_jk * momentum(j) * sqr(minusP_ik))/(ik * plusP_kn) + 
      (2 * sqrt(2) * c11 * plusP_in * plusP_jk * momentum(k) * sqr(minusP_ik))/(ik * jk * plusP_kn) - 
      (sqrt(2) * c5 * plusP_jn * minusC_ij * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c8 * plusP_jk * plusP_jn * momentum(i) * minusP_ik * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c10 * plusP_jk * plusP_jn * momentum(j) * minusP_ik * minusP_jk)/(jk * plusP_kn) + 
      (2 * sqrt(2) * c12 * plusP_jk * plusP_jn * momentum(k) * minusP_ik * minusP_jk)/(sqr(jk) * plusP_kn);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedLeftLoopCurrent(const int i, const int,
								    const int j, const int,
								    const int k, const int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_jk = plusProduct(j,k);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c6 * plusP_jk * minusC_ik)/jk) 
      - (sqrt(2) * c8 * sqr(plusP_jk) * momentum(i) * minusP_ij)/jk - 
      (sqrt(2) * c10 * sqr(plusP_jk) * momentum(j) * minusP_ij)/jk - 
      (2 * sqrt(2) * c12 * sqr(plusP_jk) * momentum(k) * minusP_ij)/sqr(jk) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ij * minusP_ij)/(jk * minusP_ik);
  }

  if ( gHel == -1 ) {
    return
      (sqrt(2) * c4 * plusP_ij * minusC_ij * minusP_ik)/(ik * plusP_jk) + 
      (sqrt(2) * c13 * minusC_kj * minusP_ik)/ik + (sqrt(2) * c13 * minusC_kj * minusP_ik)/jk + 
      (sqrt(2) * c6 * minusC_kj * minusP_ik)/jk - (sqrt(2) * c7 * plusP_ij * momentum(i)*
						   sqr(minusP_ik))/ik - 
      (sqrt(2) * c9 * plusP_ij * momentum(j) * sqr(minusP_ik))/ik - 
      (2 * sqrt(2) * c11 * plusP_ij * momentum(k) * sqr(minusP_ik))/(ik * jk);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralRightLoopCurrent(const int i,    const int,
								       const int j,    const int,
								       const int k,    const int gHel,
								       const int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kn = plusProduct(k,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kn = minusProduct(k,n);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);

  const LorentzVector<Complex> & minusC_nk = minusCurrent(n,k);
  const LorentzVector<Complex> & minusC_kn = minusCurrent(k,n);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c13 * plusP_ik * minusC_jk)/ik) - (sqrt(2) * c13 * plusP_ik * minusC_jk)/jk + 
      (sqrt(2) * c4 * plusP_ik * minusC_ji * minusP_in)/(ik * minusP_kn) + 
      (sqrt(2) * c6 * plusP_ik * minusC_nk * minusP_jk)/(jk * minusP_kn) - 
      (sqrt(2) * c7 * sqr(plusP_ik) * momentum(i) * minusP_in * minusP_jk)/(ik * minusP_kn) - 
      (sqrt(2) * c9 * sqr(plusP_ik) * momentum(j) * minusP_in * minusP_jk)/(ik * minusP_kn) - 
      (2 * sqrt(2) * c11 * sqr(plusP_ik) * momentum(k) * minusP_in * minusP_jk)/(ik * jk * minusP_kn) + 
      (sqrt(2) * c1 * plusP_ik * momentum(i) * minusP_jn)/minusP_kn + 
      (sqrt(2) * c2 * plusP_ik * momentum(j) * minusP_jn)/minusP_kn + 
      (2 * sqrt(2) * c3 * plusP_ik * momentum(k) * minusP_jn)/(jk * minusP_kn) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ji * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c8 * plusP_ik * plusP_jk * momentum(i) * minusP_jk * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c10 * plusP_ik * plusP_jk * momentum(j) * minusP_jk * minusP_jn)/(jk * minusP_kn) - 
      (2 * sqrt(2) * c12 * plusP_ik * plusP_jk * momentum(k) * minusP_jk * minusP_jn)/(sqr(jk) * minusP_kn);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c4 * plusP_in * minusC_ji * minusP_ik)/(ik * plusP_kn)) - 
      (sqrt(2) * c1 * plusP_in * momentum(i) * minusP_jk)/plusP_kn - 
      (sqrt(2) * c2 * plusP_in * momentum(j) * minusP_jk)/plusP_kn - 
      (2 * sqrt(2) * c3 * plusP_in * momentum(k) * minusP_jk)/(jk * plusP_kn) - 
      (sqrt(2) * c5 * plusP_jn * minusC_ji * minusP_jk)/(jk * plusP_kn) - 
      (sqrt(2) * c6 * plusP_ik * minusC_kn * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c7 * plusP_ik * plusP_in * momentum(i) * minusP_ik * minusP_jk)/(ik * plusP_kn) + 
      (sqrt(2) * c9 * plusP_ik * plusP_in * momentum(j) * minusP_ik * minusP_jk)/(ik * plusP_kn) + 
      (2 * sqrt(2) * c11 * plusP_ik * plusP_in * momentum(k) * minusP_ik * minusP_jk)/
      (ik * jk * plusP_kn) + 
      (sqrt(2) * c8 * plusP_ik * plusP_jn * momentum(i) * sqr(minusP_jk))/(jk * plusP_kn) + 
      (sqrt(2) * c10 * plusP_ik * plusP_jn * momentum(j) * sqr(minusP_jk))/(jk * plusP_kn) + 
      (2 * sqrt(2) * c12 * plusP_ik * plusP_jn * momentum(k) * sqr(minusP_jk))/(sqr(jk) * plusP_kn);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedRightLoopCurrent(const int i, const int,
								     const int j, const int,
								     const int k, const int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  
  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_jk = minusProduct(j,k);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c13 * plusP_ik * minusC_jk)/ik) - 
      (sqrt(2) * c13 * plusP_ik * minusC_jk)/jk - 
      (sqrt(2) * c6 * plusP_ik * minusC_jk)/jk + 
      (sqrt(2) * c7 * sqr(plusP_ik) * momentum(i) * minusP_ij)/ik + 
      (sqrt(2) * c9 * sqr(plusP_ik) * momentum(j) * minusP_ij)/ik + 
      (2 * sqrt(2) * c11 * sqr(plusP_ik) * momentum(k) * minusP_ij)/(ik * jk) - 
      (sqrt(2) * c4 * plusP_ik * minusC_ji * minusP_ij)/(ik * minusP_jk);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c5 * plusP_ij * minusC_ji * minusP_jk)/(jk * plusP_ik)) + 
      (sqrt(2) * c6 * minusC_ki * minusP_jk)/jk + 
      (sqrt(2) * c8 * plusP_ij * momentum(i) * sqr(minusP_jk))/jk + 
      (sqrt(2) * c10 * plusP_ij * momentum(j) * sqr(minusP_jk))/jk + 
      (2 * sqrt(2) * c12 * plusP_ij * momentum(k) * sqr(minusP_jk))/sqr(jk);
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftOneLoopCurrent(const int q,    const int qHel,
									 const int qbar, const int qbarHel,
									 const int g1,   const int g1Hel) {
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<2>(1,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5) * qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbargLeftLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5) * qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightOneLoopCurrent(const int q,    const int qHel,
									  const int qbar, const int qbarHel,
									  const int g1,   const int g1Hel) {

  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<2>(2,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5) * qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbargRightLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5) * qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

#ifdef CHECK_MatchboxCurrents

map<string,ofstream * >& MatchboxCurrents::checkStreams() {
  static map<string,ofstream * > theMap;
  return theMap;
}

ostream& MatchboxCurrents::checkStream(const string& id) {
  map<string,ofstream * >::iterator ret = checkStreams().find(id);
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
  double ac = abs(real(conj(c) * c));

  if ( ! isfinite(ac) ) {
    cerr << "ooops ... nan encountered in current conservation\n" << flush;
    return;
  }

  checkStream(id) << (ac > 0. ? log10(ac) : 0.) << "\n";

}

#endif // CHECK_MatchboxCurrents
