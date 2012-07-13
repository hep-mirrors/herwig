// -*- C++ -*-
//
// MatchboxCurrents.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxCurrents.h"

using namespace Herwig;

inline Complex csqr(Complex a) {
  return a*a;
}

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

  if ( getCurrent(hash<0>(3,1,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(Complex(0.,1.)*plusCurrent(q,qbar));
  }
  return cachedCurrent();
    
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightCurrent(int q,    int qHel,
								  int qbar, int qbarHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<0>(4,1,q,qHel,qbar,qbarHel)) ) {
    cacheCurrent(Complex(0.,1.)*plusCurrent(qbar,q));
  }

  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftCurrent(int q,    int qHel,
								  int qbar, int qbarHel,
								  int g,    int gHel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( gHel == 1 ) {
    if ( getCurrent(hash<0>(3,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*plusCurrent(q,qbar)
		    +(1./plusProduct(g,qbar))*plusCurrent(q,g)));
    }
    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<0>(3,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*plusCurrent(q,qbar)
		    +(1./minusProduct(q,g))*plusCurrent(g,qbar)));
    }
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
    if ( getCurrent(hash<0>(4,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,1.)*sqrt(2.)*
		   ((plusProduct(q,qbar)/(plusProduct(q,g)*plusProduct(g,qbar)))*plusCurrent(qbar,q)
		    +(1./plusProduct(q,g))*plusCurrent(qbar,g)));
    }
    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<0>(4,1,q,qHel,qbar,qbarHel,g,gHel)) ) {
      cacheCurrent(Complex(0.,-1.)*sqrt(2.)*
		   ((minusProduct(q,qbar)/(minusProduct(q,g)*minusProduct(g,qbar)))*plusCurrent(qbar,q)
		    +(1./minusProduct(g,qbar))*plusCurrent(g,q)));
    }
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggLeftCurrent(int q,    int qHel,
								   int qbar, int qbarHel,
								   int g1,   int g1Hel,
								   int g2,   int g2Hel) {
  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  int i = q; int j = qbar; int k = g1; int l = g2;

  if ( g1Hel == 1 && g2Hel == 1 ) {
    if ( getCurrent(hash<0>(5,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent((Complex(0,-2)*((invariant(i,l)*plusCurrent(i,j))/
				   (invariant(i,k) + invariant(i,l) + invariant(k,l)) - 
				   (minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,l))/
				   (invariant(i,k) + invariant(i,l) + invariant(k,l)) + 
				   (invariant(i,k)*(invariant(j,k) + invariant(k,l))*
				    (invariant(j,l)*plusCurrent(i,j) - 
				     minusProduct(j,l)*plusCurrent(i,k)*plusProduct(k,l)) + 
				    (invariant(j,k) + invariant(j,l) + invariant(k,l))*
				    minusProduct(i,k)*minusProduct(j,l)*plusProduct(i,l)*
				    (-(plusCurrent(i,j)*plusProduct(j,k)) + 
				     plusCurrent(i,l)*plusProduct(k,l)))/
				   (invariant(i,k)*invariant(j,l)*
				    (invariant(j,k) + invariant(j,l) + invariant(k,l)))))/
		   csqr(plusProduct(k,l)));
    }
    return cachedCurrent();
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    if ( getCurrent(hash<0>(5,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,1)*((invariant(i,k)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(i,l)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(j,k)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 + (invariant(j,l)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) + 
				 (minusProduct(k,l)*plusCurrent(l,j)*plusProduct(i,k))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
				 (minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) - 
				 (2.*csqr(minusProduct(j,k))*plusCurrent(i,k)*plusProduct(j,l))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,l)) - 
				 (2.*csqr(minusProduct(j,k))*plusCurrent(i,j)*
				  csqr(plusProduct(j,l)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))
				  *minusProduct(k,l)*plusProduct(k,l)) + 
				 (2.*minusProduct(i,k)*plusProduct(i,l)*
				  (invariant(j,l)*(-(minusProduct(i,k)*plusCurrent(i,j)) + 
						   minusProduct(k,l)*plusCurrent(l,j))*plusProduct(i,l) + 
				   (invariant(i,k) + invariant(i,l) + invariant(k,l))*
				   minusProduct(j,k)*plusCurrent(i,j)*plusProduct(j,l)))/
				 (invariant(i,k)*invariant(j,l)*
				  (invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,l)*plusProduct(k,l)) - 
				 (minusProduct(j,l)*plusCurrent(i,k)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (minusProduct(j,k)*plusCurrent(i,l)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l)))));
    }
    return cachedCurrent();
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    if ( getCurrent(hash<0>(5,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,1)*((invariant(i,k)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(i,l)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(j,k)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 + (invariant(j,l)*plusCurrent(i,j))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 + (minusProduct(k,l)*plusCurrent(l,j)*plusProduct(i,k))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
				 (minusProduct(k,l)*plusCurrent(k,j)*plusProduct(i,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) +
				 (2.*csqr(minusProduct(j,l))*plusCurrent(i,l)*plusProduct(j,k))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,l)) - 
				 (2.*csqr(minusProduct(j,l))*plusCurrent(i,j)*
				  csqr(plusProduct(j,k)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + 
						  invariant(k,l))*minusProduct(k,l)*plusProduct(k,l)) - 
				 (minusProduct(j,l)*plusCurrent(i,k)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 - (minusProduct(j,k)*plusCurrent(i,l)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (2.*plusProduct(i,k)*((minusProduct(i,l)*plusCurrent(i,j) + 
						       minusProduct(k,l)*plusCurrent(k,j))*
						      (invariant(j,l)*minusProduct(i,l)*plusProduct(i,k) - 
						       (invariant(i,k) + invariant(i,l) + invariant(k,l))*
						       minusProduct(j,l)*plusProduct(j,k)) + 
						      (invariant(i,k) + invariant(i,l) + invariant(k,l))*
						      minusProduct(j,l)*(minusProduct(i,l)*plusCurrent(i,l) + 
									 minusProduct(k,l)*plusCurrent(k,l))*
						      plusProduct(k,l)))/(invariant(i,k)*invariant(j,l)*
				  (invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,l)*plusProduct(k,l))));
    }
    return cachedCurrent();
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    if ( getCurrent(hash<0>(5,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,-2)*((invariant(j,k)*plusCurrent(i,j))/
				  ((invariant(j,k) + invariant(j,l) + invariant(k,l))*
				   csqr(minusProduct(k,l))) + 
				  (invariant(i,l)*(invariant(i,k)*plusCurrent(i,j) + 
						   minusProduct(k,l)*plusCurrent(l,j)*plusProduct(i,k)))/
				  (invariant(i,k)*
				   (invariant(i,k) + invariant(i,l) + invariant(k,l))*csqr(minusProduct(k,l))) - 
				  (minusProduct(i,l)*minusProduct(j,k)*plusCurrent(i,j)*plusProduct(i,k)*plusProduct(j,l))/
				  (invariant(i,k)*invariant(j,l)*csqr(minusProduct(k,l))) - 
				  (minusProduct(j,k)*plusCurrent(k,j)*plusProduct(i,k)*plusProduct(j,l))/
				  (invariant(i,k)*invariant(j,l)*minusProduct(k,l)) - 
				  (plusCurrent(i,j)*plusProduct(k,l))/
				  ((invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,l)) 
				  + (minusProduct(j,k)*plusCurrent(i,l)*plusProduct(k,l))/
				  ((invariant(j,k) + invariant(j,l) + invariant(k,l))*csqr(minusProduct(k,l))) - 
				  (plusCurrent(l,j)*plusProduct(i,k)*plusProduct(k,l))/
				  (invariant(i,k)*(invariant(i,k) + invariant(i,l) + invariant(k,l)))));
    }
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggRightCurrent(int q,    int qHel,
								    int qbar, int qbarHel,
								    int g1,   int g1Hel,
								    int g2,   int g2Hel) {

  static LorentzVector<Complex> czero(0.,0.,0.,0.);
  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  int i = q; int j = qbar; int k = g1; int l = g2;

  if ( g1Hel == 1 && g2Hel == 1 ) {
    if ( getCurrent(hash<0>(6,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,2)*(-(((invariant(i,l)*plusCurrent(j,i))/
				    (invariant(i,k) + invariant(i,l) + invariant(k,l)) + 
				    (invariant(k,l)*plusCurrent(j,i))/
				    (invariant(i,k) + invariant(i,l) + invariant(k,l)) + 
				    (invariant(j,k)*plusCurrent(j,i) + 
				     minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,k))/
				    (invariant(j,k) + invariant(j,l) + invariant(k,l)))/
				   csqr(plusProduct(k,l))) + 
				 (minusProduct(i,k)*((minusProduct(k,l)*plusCurrent(j,l))/
						     (invariant(i,k) + invariant(i,l) + invariant(k,l)) + 
						     (-((invariant(i,l)*plusCurrent(j,l)*plusProduct(k,l))/
							(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
						      (minusProduct(j,l)*plusProduct(j,k)*
						       (plusCurrent(j,i)*plusProduct(i,l) + 
							plusCurrent(j,k)*plusProduct(k,l)))/invariant(j,l))/
						     csqr(plusProduct(k,l))))/invariant(i,k)));
    }
    return cachedCurrent();
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    if ( getCurrent(hash<0>(6,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,1)*((invariant(i,k)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(i,l)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(j,k)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 + (invariant(j,l)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 - (minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,k))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (minusProduct(k,l)*plusCurrent(k,i)*plusProduct(j,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 - (2.*csqr(minusProduct(j,k))*plusCurrent(j,i)*
				    csqr(plusProduct(j,l)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))
				  *minusProduct(k,l)*plusProduct(k,l)) + 
				 (2.*minusProduct(j,k)*plusCurrent(l,i)*csqr(plusProduct(j,l)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,l)) + 
				 (minusProduct(i,l)*plusCurrent(j,k)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
				 (minusProduct(i,k)*plusCurrent(j,l)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) - 
				 (2.*minusProduct(i,k)*(invariant(j,l)*minusProduct(i,k)*plusProduct(i,l)*
						       (plusCurrent(j,i)*plusProduct(i,l) + plusCurrent(j,k)*plusProduct(k,l)) 
						       - (invariant(i,k) + invariant(i,l) + invariant(k,l))*plusProduct(j,l)*
						       (minusProduct(j,k)*(plusCurrent(j,i)*plusProduct(i,l) + 
									   plusCurrent(j,k)*plusProduct(k,l)) - 
							minusProduct(k,l)*(plusCurrent(l,i)*plusProduct(i,l) + 
									   plusCurrent(l,k)*plusProduct(k,l)))))/
				 (invariant(i,k)*invariant(j,l)*
				  (invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(k,l)*plusProduct(k,l))));
    }
    return cachedCurrent();
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    if ( getCurrent(hash<0>(6,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,1)*((invariant(i,k)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(i,l)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) 
				 - (invariant(j,k)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) 
				 + (invariant(j,l)*plusCurrent(j,i))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (minusProduct(k,l)*plusCurrent(l,i)*plusProduct(j,k))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (minusProduct(k,l)*plusCurrent(k,i)*plusProduct(j,l))/
				 (invariant(k,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))) - 
				 (2.*csqr(minusProduct(j,l))*plusCurrent(j,i)*csqr(plusProduct(j,k)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*
				  minusProduct(k,l)*plusProduct(k,l)) - 
				 (2.*minusProduct(j,l)*plusCurrent(k,i)*csqr(plusProduct(j,k)))/
				 (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*plusProduct(k,l)) + 
				 (minusProduct(i,l)*plusCurrent(j,k)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
				 (minusProduct(i,k)*plusCurrent(j,l)*plusProduct(k,l))/
				 (invariant(k,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))) + 
				 (2.*minusProduct(i,l)*plusProduct(i,k)*
				  ((invariant(i,k) + invariant(i,l) + invariant(k,l))*minusProduct(j,l)*
				   plusCurrent(j,i)*plusProduct(j,k) + invariant(j,l)*minusProduct(i,l)*
				   (-(plusCurrent(j,i)*plusProduct(i,k)) + plusCurrent(j,l)*plusProduct(k,l))))/
				 (invariant(i,k)*invariant(j,l)*(invariant(i,k) + invariant(i,l) + invariant(k,l))*
				  minusProduct(k,l)*plusProduct(k,l))));
    }
    return cachedCurrent();
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    if ( getCurrent(hash<0>(6,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
      cacheCurrent(Complex(0,-2)*((invariant(i,l)*plusCurrent(j,i))/
				  ((invariant(i,k) + invariant(i,l) + invariant(k,l))*
				   csqr(minusProduct(k,l))) - 
				  (minusProduct(i,l)*minusProduct(j,k)*plusCurrent(j,i)*plusProduct(i,k)*plusProduct(j,l))/
				  (invariant(i,k)*invariant(j,l)*csqr(minusProduct(k,l))) + 
				  (minusProduct(i,l)*plusCurrent(l,i)*plusProduct(i,k)*plusProduct(j,l))/
				  (invariant(i,k)*invariant(j,l)*minusProduct(k,l)) + 
				  (invariant(j,k)*(invariant(j,l)*plusCurrent(j,i) - 
						   minusProduct(k,l)*plusCurrent(k,i)*plusProduct(j,l)))/
				  (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l))*
				   csqr(minusProduct(k,l))) - 
				  (plusCurrent(j,i)*plusProduct(k,l))/
				  ((invariant(j,k) + invariant(j,l) + invariant(k,l))*minusProduct(k,l)) - 
				  (minusProduct(i,l)*plusCurrent(j,k)*plusProduct(k,l))/
				  ((invariant(i,k) + invariant(i,l) + invariant(k,l))*
				   csqr(minusProduct(k,l))) + 
				  (plusCurrent(k,i)*plusProduct(j,l)*plusProduct(k,l))/
				  (invariant(j,l)*(invariant(j,k) + invariant(j,l) + invariant(k,l)))));
    }
    return cachedCurrent();
  }

  return czero;

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
    if ( getCurrent(hash<0>(7,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(k,i)*minusProduct(i,l)*plusCurrent(i,j)+
		     plusProduct(i,k)*minusProduct(l,k)*plusCurrent(k,j))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,k)*minusProduct(l,j)*plusCurrent(i,j)+
		     plusProduct(l,k)*minusProduct(l,j)*plusCurrent(i,l))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<0>(7,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(l,i)*minusProduct(i,k)*plusCurrent(i,j)+
		     plusProduct(i,l)*minusProduct(k,l)*plusCurrent(l,j))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,l)*minusProduct(k,j)*plusCurrent(i,j)+
		     plusProduct(k,l)*minusProduct(k,j)*plusCurrent(i,k))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
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
    if ( getCurrent(hash<0>(8,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(k,i)*minusProduct(i,l)*plusCurrent(j,i)+
		     plusProduct(l,k)*minusProduct(l,i)*plusCurrent(j,l))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,k)*minusProduct(l,j)*plusCurrent(j,i)+
		     plusProduct(j,k)*minusProduct(l,k)*plusCurrent(k,i))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<0>(8,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/invariant(k,l))*
		   ((plusProduct(l,i)*minusProduct(i,k)*plusCurrent(j,i)+
		     plusProduct(k,l)*minusProduct(k,i)*plusCurrent(j,k))/
		    (invariant(k,l)+invariant(i,l)+invariant(i,k))-
		    (plusProduct(j,l)*minusProduct(k,j)*plusCurrent(j,i)+
		     plusProduct(j,l)*minusProduct(k,l)*plusCurrent(l,i))/
		    (invariant(k,l)+invariant(j,l)+invariant(j,k))));
    }
    return cachedCurrent();
  }

  return czero;

}
