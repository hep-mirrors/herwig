#include "ThePEG/Config/Constants.h"

namespace Herwig {
using namespace ThePEG;


  /** \ingroup Utilities
   * This is a namespace which provides the loop functions
   * from NPB297 (1988) 221-243 which are used to
   * calculate Higgs production with an additional jet
   * and the real correction to \f$h^0\to gg\f$.
   */
 namespace HiggsLoopFunctions {
   
   /**
    *  Epsilon parameter
    */
   const Complex epsi = Complex(0.,-1.e-20);
   
   /**
    *  The \f$W_1(s)\f$ function of NPB297 (1988) 221-243.
    * @param s   The invariant
    * @param mf2 The fermion mass squared
    */
   Complex W1(Energy2 s,Energy2 mf2) {
     double root = sqrt(abs(1.-4.*mf2/s));
     if(s<ZERO)     return 2.*root*asinh(0.5*sqrt(-s/mf2));
     else if(s<4.*mf2) return 2.*root*asin(0.5*sqrt( s/mf2));
     else              return root*(2.*acosh(0.5*sqrt(s/mf2))
				    -Constants::pi*Complex(0.,1.));
   }
   
   /**
    *  The \f$W_2(s)\f$ function of NPB297 (1988) 221-243.
    * @param s   The invariant
    * @param mf2 The fermion mass squared
    */
   Complex W2(Energy2 s,Energy2 mf2) {
     double root=0.5*sqrt(abs(s)/mf2);
     if(s<ZERO)     return 4.*sqr(asinh(root));
     else if(s<4.*mf2) return -4.*sqr(asin(root));
     else              return 4.*sqr(acosh(root))-sqr(Constants::pi)
			 -4.*Constants::pi*acosh(root)*Complex(0.,1.);
   }
   
   /**
    * The \f$I_3(s,t,u,v)\f$ function of NPB297 (1988) 221-243.
    * @param s The \f$s\f$ invariant
    * @param t The \f$t\f$ invariant
    * @param u The \f$u\f$ invariant
    * @param v The \f$v\f$ invariant
    * @param mf2 The fermion mass squared
    */
   Complex I3(Energy2 s, Energy2 t, Energy2 u, Energy2 v, Energy2 mf2) {
     double ratio=(4.*mf2*t/(u*s)),root(sqrt(1+ratio));
     if(v==ZERO) return 0.;
     Complex y=0.5*(1.+sqrt(1.-4.*(mf2+epsi*MeV*MeV)/v));
     Complex xp=0.5*(1.+root),xm=0.5*(1.-root);
     Complex output = 
       Math::Li2(xm/(xm-y))-Math::Li2(xp/(xp-y))+
       Math::Li2(xm/(y-xp))-Math::Li2(xp/(y-xm))+
       log(-xm/xp)*log(1.-epsi-v/mf2*xp*xm);
     return output*2./root;
   }
   
   /**
    * The \f$W_3(s,t,u,v)\f$ function of NPB297 (1988) 221-243.
    * @param s   The \f$s\f$ invariant
    * @param t   The \f$t\f$ invariant
    * @param u   The \f$u\f$ invariant
    * @param v   The \f$u\f$ invariant
    * @param mf2 The fermion mass squared.
    */
   Complex W3(Energy2 s, Energy2 t, Energy2 u, Energy2 v, Energy2 mf2) {
     return I3(s,t,u,v,mf2)-I3(s,t,u,s,mf2)-I3(s,t,u,u,mf2);
   }
   
   /**
    * The \f$b_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
    * @param s   The \f$s\f$ invariant
    * @param t   The \f$t\f$ invariant
    * @param u   The \f$u\f$ invariant
    * @param mf2 The fermion mass squared.
    */
   Complex b2(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) {
     Energy2 mh2(s+u+t);
     complex<Energy2> output=s*(u-s)/(s+u)+2.*u*t*(u+2.*s)/sqr(s+u)*(W1(t,mf2)-W1(mh2,mf2))
       +(mf2-0.25*s)*(0.5*(W2(s,mf2)+W2(mh2,mf2))-W2(t,mf2)+W3(s,t,u,mh2,mf2))
       +sqr(s)*(2.*mf2/sqr(s+u)-0.5/(s+u))*(W2(t,mf2)-W2(mh2,mf2))
       +0.5*u*t/s*(W2(mh2,mf2)-2.*W2(t,mf2))
       +0.125*(s-12.*mf2-4.*u*t/s)*W3(t,s,u,mh2,mf2);
     return output*mf2/sqr(mh2);
   }
   
   /**
    * The \f$b_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
    * @param s   The \f$s\f$ invariant
    * @param t   The \f$t\f$ invariant
    * @param u   The \f$u\f$ invariant
    * @param mf2 The fermion mass squared.
    */
   Complex b4(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) {
     Energy2 mh2(s+t+u);
     return mf2/mh2*(-2./3.+(mf2/mh2-0.25)*(W2(t,mf2)-W2(mh2,mf2)+W3(s,t,u,mh2,mf2)));
   }
   
   /**
    * The \f$A_1(m_h^2)\f$ function of NPB297 (1988) 221-243.
    * @param mh2 The Higgs mass squared
    * @param mf2 The fermion mass squared.
    */
   Complex A1(Energy2 mh2, Energy2 mf2) {
     return mf2/mh2*(4.-W2(mh2,mf2)*(1.-4.*mf2/mh2));
   }
   
   /**
    * The \f$A_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
    * @param s   The \f$s\f$ invariant
    * @param t   The \f$t\f$ invariant
    * @param u   The \f$u\f$ invariant
    * @param mf2 The fermion mass squared.
    */
   Complex A2(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) {
     return b2(s,t,u,mf2)+b2(s,u,t,mf2);
   }
   
   /**
    * The \f$A_4(s,t,u)\f$ function of NPB297 (1988) 221-243.
    * @param s   The \f$s\f$ invariant
    * @param t   The \f$t\f$ invariant
    * @param u   The \f$u\f$ invariant
    * @param mf2 The fermion mass squared.
    */
   Complex A4(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) {
     return b4(s,t,u,mf2)+b4(u,s,t,mf2)+b4(t,u,s,mf2);
   }
 }
}
