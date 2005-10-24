// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVLoopVertex class.
//

#include "SVVLoopVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
  namespace Helicity {
    using namespace ThePEG;
    using namespace Looptools;
    
    SVVLoopVertex::~SVVLoopVertex() {}
    
    void SVVLoopVertex::persistentOutput(PersistentOStream & os) const {}
    
    void SVVLoopVertex::persistentInput(PersistentIStream & is, int) {} 
   
    ClassDescription<SVVLoopVertex> SVVLoopVertex::initSVVLoopVertex;
    // Definition of the static class description member.
    
    void SVVLoopVertex::Init() {
      
      static ClassDocumentation<SVVLoopVertex> documentation
	("The SVVLoopVertex class calculates the tenosr integral"
	 " coefficients using Looptools.");
     }
  
    void SVVLoopVertex::setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
				    tcPDPtr part3) {
      //Kinematic invariants
      Energy2 ps2 = invariant(0,0);
      Energy2 pv1s = invariant(1,1);
      Energy2 pv2s = invariant(2,2);
      ffini();
      Complex a(0.),b(0.),c(0.),d(0.),e(0.),f(0.);
      for(unsigned int i=0;i<masses.size();++i) {
	Energy2 mls = masses[i]*masses[i];
	if(type[i] == PDT::Spin1Half) {
	  Complex C0 = C0i(cc0,pv1s,pv2s,ps2,mls,mls,mls); 
	  int theC = Cget(ps2,pv2s,pv1s,
			mls,mls,mls);
	  Complex C1 = Cval(cc1,theC);Complex C2 = Cval(cc2,theC);
	  Complex C00 = Cval(cc00,theC);Complex C11 = Cval(cc11,theC);
	  Complex C12 = Cval(cc12,theC);Complex C21 = C12;
	  Complex C22 = Cval(cc22,theC);
	  a += 2.*(left[i]+right[i])*masses[i]*(-2.*B0(ps2,mls,mls)+ 
						C0*(pv1s + pv2s - ps2)
						+ 8.*C00);
	  b += 8.*(left[i]+right[i])*masses[i]*(C0 + 3.*C1 +3.*C2 +
						2.*(C11 + 2.*C12 + C22));
	  c += 4.*(left[i]+right[i])*masses[i]*(C0 +2.*(2.*C1+C2 + 
							2.*(C11 +C12)));
	  d += 4.*(left[i]+right[i])*masses[i]*(C0 + 4.*(C1+C11+C12));
	  e += 8.*(left[i]+right[i])*masses[i]*(C1 + 2.*C11);
	  f += 4.*(left[i]-right[i])*masses[i]*C0;
	}
 	
	else if(type[i] == PDT::Spin1) {
	  int theC = Cget(ps2,pv2s,pv1s,mls,mls,mls);
	  Complex C1 = Cval(cc1,theC);Complex C2 = Cval(cc2,theC);
	  Complex C00 = Cval(cc00,theC);Complex C11 = Cval(cc11,theC);
	  Complex C12 = Cval(cc12,theC);Complex C21 = C12;
	  Complex C22 = Cval(cc22,theC);
	  	  
	  /**
	   * vector type can contain different types of particle 
	   * and hence the coupling is different
	   * Here left[i] is used for the coupling of the ith 
	   * type rather than creating another
	   * vector to hold them.
	   */
	 Energy2  pv12 = pv1s*pv2s;
	 
	  a += 0.5*left[i]*(B0(ps2,mls,mls)*(2.*(mls*mls)*(-6.*mls + pv1s + pv2s) + mls*(-2.*mls + pv1s + pv2s)*ps2) + 2.*(8.*pow(mls,3)*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv1s - 2.*(mls*mls)*C0(pv1s,pv2s,ps2,mls,mls,mls)*(pv1s*pv1s) + 2.*mls*B00(pv1s,mls,mls)*pv2s + 8.*pow(mls,3)*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv2s + mls*B0(pv1s,mls,mls)*pv12 + mls*B0(pv2s,mls,mls)*pv12 - B00(pv1s,mls,mls)*pv12 - 2.*(mls*mls)*C0(pv1s,pv2s,ps2,mls,mls,mls)*(pv2s*pv2s) + 
          B00(pv2s,mls,mls)*pv1s*(2.*mls - pv2s) - mls*A0(mls)*(pv1s + pv2s) - 8.*pow(mls,3)*C0(pv1s,pv2s,ps2,mls,mls,mls)*ps2 + 
          2.*(mls*mls)*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv1s*ps2 + 2.*(mls*mls)*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv2s*ps2 - 
          mls*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv12*ps2 + (24.*pow(mls,3) + 2.*mls*pv12 - 4.*(mls*mls)*(pv1s + pv2s) + 
							(2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C00 ))/pow(mls,3);

	  b +=  -0.25*left[i]*(8.*mls*B0(pv2s,mls,mls)*pv2s - 4.*B11(pv1s,mls,mls)*(2.*mls - pv1s)*pv2s + 
			       2.*B1(pv1s,mls,mls)*(2.*mls - pv1s)*(4.*mls - 3.*pv2s) + 4.*B00(pv2s,mls,mls)*(2.*mls - pv2s) - 2.*A0(mls)*(2.*mls + pv2s) + 
			       2.*B0(pv1s,mls,mls)*(4.*(mls*mls) - 2.*mls*pv1s + pv12) - 2.*mls*B0(ps2,mls,mls)*(pv2s - ps2) - 4.*mls*B11(ps2,mls,mls)*(2.*mls + ps2) + 
       B1(ps2,mls,mls)*(2.*mls - pv2s)*(2.*mls + ps2) - B1(ps2,mls,mls)*(6.*mls - pv2s)*(2.*mls + ps2) - B0(ps2,mls,mls)*(8.*(mls*mls) + (2.*mls - pv2s)*ps2) - 2.*C0(pv1s,pv2s,ps2,mls,mls,mls)*(2.*mls*(12.*(mls*mls) + pv2s*(pv1s + pv2s) - 2.*mls*(pv1s + 3.*pv2s)) + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2) + 4.*mls*(2.*mls*B0(ps2,mls,mls) + 
          (B1(ps2,mls,mls) + B11(ps2,mls,mls))*(2.*mls + ps2)) - 
       2.*(2.*mls*(36.*(mls*mls) + 2.*pv12 + (pv2s*pv2s) - 6.*mls*(pv1s + pv2s)) + 
          (12.*(mls*mls) + 3.*pv12 - 2.*mls*(3.*pv1s + 4.*pv2s))*ps2)*C1 - 2.*(2.*mls*(36.*(mls*mls) + 2.*pv12 + (pv2s*pv2s) - 6.*mls*(pv1s + pv2s)) + 
          (12.*(mls*mls) + 3.*pv12 - 2.*mls*(3.*pv1s + 4.*pv2s))*ps2)*C2 - 4.*(24.*pow(mls,3) + 2.*mls*pv12 - 4.*(mls*mls)*(pv1s + pv2s) + 
          (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C11 - 8.*(24.*pow(mls,3) + 2.*mls*pv12 - 4.*(mls*mls)*(pv1s + pv2s) + 
          (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C12 - 4.*(24.*pow(mls,3) + 2.*mls*pv12 - 4.*(mls*mls)*(pv1s + pv2s) + 
							 (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C22)/pow(mls,3);

	  c+= (-1/24.)*left[i]*(-12*(mls*mls) + 8*mls*pv1s - (pv1s*pv1s) + 48*B00(pv1s,mls,mls)*(2*mls - pv1s) + 
       24*B00(pv2s,mls,mls)*(2*mls - pv2s) - 6*A0(mls)*(6*mls + pv1s - 2*ps2) - 6*B1(pv1s,mls,mls)*(2*mls - pv1s)*(2*mls - pv1s + 3*pv2s - ps2) - 
       24*mls*B11(ps2,mls,mls)*(2*mls + ps2) - 3*B1(ps2,mls,mls)*(4*mls - pv1s + pv2s - ps2)*(2*mls + ps2) - 
       3*B1(ps2,mls,mls)*(2*mls + ps2)*(4*mls + pv1s - pv2s + ps2) + 6*B1(pv2s,mls,mls)*(2*mls - pv2s)*(-3*pv1s + pv2s + ps2) + 
       12*B11(pv2s,mls,mls)*(2*mls - pv2s)*(-pv1s + pv2s + ps2) + 6*B11(pv1s,mls,mls)*(2*mls - pv1s)*(3*pv1s - 2*pv2s + 2*ps2) + 
       6*B0(pv2s,mls,mls)*(2*mls*pv1s + (4*mls + pv1s)*pv2s - 4*mls*ps2) - 6*B0(pv1s,mls,mls)*(2*(mls*mls) - 5*mls*pv1s - (2*mls + pv1s)*pv2s + 
          4*mls*ps2) - 3*B0(ps2,mls,mls)*(2*mls*(4*mls + pv1s + 2*pv2s) - (4*mls + pv1s)*ps2) - 3*B0(ps2,mls,mls)*(2*mls*(4*mls + 2*pv1s + pv2s) - (2*mls + pv2s)*ps2 + 
          (ps2*ps2)) - 6*C0(pv1s,pv2s,ps2,mls,mls,mls)*(2*mls*(12*(mls*mls) - (pv1s*pv1s) + pv12 + (pv2s*pv2s) - 6*mls*(pv1s + pv2s)) + 
          (12*(mls*mls) + pv12 + 2*mls*(pv1s - pv2s))*ps2 - 2*mls*(ps2*ps2)) + 24*mls*(2*mls*B0(ps2,mls,mls) + 
          (B1(ps2,mls,mls) + B11(ps2,mls,mls))*(2*mls + ps2)) - 24*(mls*(24*(mls*mls) - (pv1s*pv1s) + 2*pv12 + (pv2s*pv2s) - 4*mls*(pv1s + pv2s)) + 
          (4*(mls*mls) + pv12 - mls*(pv1s + 3*pv2s))*ps2)*C1 - 12*(24*pow(mls,3) - 2*mls*(pv1s*pv1s) - 4*(mls*mls)*(pv1s + pv2s) + 
          (4*(mls*mls) - 2*mls*pv2s + pv12)*ps2)*C2 - 24*(24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s) + 
          (2*mls - pv1s)*(2*mls - pv2s)*ps2)* C11 -  24*(24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s) + 
          (2*mls - pv1s)*(2*mls - pv2s)*ps2)*C12 + 72*(2*mls - pv1s)*(-0.125*mls - 0.125*A0(mls) + (0.25*mls*B1(pv1s,mls,mls)) + (0.125*pv1s*B11(pv1s,mls,mls)))  - 
			12*(-2*mls + pv1s)*(0.25*mls - 0.25*A0(mls) - (0.5*mls*B1(pv1s,mls,mls)) - (0.75*pv1s*B11(pv1s,mls,mls))) )/pow(mls,3);

	  d+= -1.*left[i]*(-2*(mls*mls)*B0(ps2,mls,mls) - 2*mls*C0(pv1s,pv2s,ps2,mls,mls,mls)*(8*(mls*mls) + pv12 - 2*mls*(pv1s + pv2s)) - 
       mls*B1(ps2,mls,mls)*(2*mls + ps2) - mls*B11(ps2,mls,mls)*(2*mls + ps2) + mls*(2*mls*B0(ps2,mls,mls) + (B1(ps2,mls,mls) + B11(ps2,mls,mls))*(2*mls + ps2)) - (24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s) +  (2*mls - pv1s)*(2*mls - pv2s)*ps2)* C1 -  2*mls*pv12*C2 - (24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s) + (2*mls - pv1s)*(2*mls - pv2s)*ps2)*C11 - (24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s)+(2*mls - pv1s)*(2*mls - pv2s)*ps2)*C12 )/pow(mls,3);

	  e+= -0.25*left[i]*(8*mls*B0(pv1s,mls,mls)*pv1s + 4*B00(pv1s,mls,mls)*(2*mls - pv1s) - 
       2*A0(mls)*(2*mls + pv1s) - 4*B11(pv2s,mls,mls)*pv1s*(2*mls - pv2s) + 2*B1(pv2s,mls,mls)*(4*mls - 3*pv1s)*(2*mls - pv2s) + 
       2*B0(pv2s,mls,mls)*(4*(mls*mls) - 2*mls*pv2s + pv12) - 2*mls*B0(ps2,mls,mls)*(pv1s - ps2) + 
       4*mls*C0(pv1s,pv2s,ps2,mls,mls,mls)*pv1s*(4*mls - pv2s - ps2) - 4*mls*B11(ps2,mls,mls)*(2*mls + ps2) + 
       B1(ps2,mls,mls)*(2*mls - pv1s)*(2*mls + ps2) - B1(ps2,mls,mls)*(6*mls - pv1s)*(2*mls + ps2) - 
       B0(ps2,mls,mls)*(8*(mls*mls) + (2*mls - pv1s)*ps2) + 4*mls*(2*mls*B0(ps2,mls,mls) + 
          (B1(ps2,mls,mls) + B11(ps2,mls,mls))*(2*mls + ps2)) - 2*(2*mls*(12*(mls*mls) - pv1s + 2*pv12 - 2*mls*(pv1s + pv2s)) + 
          (4*(mls*mls) - 2*mls*pv2s + pv12)*ps2)*C1 - 4*(24*pow(mls,3) + 2*mls*pv12 - 4*(mls*mls)*(pv1s + pv2s) + 
          (2*mls - pv1s)*(2*mls - pv2s)*ps2)*C11 )/pow(mls,3);
	}    
	else if(type[i] == PDT::Spin0) {
	  int theC = Cget(ps2,pv2s,pv1s,
			  mls,mls,mls);
	  Complex C1 = Cval(cc1,theC);Complex C2 = Cval(cc2,theC);
	  Complex C00 = Cval(cc00,theC);Complex C11 = Cval(cc11,theC);
	  Complex C12 = Cval(cc12,theC);Complex C21 = C12;
	  Complex C22 = Cval(cc22,theC);
	  
	  /**
	   * vector type can contain different types of particle 
	   * and hence the coupling is different
	   * Here left[i] is used for the coupling of the ith 
	   * type rather than creating another
	   * vector to hold them.
	  */
	  a += left[i]*(B0(ps2,mls,mls) - 4.*C00);

	  b += -2.*left[i]*(C0(pv1s,pv2s,ps2,mls,mls,mls)+ 3.*C1 + 3.*C2 +2.*(C11 + 2.*C12 + C22 ));

	  c += -1.*left[i]*(C0(pv1s,pv2s,ps2,mls,mls,mls)+ 2.*(2.*C1+ C2 + 2.*(C11 +C12)));

	  d+=-4.*left[i]*(C1 +C11 + C12);

	  e+=-2.*left[i]*(C1 + 2.*C11);
	}
	else {
	  throw HelicityConsistencyError() << "SVVLoopVertex::setCoupling "
					   << "Incorrect particle in SVV loop."
					   << Exception::warning;
	}
      }
      //Looptools defines integrals differently
      double fact = 1/16./Constants::pi/Constants::pi;
      a00(fact*a);
      a11(fact*b);
      a12(fact*c);
      a21(fact*d);
      a22(fact*e);
      aEp(fact*f);
    }
  } 
}
