// test program to check looptool linking against c++ code.

#include <iostream>
#include <complex>

#include "Herwig/Looptools/clooptools.h"

namespace LT = Herwig::Looptools;
using namespace std;

int main() {

  { 
    double ps2 = 1.2268e+10;
    double pv1s = -4.76837e-07;
    double pv2s = -4.76837e-07;
    double mls = 2.65501e+11;
    
    LT::ltini();
     LT::clearcache();

    int theC = LT::Cget(ps2,pv2s,pv1s,
		     mls,mls,mls);
    cerr << "theC: " << theC << '\n';
    complex<double> C1 = LT::Cval(LT::cc1,theC);
    LT::ltexi();
    cout << "after: " << C1 << endl;
  
    cerr << "int : " << sizeof(int) 
	 << " long : " << sizeof(long) << '\n';
  }

  // ##################################################

  LT::ltini();

  typedef complex<double> Complex;

    //Kinematic invariants
  double ps2 = 1.2268e+10;
  double pv1s = -4.76837e-07;
  double pv2s = -4.76837e-07;
  //cerr << ps2 << ' ' << pv1s << ' ' << pv2s << ' ' << 6 << '\n';
  Complex a(0.),b(0.),c(0.),d(0.),e(0.),f(0.);
  for(unsigned int i = 0; i< 1; ++i) {
    double lmass = 515268;
    //cerr << lmass << '\n';
    double mls = lmass * lmass;
    Complex lc = 803.477;
    //cerr << lc << '\n'; 
    LT::clearcache();
    {
      Complex C0 =  LT::C0i(LT::cc0,pv1s,pv2s,ps2,mls,mls,mls); 
      int  theC = LT::Cget(ps2,pv2s,pv1s,    mls,mls,mls);
      Complex C1  = LT::Cval(LT::cc1,theC);
      Complex C2  = LT::Cval(LT::cc2,theC);
      Complex C00 = LT::Cval(LT::cc00,theC);
      Complex C11 = LT::Cval(LT::cc11,theC);
      Complex C12 = LT::Cval(LT::cc12,theC);
      Complex C22 = LT::Cval(LT::cc22,theC);
      Complex lpr = lc + 803.477;
      //cerr << lpr << '\n';

      a +=  4.*lpr*lmass*(-2.*LT::B0(ps2,mls,mls)
			  + C0*(pv1s + pv2s - ps2) 
			  + 8.*C00)/ps2;
      b +=  8.*lpr*lmass*(C0 + 3.*C1 +3.*C2 + 2.*(C11 + 2.*C12 + C22)); 
      c +=  4.*lpr*lmass*(C0 +2.*(2.*C1+C2 + 2.*(C11 +C12)));
      d +=  4.*lpr*lmass*(C0 + 4.*(C1+C11+C12));
      e +=  8.*lpr*lmass*(C1 + 2.*C11);
      f +=  4.*(lc - 803.477)*lmass*C0;


      cerr << a << ' ' 
	   << b << ' '
	   << c << ' '
	   << d << ' '
	   << e << ' '
	   << f << '\n';

    } LT::clearcache();
    {
      int theC = LT::Cget(ps2,pv2s,pv1s,mls,mls,mls);
      Complex C1 = LT::Cval(LT::cc1,theC);
      Complex C2 = LT::Cval(LT::cc2,theC);
      Complex C00 = LT::Cval(LT::cc00,theC);
      Complex C11 = LT::Cval(LT::cc11,theC);
      Complex C12 = LT::Cval(LT::cc12,theC);
      Complex C22 = LT::Cval(LT::cc22,theC);
      
      /**
       * vector type can contain different types of particle 
       * and hence the coupling is different
       * Here left is used for the coupling of the ith 
       * type rather than creating another
       * vector to hold them.
       */
      double pv12 = pv1s*pv2s;
      Complex 
	C0A(LT::C0(pv1s,pv2s,ps2,mls,mls,mls)),A0A(LT::A0(mls)),
	B0A(LT::B0(ps2 ,mls,mls)),
	B1A(LT::B1(ps2 ,mls,mls)),B11A(LT::B11(ps2 ,mls,mls)),
	B0B(LT::B0(pv1s,mls,mls)),B00B(LT::B00(pv1s,mls,mls)),
	B1B(LT::B1(pv1s,mls,mls)),B11B(LT::B11(pv1s,mls,mls)),
	B0C(LT::B0(pv2s,mls,mls)),B00C(LT::B00(pv2s,mls,mls)),
	B1C(LT::B1(pv2s,mls,mls)),B11C(LT::B11(pv2s,mls,mls));
      double mls2(mls*mls),mls3(mls2*mls);
      // coefficient
      a += 
	0.5*lc*(B0A*(2.*mls2*(-6.*mls + pv1s + pv2s) 
		       + mls*(-2.*mls + pv1s + pv2s)*ps2) 
	+ 2.*(8.*mls3*C0A*pv1s - 2.*mls2*C0A*(pv1s*pv1s) 
	      + 2.*mls*B00B*pv2s + 8.*mls3*C0A*pv2s 
	      + mls*B0B*pv12 + mls*B0C*pv12 - B00B*pv12 
	      - 2.*mls2*C0A*(pv2s*pv2s) 
	      + B00C*pv1s*(2.*mls - pv2s) 
	      - mls*A0A*(pv1s + pv2s) - 8.*mls3*C0A*ps2 
	      + 2.*mls2*C0A*pv1s*ps2 + 2.*mls2*C0A*pv2s*ps2 
	      - mls*C0A*pv12*ps2 
	      + (24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
		 + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C00 ) )/mls3*2./ps2;

      b += -0.25*lc*
	( 8.*mls*B0C*pv2s - 4.*B11B*(2.*mls - pv1s)*pv2s 
	  + 2.*B1B*(2.*mls - pv1s)*(4.*mls - 3.*pv2s) + 4.*B00C*(2.*mls - pv2s) 
	  - 2.*A0A*(2.*mls + pv2s) + 2.*B0B*(4.*mls2 - 2.*mls*pv1s + pv12) 
	  - 2.*mls*B0A*(pv2s - ps2) - 4.*mls*B11A*(2.*mls + ps2) 
	  + B1A*(2.*mls - pv2s)*(2.*mls + ps2) 
	  - B1A*(6.*mls - pv2s)*(2.*mls + ps2) 
	  - B0A*(8.*mls2 + (2.*mls - pv2s)*ps2) 
	  - 2.*C0A*(2.*mls*(12.*mls2 + pv2s*(pv1s + pv2s) - 2.*mls*(pv1s + 3.*pv2s)) 
		    + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2) 
	  + 4.*mls*(2.*mls*B0A + (B1A + B11A)*(2.*mls + ps2)) 
	  - 2.*(2.*mls*(36.*mls2 + 2.*pv12 + (pv2s*pv2s) - 6.*mls*(pv1s + pv2s)) 
		+ (12.*mls2 + 3.*pv12 - 2.*mls*(3.*pv1s + 4.*pv2s))*ps2)*C1 
	  - 2.*(2.*mls*(36.*mls2 + 2.*pv12 + (pv2s*pv2s) - 6.*mls*(pv1s + pv2s)) 
		+ (12.*mls2 + 3.*pv12 - 2.*mls*(3.*pv1s + 4.*pv2s))*ps2)*C2 
	  - 4.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
		+ (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C11 
	  - 8.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
		+ (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C12 
	  - 4.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
		+ (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C22 )/mls3;

      c+= -lc* 
	(-12.*mls2 + 8.*mls*pv1s - (pv1s*pv1s) + 48.*B00B*(2.*mls - pv1s) 
	 + 24.*B00C*(2.*mls - pv2s) - 6.*A0A*(6.*mls + pv1s - 2.*ps2) 
	 - 6.*B1B*(2.*mls - pv1s)*(2.*mls - pv1s + 3.*pv2s - ps2) 
	 - 24.*mls*B11A*(2.*mls + ps2) 
	 - 3.*B1A*(4.*mls - pv1s + pv2s - ps2)*(2.*mls + ps2) 
	 - 3.*B1A*(2.*mls + ps2)*(4.*mls + pv1s - pv2s + ps2) 
	 + 6.*B1C*(2.*mls - pv2s)*(-3.*pv1s + pv2s + ps2) 
	 + 12.*B11C*(2.*mls - pv2s)*(-pv1s + pv2s + ps2) 
	 + 6.*B11B*(2.*mls - pv1s)*(3.*pv1s - 2.*pv2s + 2.*ps2) 
	 + 6.*B0C*(2.*mls*pv1s + (4.*mls + pv1s)*pv2s - 4.*mls*ps2) 
	 - 6.*B0B*(2.*mls2 - 5.*mls*pv1s - (2.*mls + pv1s)*pv2s 
		  + 4.*mls*ps2) 
	 - 3.*B0A*(2.*mls*(4.*mls + pv1s + 2.*pv2s) - (4.*mls + pv1s)*ps2) 
	 - 3.*B0A*(2.*mls*(4.*mls + 2.*pv1s + pv2s) - (2.*mls + pv2s)*ps2 
		  + (ps2*ps2)) 
	 - 6.*C0A*(2.*mls*(12.*mls2 - (pv1s*pv1s) + pv12 + (pv2s*pv2s) 
			 - 6.*mls*(pv1s + pv2s)) 
		  + (12.*mls2 + pv12 + 2.*mls*(pv1s - pv2s))*ps2 - 2.*mls*(ps2*ps2)) 
	 + 24.*mls*(2.*mls*B0A + (B1A + B11A)*(2.*mls + ps2)) 
	 - 24.*(mls*(24.*mls2 - (pv1s*pv1s) + 2.*pv12 + (pv2s*pv2s) 
		    - 4.*mls*(pv1s + pv2s)) 
	       + (4.*mls2 + pv12 - mls*(pv1s + 3.*pv2s))*ps2)*C1 
	 - 12.*(24.*mls3 - 2.*mls*(pv1s*pv1s) - 4.*mls2*(pv1s + pv2s) 
	       + (4.*mls2 - 2.*mls*pv2s + pv12)*ps2)*C2 
	 - 24.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
	       +(2.*mls - pv1s)*(2.*mls - pv2s)*ps2)* C11 
	 -  24.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
		+ (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C12 
	 + 72.*(2.*mls - pv1s)*(-0.125*mls - 0.125*A0A + (0.25*mls*B1B) 
			      + (0.125*pv1s*B11B)) 
	 - 12.*(-2.*mls + pv1s)*(0.25*mls - 0.25*A0A - (0.5*mls*B1B) 
			       - (0.75*pv1s*B11B)) )/24./mls3;
      
      d+= -lc*
	(-2.*mls2*B0A - 2.*mls*C0A*(8.*mls2 + pv12 - 2.*mls*(pv1s + pv2s)) 	
	 - mls*B1A*(2.*mls + ps2) - mls*B11A*(2.*mls + ps2) 
	 + mls*(2.*mls*B0A + (B1A + B11A)*(2.*mls + ps2)) 
	 - (24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
	    + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)* C1 -  2.*mls*pv12*C2 
	 - (24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
	    + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C11 
	 - (24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s)
	    +(2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C12 )/mls3;
      
      e+= -0.25*lc*
	(8.*mls*B0B*pv1s + 4.*B00B*(2*mls - pv1s) - 2.*A0A*(2.*mls + pv1s) 
	 - 4.*B11C*pv1s*(2.*mls - pv2s) + 2.*B1C*(4.*mls - 3.*pv1s)*(2.*mls - pv2s) 
	 + 2.*B0C*(4.*mls2 - 2.*mls*pv2s + pv12) - 2.*mls*B0A*(pv1s - ps2) 
	 + 4.*mls*C0A*pv1s*(4.*mls - pv2s - ps2) - 4.*mls*B11A*(2.*mls + ps2) 
	 + B1A*(2.*mls - pv1s)*(2.*mls + ps2) - B1A*(6.*mls - pv1s)*(2.*mls + ps2) 
	 - B0A*(8.*mls2 + (2.*mls - pv1s)*ps2) 
	 + 4.*mls*(2.*mls*B0A + (B1A + B11A)*(2.*mls + ps2)) 
	 - 2*(2.*mls*(12.*mls2 - pv1s + 2.*pv12 - 2.*mls*(pv1s + pv2s)) 
	      + (4.*mls2 - 2.*mls*pv2s + pv12)*ps2)*C1 
	 - 4.*(24.*mls3 + 2.*mls*pv12 - 4.*mls2*(pv1s + pv2s) 
	      + (2.*mls - pv1s)*(2.*mls - pv2s)*ps2)*C11 )/mls3;


      cerr << a << ' ' 
	   << b << ' '
	   << c << ' '
	   << d << ' '
	   << e << ' '
	   << f << '\n';
    }     LT::clearcache();
    {
      int theC = LT::Cget(ps2,pv2s,pv1s,
		      mls,mls,mls);
      Complex C1 = LT::Cval(LT::cc1,theC);
      Complex C2 = LT::Cval(LT::cc2,theC);
      Complex C00 = LT::Cval(LT::cc00,theC);
      Complex C11 = LT::Cval(LT::cc11,theC);
      Complex C12 = LT::Cval(LT::cc12,theC);
      Complex C22 = LT::Cval(LT::cc22,theC);
      Complex Cz = LT::C0(pv1s,pv2s,ps2,mls,mls,mls);
      /**
       * vector type can contain different types of particle 
       * and hence the coupling is different
       * Here left[i] is used for the coupling of the ith 
       * type rather than creating another
       * vector to hold them.
       */
      a +=  4.*lc*(LT::B0(ps2,mls,mls) - 4.*C00)/ps2;
      b += -4.*lc*(Cz + 3.*C1 + 3.*C2 +2.*(C11 + 2.*C12 + C22 ));
      c += -2.*lc*(Cz + 2.*(2.*C1+ C2 + 2.*(C11 +C12)));
      d += -8.*lc*(C1 +C11 + C12);
      e += -4.*lc*(C1 + 2.*C11);

      cerr << a << ' ' 
	   << b << ' '
	   << c << ' '
	   << d << ' '
	   << e << ' '
	   << f << '\n';
    }

  }

 LT::clearcache();

  LT::ltexi();










  return 0;
}
