/* -*- C++ -*-
   clooptools.h
   the C/C++ header file with all definitions for LoopTools
   this file is part of LoopTools
   last modified 9 Dec 10 th
   dgrell 2011-01-21 for Herwig: moved definitions and extern declarations to our clooptools.cc
*/


#ifndef HERWIG_clooptools_h_
#define HERWIG_clooptools_h_

//#define cachelookup_ ljcachelookup_

/** complex defn for Looptools
 */
struct DOUBLE_COMPLEX {
  /** 
   * Real part
   */
  double re;

  /**
   * Imaginary part
   */ 
  double im; 
};
typedef const DOUBLE_COMPLEX CDOUBLE_COMPLEX;

#include <complex>
typedef std::complex<double> double_complex;

/****************************************************************/

#define AARGS(t) t(m)

#define BARGS(t) t(p), t(m1), t(m2)

#define CARGS(t) t(p1), t(p2), t(p1p2), t(m1), t(m2), t(m3)

#define DARGS(t) t(p1), t(p2), t(p3), t(p4), t(p1p2), t(p2p3),	\
    t(m1), t(m2), t(m3), t(m4)

#define EARGS(t) t(p1), t(p2), t(p3), t(p4), t(p5),	\
    t(p1p2), t(p2p3), t(p3p4), t(p4p5), t(p5p1),	\
    t(m1), t(m2), t(m3), t(m4), t(m5)

#define _Cr_(v) const double v
#define _Cc_(v) const double_complex v
#define _Fr_(v) const double *v
#define _Fc_(v) CDOUBLE_COMPLEX *v
#define _Frp_(v) &v
#define _Fcp_(v) (CDOUBLE_COMPLEX *)&v
#define _Id_(v) v

/****************************************************************/

namespace Herwig {
  namespace Looptools {

    enum {
      bb0, bb1, bb00, bb11, bb001, bb111, dbb0, dbb1, dbb00, dbb11,
      Nbb
    };

    enum {
      cc0, cc1, cc2, cc00, cc11, cc12, cc22, cc001, cc002, cc111, cc112,
      cc122, cc222, cc0000, cc0011, cc0012, cc0022, cc1111, cc1112, cc1122,
      cc1222, cc2222,
      Ncc
    };

    enum {
      dd0, dd1, dd2, dd3, dd00, dd11, dd12, dd13, dd22, dd23, dd33,
      dd001, dd002, dd003, dd111, dd112, dd113, dd122, dd123, dd133, dd222,
      dd223, dd233, dd333, dd0000, dd0011, dd0012, dd0013, dd0022, dd0023,
      dd0033, dd1111, dd1112, dd1113, dd1122, dd1123, dd1133, dd1222, 
      dd1223, dd1233, dd1333, dd2222, dd2223, dd2233, dd2333, dd3333, 
      dd00001, dd00002, dd00003, dd00111, dd00112, dd00113, dd00122, 
      dd00123, dd00133, dd00222, dd00223, dd00233, dd00333, dd11111, 
      dd11112, dd11113, dd11122, dd11123, dd11133, dd11222, dd11223, 
      dd11233, dd11333, dd12222, dd12223, dd12233, dd12333, dd13333, 
      dd22222, dd22223, dd22233, dd22333, dd23333, dd33333,
      Ndd
    };

    enum {
      ee0, ee1, ee2, ee3, ee4, ee00, ee11, ee12, ee13, ee14, ee22, ee23, 
      ee24, ee33, ee34, ee44, ee001, ee002, ee003, ee004, ee111, ee112, 
      ee113, ee114, ee122, ee123, ee124, ee133, ee134, ee144, ee222,
      ee223, ee224, ee233, ee234, ee244, ee333, ee334, ee344, ee444,
      ee0000, ee0011, ee0012, ee0013, ee0014, ee0022, ee0023, ee0024,
      ee0033, ee0034, ee0044, ee1111, ee1112, ee1113, ee1114, ee1122, 
      ee1123, ee1124, ee1133, ee1134, ee1144, ee1222, ee1223, ee1224,
      ee1233, ee1234, ee1244, ee1333, ee1334, ee1344, ee1444, ee2222,
      ee2223, ee2224, ee2233, ee2234, ee2244, ee2333, ee2334, ee2344,
      ee2444, ee3333, ee3334, ee3344, ee3444, ee4444,
      Nee
    };

    enum {
      KeyA0 = 1,
      KeyBget = 1<<2,
      KeyC0 = 1<<4,
      KeyD0 = 1<<6,
      KeyE0 = 1<<8,
      KeyEget = 1<<10,
      KeyEgetC = 1<<12,
      KeyAll = KeyA0 + KeyBget + KeyC0 + KeyD0 + KeyE0 + KeyEget + KeyEgetC
    };

    enum {
      DebugB = 1,
      DebugC = 1<<1,
      DebugD = 1<<2,
      DebugE = 1<<3,
      DebugAll = DebugB + DebugC + DebugD + DebugE
    };

    double_complex ToComplex(DOUBLE_COMPLEX c);

    /**
     *  Looptools initialisation
     */
    void ltini(std::string logfilename = std::string("Looptools.log"));

    /**
     *  Looptools termination
     */
    void ltexi(std::string logfilename = std::string("Looptools.log"));


    double_complex A0(AARGS(_Cr_));
    double_complex A0C(AARGS(_Cc_));
    double_complex A00(AARGS(_Cr_));

    double_complex A00C(AARGS(_Cc_));

    /****************************************************************/

    long Bget(BARGS(_Cr_));
    long BgetC(BARGS(_Cc_));
    DOUBLE_COMPLEX *Bcache(const long integral);
    DOUBLE_COMPLEX *BcacheC(const long integral);

    double_complex Bval(const int i, const long integral);

    double_complex BvalC(const int i, const long integral);

    double_complex B0i(const int i, BARGS(_Cr_));

    double_complex B0iC(const int i, BARGS(_Cc_));

    double_complex B0(BARGS(_Cr_));
    double_complex B1(BARGS(_Cr_));
    double_complex B00(BARGS(_Cr_));
    double_complex B11(BARGS(_Cr_));
    double_complex B001(BARGS(_Cr_));
    double_complex B111(BARGS(_Cr_));
    double_complex DB0(BARGS(_Cr_));
    double_complex DB1(BARGS(_Cr_));
    double_complex DB00(BARGS(_Cr_));
    double_complex DB11(BARGS(_Cr_));

    double_complex B0C(BARGS(_Cc_));
    double_complex B1C(BARGS(_Cc_));
    double_complex B00C(BARGS(_Cc_));
    double_complex B11C(BARGS(_Cc_));
    double_complex B001C(BARGS(_Cc_));
    double_complex B111C(BARGS(_Cc_));
    double_complex DB0C(BARGS(_Cc_));
    double_complex DB1C(BARGS(_Cc_));
    double_complex DB00C(BARGS(_Cc_));
    double_complex DB11C(BARGS(_Cc_));

    /****************************************************************/

    double_complex C0(CARGS(_Cr_));

    double_complex C0C(CARGS(_Cc_));

    long Cget(CARGS(_Cr_));

    long CgetC(CARGS(_Cc_));

    DOUBLE_COMPLEX *Ccache(const long integral);

    DOUBLE_COMPLEX *CcacheC(const long integral);

    double_complex Cval(const int i, const long integral);

    double_complex CvalC(const int i, const long integral);

    double_complex C0i(const int i, CARGS(_Cr_));

    double_complex C0iC(const int i, CARGS(_Cc_));

    /****************************************************************/

    double_complex D0(DARGS(_Cr_));

    double_complex D0C(DARGS(_Cc_));

    long Dget(DARGS(_Cr_));

    long DgetC(DARGS(_Cc_));

    DOUBLE_COMPLEX *Dcache(const long integral);

    DOUBLE_COMPLEX *DcacheC(const long integral);

    double_complex Dval(const int i, const long integral);

    double_complex DvalC(const int i, const long integral);

    double_complex D0i(const int i, DARGS(_Cr_));

    double_complex D0iC(const int i, DARGS(_Cc_));

    /****************************************************************/

    double_complex E0(EARGS(_Cr_));

    double_complex E0C(EARGS(_Cc_));

    long Eget(EARGS(_Cr_));

    long EgetC(EARGS(_Cc_));

    DOUBLE_COMPLEX *Ecache(const long integral);

    DOUBLE_COMPLEX *EcacheC(const long integral);

    double_complex Eval(const int i, const long integral);

    double_complex EvalC(const int i, const long integral);

    double_complex E0i(const int i, EARGS(_Cr_));

    double_complex E0iC(const int i, EARGS(_Cc_));

    /****************************************************************/

    double_complex Li2(const double x);

    double_complex Li2C(const double_complex x);

    /****************************************************************/


    void setmudim(const double mudim);

    double getmudim();


    void setdelta(const double delta);

    double getdelta();


    void setlambda(const double lambda);
    double getlambda();


    void setminmass(const double minmass);
    double getminmass();


    void setmaxdev(const double maxdev);
    double getmaxdev();


    void setwarndigits(const long warndigits);
    long getwarndigits();


    void seterrdigits(const long errdigits);
    long geterrdigits();


    void setversionkey(const long versionkey);
    long getversionkey();

    void setdebugkey(const long debugkey);
    long getdebugkey();


    void setdebugrange(const long debugfrom, const long debugto);

    void setcmpbits(const long cmpbits);

    long getcmpbits();

    void clearcache();
    void markcache();
    void restorecache();

  } // namespace Looptools
} // namespace Herwig

#endif

