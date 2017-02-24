/* -*- C++ -*-
   clooptools.h
   the C/C++ header file with all definitions for LoopTools
   this file is part of LoopTools
   last modified 2 Sep 14 th
   dgrell 2011-01-21 for Herwig: moved definitions and extern declarations to our clooptools.cc
   dgrell 2016-09-08 for Herwig: update to Looptols-2.13
*/


#ifndef HERWIG_clooptools_h_
#define HERWIG_clooptools_h_

#include <complex>

/****************************************************************/

// Adapted version of ftypes.h


typedef double RealType;
typedef double REAL;

typedef int INTEGER;
typedef const REAL CREAL;

struct COMPLEX { REAL re, im; };
typedef const COMPLEX CCOMPLEX;

typedef std::complex<RealType> ComplexType;

typedef const RealType cRealType;
typedef const ComplexType cComplexType;

//#define cachelookup_ ljcachelookup_

// typedef const COMPLEX CCOMPLEX;

typedef long long int memindex;

/****************************************************************/

#define AARGS(t) t(m)

#define BARGS(t) t(p), t(m1), t(m2)

#define CARGS(t) t(p1), t(p2), t(p1p2), t(m1), t(m2), t(m3)

#define DARGS(t) t(p1), t(p2), t(p3), t(p4), t(p1p2), t(p2p3),	\
    t(m1), t(m2), t(m3), t(m4)

#define EARGS(t) t(p1), t(p2), t(p3), t(p4), t(p5),	\
    t(p1p2), t(p2p3), t(p3p4), t(p4p5), t(p5p1),	\
    t(m1), t(m2), t(m3), t(m4), t(m5)

/****************************************************************/

#define _lt_Cr_(v) cRealType v
#define _lt_Cc_(v) cComplexType v
#define _lt_Fr_(v) CREAL *v
#define _lt_Fc_(v) CCOMPLEX *v
#define _lt_Id_(v) v
#define _lt_Frp_(v) &v
#define _lt_Fcp_(v) (CCOMPLEX *)&v
#define _lt_Fap_(v) (COMPLEX *)v

// #define _lt_Frd_(f)
// #define _lt_Fcd_(f)
// #define _lt_Fad_(v,n)
// #define _lt_Fax_(v,n)

/****************************************************************/

namespace Herwig {
  namespace Looptools {

    enum {
      aa0 = 0, 
      aa00 = 3, 
      Naa = 6
    };

    enum {
      bb0 = 0,
      bb1 = 3,
      bb00 = 6,
      bb11 = 9,
      bb001 = 12,
      bb111 = 15,
      dbb0 = 18,
      dbb1 = 21,
      dbb00 = 24,
      dbb11 = 27,
      dbb001 = 30,
      Nbb = 33
    };

    enum {
      cc0 = 0,
      cc1 = 3,
      cc2 = 6,
      cc00 = 9,
      cc11 = 12,
      cc12 = 15,
      cc22 = 18,
      cc001 = 21,
      cc002 = 24,
      cc111 = 27,
      cc112 = 30,
      cc122 = 33,
      cc222 = 36,
      cc0000 = 39,
      cc0011 = 42,
      cc0012 = 45,
      cc0022 = 48,
      cc1111 = 51,
      cc1112 = 54,
      cc1122 = 57,
      cc1222 = 60,
      cc2222 = 63,
      Ncc = 66
    };

    enum {
      dd0 = 0,
      dd1 = 3,
      dd2 = 6,
      dd3 = 9,
      dd00 = 12,
      dd11 = 15,
      dd12 = 18,
      dd13 = 21,
      dd22 = 24,
      dd23 = 27,
      dd33 = 30,
      dd001 = 33,
      dd002 = 36,
      dd003 = 39,
      dd111 = 42,
      dd112 = 45,
      dd113 = 48,
      dd122 = 51,
      dd123 = 54,
      dd133 = 57,
      dd222 = 60,
      dd223 = 63,
      dd233 = 66,
      dd333 = 69,
      dd0000 = 72,
      dd0011 = 75,
      dd0012 = 78,
      dd0013 = 81,
      dd0022 = 84,
      dd0023 = 87,
      dd0033 = 90,
      dd1111 = 93,
      dd1112 = 96,
      dd1113 = 99,
      dd1122 = 102,
      dd1123 = 105,
      dd1133 = 108,
      dd1222 = 111,
      dd1223 = 114,
      dd1233 = 117,
      dd1333 = 120,
      dd2222 = 123,
      dd2223 = 126,
      dd2233 = 129,
      dd2333 = 132,
      dd3333 = 135,
      dd00001 = 138,
      dd00002 = 141,
      dd00003 = 144,
      dd00111 = 147,
      dd00112 = 150,
      dd00113 = 153,
      dd00122 = 156,
      dd00123 = 159,
      dd00133 = 162,
      dd00222 = 165,
      dd00223 = 168,
      dd00233 = 171,
      dd00333 = 174,
      dd11111 = 177,
      dd11112 = 180,
      dd11113 = 183,
      dd11122 = 186,
      dd11123 = 189,
      dd11133 = 192,
      dd11222 = 195,
      dd11223 = 198,
      dd11233 = 201,
      dd11333 = 204,
      dd12222 = 207,
      dd12223 = 210,
      dd12233 = 213,
      dd12333 = 216,
      dd13333 = 219,
      dd22222 = 222,
      dd22223 = 225,
      dd22233 = 228,
      dd22333 = 231,
      dd23333 = 234,
      dd33333 = 237,
      Ndd = 240
    };

    enum {
      ee0 = 0,
      ee1 = 3,
      ee2 = 6,
      ee3 = 9,
      ee4 = 12,
      ee00 = 15,
      ee11 = 18,
      ee12 = 21,
      ee13 = 24,
      ee14 = 27,
      ee22 = 30,
      ee23 = 33,
      ee24 = 36,
      ee33 = 39,
      ee34 = 42,
      ee44 = 45,
      ee001 = 48,
      ee002 = 51,
      ee003 = 54,
      ee004 = 57,
      ee111 = 60,
      ee112 = 63,
      ee113 = 66,
      ee114 = 69,
      ee122 = 72,
      ee123 = 75,
      ee124 = 78,
      ee133 = 81,
      ee134 = 84,
      ee144 = 87,
      ee222 = 90,
      ee223 = 93,
      ee224 = 96,
      ee233 = 99,
      ee234 = 102,
      ee244 = 105,
      ee333 = 108,
      ee334 = 111,
      ee344 = 114,
      ee444 = 117,
      ee0000 = 120,
      ee0011 = 123,
      ee0012 = 126,
      ee0013 = 129,
      ee0014 = 132,
      ee0022 = 135,
      ee0023 = 138,
      ee0024 = 141,
      ee0033 = 144,
      ee0034 = 147,
      ee0044 = 150,
      ee1111 = 153,
      ee1112 = 156,
      ee1113 = 159,
      ee1114 = 162,
      ee1122 = 165,
      ee1123 = 168,
      ee1124 = 171,
      ee1133 = 174,
      ee1134 = 177,
      ee1144 = 180,
      ee1222 = 183,
      ee1223 = 186,
      ee1224 = 189,
      ee1233 = 192,
      ee1234 = 195,
      ee1244 = 198,
      ee1333 = 201,
      ee1334 = 204,
      ee1344 = 207,
      ee1444 = 210,
      ee2222 = 213,
      ee2223 = 216,
      ee2224 = 219,
      ee2233 = 222,
      ee2234 = 225,
      ee2244 = 228,
      ee2333 = 231,
      ee2334 = 234,
      ee2344 = 237,
      ee2444 = 240,
      ee3333 = 243,
      ee3334 = 246,
      ee3344 = 249,
      ee3444 = 252,
      ee4444 = 255,
      Nee = 258
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
      DebugA = 1,
      DebugB = 1<<1,
      DebugC = 1<<2,
      DebugD = 1<<3,
      DebugE = 1<<4,
      DebugAll = DebugA + DebugB + DebugC + DebugD + DebugE
    };

    ComplexType ToComplex(COMPLEX c);

    /**
     *  Looptools initialisation
     */
    void ltini(std::string logfilename = std::string("Looptools.log"));

    /**
     *  Looptools termination
     */
    void ltexi(std::string logfilename = std::string("Looptools.log"));

    /****************************************************************/

    memindex Aget(AARGS(_lt_Cr_));
    memindex AgetC(AARGS(_lt_Cc_));

    void Aput(ComplexType *res, AARGS(_lt_Cr_));
    void AputC(ComplexType *res, AARGS(_lt_Cc_));

    void Aputnocache(ComplexType *res, AARGS(_lt_Cr_));
    void AputnocacheC(ComplexType *res, AARGS(_lt_Cc_));

    COMPLEX *Acache(const memindex integral);
    COMPLEX *AcacheC(const memindex integral);

    ComplexType Aval(const int i, const memindex integral);
    ComplexType AvalC(const int i, const memindex integral);

    ComplexType A0i(const int i, AARGS(_lt_Cr_));
    ComplexType A0iC(const int i, AARGS(_lt_Cc_));

    ComplexType A0(AARGS(_lt_Cr_));
    ComplexType A00(AARGS(_lt_Cr_));

    ComplexType A0C(AARGS(_lt_Cc_));
    ComplexType A00C(AARGS(_lt_Cc_));

    /****************************************************************/

    memindex Bget(BARGS(_lt_Cr_));
    memindex BgetC(BARGS(_lt_Cc_));

    void Bput(ComplexType *res, BARGS(_lt_Cr_));
    void BputC(ComplexType *res, BARGS(_lt_Cc_));

    void Bputnocache(ComplexType *res, BARGS(_lt_Cr_));
    void BputnocacheC(ComplexType *res, BARGS(_lt_Cc_));

    COMPLEX *Bcache(const memindex integral);
    COMPLEX *BcacheC(const memindex integral);

    ComplexType Bval(const int i, const memindex integral);
    ComplexType BvalC(const int i, const memindex integral);

    ComplexType B0i(const int i, BARGS(_lt_Cr_));
    ComplexType B0iC(const int i, BARGS(_lt_Cc_));

    ComplexType B0(BARGS(_lt_Cr_));
    ComplexType B1(BARGS(_lt_Cr_));
    ComplexType B00(BARGS(_lt_Cr_));
    ComplexType B11(BARGS(_lt_Cr_));
    ComplexType B001(BARGS(_lt_Cr_));
    ComplexType B111(BARGS(_lt_Cr_));
    ComplexType DB0(BARGS(_lt_Cr_));
    ComplexType DB1(BARGS(_lt_Cr_));
    ComplexType DB00(BARGS(_lt_Cr_));
    ComplexType DB11(BARGS(_lt_Cr_));

    ComplexType B0C(BARGS(_lt_Cc_));
    ComplexType B1C(BARGS(_lt_Cc_));
    ComplexType B00C(BARGS(_lt_Cc_));
    ComplexType B11C(BARGS(_lt_Cc_));
    ComplexType B001C(BARGS(_lt_Cc_));
    ComplexType B111C(BARGS(_lt_Cc_));
    ComplexType DB0C(BARGS(_lt_Cc_));
    ComplexType DB1C(BARGS(_lt_Cc_));
    ComplexType DB00C(BARGS(_lt_Cc_));
    ComplexType DB11C(BARGS(_lt_Cc_));

    /****************************************************************/

    memindex Cget(CARGS(_lt_Cr_));
    memindex CgetC(CARGS(_lt_Cc_));

    void Cput(ComplexType *res, CARGS(_lt_Cr_));
    void CputC(ComplexType *res, CARGS(_lt_Cc_));

    void C0nocache(ComplexType *res, CARGS(_lt_Cr_));
    void C0nocacheC(ComplexType *res, CARGS(_lt_Cc_));

    COMPLEX *Ccache(const memindex integral);
    COMPLEX *CcacheC(const memindex integral);

    ComplexType Cval(const int i, const memindex integral);
    ComplexType CvalC(const int i, const memindex integral);

    ComplexType C0i(const int i, CARGS(_lt_Cr_));
    ComplexType C0iC(const int i, CARGS(_lt_Cc_));

    ComplexType C0(CARGS(_lt_Cr_));
    ComplexType C0C(CARGS(_lt_Cc_));

    /****************************************************************/

    memindex Dget(DARGS(_lt_Cr_));
    memindex DgetC(DARGS(_lt_Cc_));

    void Dput(ComplexType *res, DARGS(_lt_Cr_));
    void DputC(ComplexType *res, DARGS(_lt_Cc_));

    void D0nocache(ComplexType *res, DARGS(_lt_Cr_));
    void D0nocacheC(ComplexType *res, DARGS(_lt_Cc_));

    COMPLEX *Dcache(const memindex integral);
    COMPLEX *DcacheC(const memindex integral);

    ComplexType Dval(const int i, const memindex integral);
    ComplexType DvalC(const int i, const memindex integral);

    ComplexType D0i(const int i, DARGS(_lt_Cr_));
    ComplexType D0iC(const int i, DARGS(_lt_Cc_));

    ComplexType D0(DARGS(_lt_Cr_));
    ComplexType D0C(DARGS(_lt_Cc_));

    /****************************************************************/

    memindex Eget(EARGS(_lt_Cr_));
    memindex EgetC(EARGS(_lt_Cc_));

    void Eput(ComplexType *res, EARGS(_lt_Cr_));
    void EputC(ComplexType *res, EARGS(_lt_Cc_));

    void E0nocache(ComplexType *res, EARGS(_lt_Cr_));
    void E0nocacheC(ComplexType *res, EARGS(_lt_Cc_));

    COMPLEX *Ecache(const memindex integral);
    COMPLEX *EcacheC(const memindex integral);

    ComplexType Eval(const int i, const memindex integral);
    ComplexType EvalC(const int i, const memindex integral);

    ComplexType E0i(const int i, EARGS(_lt_Cr_));
    ComplexType E0iC(const int i, EARGS(_lt_Cc_));

    ComplexType E0(EARGS(_lt_Cr_));
    ComplexType E0C(EARGS(_lt_Cc_));

    /****************************************************************/

    ComplexType Li2(_lt_Cr_(x));
    ComplexType Li2C(_lt_Cc_(x));
    ComplexType Li2omx(_lt_Cr_(x));
    ComplexType Li2omxC(_lt_Cc_(x));

    /****************************************************************/

    void clearcache();
    void markcache();
    void restorecache();

    void setmudim(cRealType mudim);
    RealType getmudim();

    void setdelta(cRealType delta);
    RealType getdelta();

    void setuvdiv(cRealType uvdiv);
    RealType getuvdiv();

    void setlambda(cRealType lambda);
    RealType getlambda();

    void setminmass(cRealType minmass);
    RealType getminmass();

    void setmaxdev(cRealType maxdev);
    RealType getmaxdev();

    void setwarndigits(const int warndigits);
    int getwarndigits();

    void seterrdigits(const int errdigits);
    int geterrdigits();

    void setversionkey(const int versionkey);
    int getversionkey();

    void setdebugkey(const int debugkey);
    int getdebugkey();

    void setdebugrange(const int debugfrom, const int debugto);

    void setcmpbits(const int cmpbits);
    int getcmpbits();

    void setdiffeps(cRealType diffeps);
    RealType getdiffeps();

    void setzeroeps(cRealType zeroeps);
    RealType getzeroeps();

  } // namespace Looptools
} // namespace Herwig

#endif

