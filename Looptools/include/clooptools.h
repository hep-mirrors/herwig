/* -*- C++ -*-
	clooptools.h
		the C/C++ header file with all definitions for LoopTools
		this file is part of LoopTools
		last modified 21 Dec 06 th
		dgrell 2008-01 for Herwig++
*/


#ifndef HERWIG_clooptools_h_
#define HERWIG_clooptools_h_

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

extern "C" {

extern void a0sub_(DOUBLE_COMPLEX *result, const double *m);
extern void a00sub_(DOUBLE_COMPLEX *result, const double *m);

extern long bget_(const double *p, const double *m1, const double *m2);

extern void c0sub_(DOUBLE_COMPLEX *result, 
		   const double *p1, const double *p2, 
		   const double *p1p2, 
		   const double *m1, const double *m2, const double *m3);
extern long cget_(const double *p1, const double *p2, 
		  const double *p1p2, 
		  const double *m1, const double *m2, const double *m3);

extern void d0sub_(DOUBLE_COMPLEX *result, 
		   const double *p1, const double *p2, 
		   const double *p3, const double *p4, 
		   const double *p1p2, const double *p2p3, 
		   const double *m1, const double *m2, 
		   const double *m3, const double *m4);
extern long dget_(const double *p1, const double *p2, 
		  const double *p3, const double *p4, 
		  const double *p1p2, const double *p2p3, 
		  const double *m1, const double *m2, 
		  const double *m3, const double *m4);

extern void e0sub_(DOUBLE_COMPLEX *result, 
		   const double *p1, const double *p2, 
		   const double *p3, const double *p4, const double *p5, 
		   const double *p1p2, const double *p2p3, 
		   const double *p3p4, const double *p4p5, 
		   const double *p5p1, 
		   const double *m1, const double *m2, 
		   const double *m3, const double *m4, const double *m5);
extern long eget_(const double *p1, const double *p2, 
		  const double *p3, const double *p4, const double *p5, 
		  const double *p1p2, const double *p2p3, 
		  const double *p3p4, const double *p4p5, 
		  const double *p5p1, 
		  const double *m1, const double *m2, 
		  const double *m3, const double *m4, const double *m5);


extern void li2sub_(DOUBLE_COMPLEX *result, const double *x);
extern void li2csub_(DOUBLE_COMPLEX *result, CDOUBLE_COMPLEX *x);

extern void ffini_(void);
extern void ffexi_(void);

extern void clearcache_(void);
extern void markcache_(void);
extern void restorecache_(void);

extern struct {		/* MUST match common block ltvars in lt.h! */
  DOUBLE_COMPLEX cache[8][2];
  DOUBLE_COMPLEX savedptr[8];
  double maxdev;
  long serial, warndigits, errdigits, versionkey;
  long debugkey, debugfrom, debugto;
} ltvars_;

extern struct {		/* MUST match common block ffregul in ff.h! */
  double mudim, delta, lambda;
} ffregul_;

}

/****************************************************************/

namespace Herwig {
  namespace Looptools {

    inline double_complex ToComplex(DOUBLE_COMPLEX c) {
      return double_complex(c.re, c.im);
    }

    /**
     *  Looptools initialisation
     */
    void ffini(std::string logfilename = std::string("Looptools.log"));

    /**
     *  Looptools termination
     */
    void ffexi(std::string logfilename = std::string("Looptools.log"));

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
  KeyALL = KeyA0 + KeyBget + KeyC0 + KeyD0 + KeyE0 + KeyEget + KeyEgetC
};

enum {
  DebugB = 1,
  DebugC = 1<<1,
  DebugD = 1<<2,
  DebugE = 1<<3,
  DebugAll = DebugB + DebugC + DebugD + DebugE
};


inline double_complex A0(const double m)
{
  DOUBLE_COMPLEX result;
  a0sub_(&result, &m);
  return ToComplex(result);
}

inline double_complex A00(const double m)
{
  DOUBLE_COMPLEX result;
  a00sub_(&result, &m);
  return ToComplex(result);
}

/****************************************************************/

inline long Bget(const double p, const double m1, const double m2)
{
  return bget_(&p, &m1, &m2);
}

inline DOUBLE_COMPLEX *Bcache(const long integral)
  { return &ltvars_.cache[0][integral]; }

inline DOUBLE_COMPLEX *BcacheC(const long integral)
  { return &ltvars_.cache[1][integral]; }

inline double_complex Bval(const long i, const long integral)
  { return ToComplex(Bcache(integral)[i]); }

inline double_complex BvalC(const long i, const long integral)
  { return ToComplex(BcacheC(integral)[i]); }

inline double_complex B0i(const long i, const double p, 
			  const double m1, const double m2)
  { return Bval(i, Bget(p, m1, m2)); }

inline double_complex B0(const double p, const double m1, const double m2)
  { return B0i(bb0, p, m1, m2); }
inline double_complex B1(const double p, const double m1, const double m2)
  { return B0i(bb1, p, m1, m2); }
inline double_complex B00(const double p, const double m1, const double m2)
  { return B0i(bb00, p, m1, m2); }
inline double_complex B11(const double p, const double m1, const double m2)
  { return B0i(bb11, p, m1, m2); }
inline double_complex B001(const double p, const double m1, const double m2)
  { return B0i(bb001, p, m1, m2); }
inline double_complex B111(const double p, const double m1, const double m2)
  { return B0i(bb111, p, m1, m2); }
inline double_complex DB0(const double p, const double m1, const double m2)
  { return B0i(dbb0, p, m1, m2); }
inline double_complex DB1(const double p, const double m1, const double m2)
  { return B0i(dbb1, p, m1, m2); }
inline double_complex DB00(const double p, const double m1, const double m2)
  { return B0i(dbb00, p, m1, m2); }
inline double_complex DB11(const double p, const double m1, const double m2)
  { return B0i(dbb11, p, m1, m2); }

/****************************************************************/

inline double_complex C0(const double p1, const double p2, 
			 const double p1p2, 
			 const double m1, const double m2, const double m3)
{
  DOUBLE_COMPLEX result;
  c0sub_(&result, &p1, &p2, &p1p2, &m1, &m2, &m3);
  return ToComplex(result);
}

inline long Cget(const double p1, const double p2, 
		 const double p1p2, const double m1, 
		 const double m2, const double m3)
{
  return cget_(&p1, &p2, &p1p2, &m1, &m2, &m3);
}

inline DOUBLE_COMPLEX *Ccache(const long integral)
  { return &ltvars_.cache[2][integral]; }

inline DOUBLE_COMPLEX *CcacheC(const long integral)
  { return &ltvars_.cache[3][integral]; }

inline double_complex Cval(const long i, const long integral)
  { return ToComplex(Ccache(integral)[i]); }

inline double_complex CvalC(const long i, const long integral)
  { return ToComplex(CcacheC(integral)[i]); }

inline double_complex C0i(const long i, 
			  const double p1, const double p2, 
			  const double p1p2, const double m1, 
			  const double m2, const double m3)
  { return Cval(i, Cget(p1, p2, p1p2, m1, m2, m3)); }

/****************************************************************/

inline double_complex D0(const double p1, const double p2, 
			 const double p3, const double p4, 
			 const double p1p2, const double p2p3, 
			 const double m1, const double m2, 
			 const double m3, const double m4)
{
  DOUBLE_COMPLEX result;
  d0sub_(&result, &p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
  return ToComplex(result);
}

inline long Dget(const double p1, const double p2, 
		 const double p3, const double p4, 
		 const double p1p2, const double p2p3, 
		 const double m1, const double m2, 
		 const double m3, const double m4)
{
  return dget_(&p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
}

inline DOUBLE_COMPLEX *Dcache(const long integral)
  { return &ltvars_.cache[4][integral]; }

inline DOUBLE_COMPLEX *DcacheC(const long integral)
  { return &ltvars_.cache[5][integral]; }

inline double_complex Dval(const long i, const long integral)
  { return ToComplex(Dcache(integral)[i]); }

inline double_complex DvalC(const long i, const long integral)
  { return ToComplex(DcacheC(integral)[i]); }

inline double_complex D0i(const long i, 
			  const double p1, const double p2, 
			  const double p3, const double p4, 
			  const double p1p2, const double p2p3, 
			  const double m1, const double m2, 
			  const double m3, const double m4)
  { return Dval(i, Dget(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)); }

/****************************************************************/

inline double_complex E0(const double p1, const double p2, 
			 const double p3, const double p4, const double p5, 
			 const double p1p2, const double p2p3, 
			 const double p3p4, const double p4p5, 
			 const double p5p1, 
			 const double m1, const double m2, 
			 const double m3, const double m4, const double m5)
{
  DOUBLE_COMPLEX result;
  e0sub_(&result, &p1, &p2, &p3, &p4, &p5, &p1p2, &p2p3, &p3p4, &p4p5, &p5p1, &m1, &m2, &m3, &m4, &m5);
  return ToComplex(result);
}

inline long Eget(const double p1, const double p2, 
		 const double p3, const double p4, 
		 const double p5, 
		 const double p1p2, const double p2p3, 
		 const double p3p4, const double p4p5, 
		 const double p5p1, 
		 const double m1, const double m2, 
		 const double m3, const double m4, const double m5)
{
  return eget_(&p1, &p2, &p3, &p4, &p5, &p1p2, &p2p3, &p3p4, &p4p5, &p5p1, &m1, &m2, &m3, &m4, &m5);
}

inline DOUBLE_COMPLEX *Ecache(const long integral)
  { return &ltvars_.cache[6][integral]; }

inline DOUBLE_COMPLEX *EcacheC(const long integral)
  { return &ltvars_.cache[7][integral]; }

inline double_complex Eval(const long i, const long integral)
  { return ToComplex(Ecache(integral)[i]); }

inline double_complex EvalC(const long i, const long integral)
  { return ToComplex(EcacheC(integral)[i]); }

inline double_complex E0i(const long i, 
			  const double p1, const double p2, 
			  const double p3, const double p4, 
			  const double p5, 
			  const double p1p2, const double p2p3, 
			  const double p3p4, const double p4p5, 
			  const double p5p1, 
			  const double m1, const double m2, 
			  const double m3, const double m4, const double m5)
{ return Eval(i, Eget(p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)); }

/****************************************************************/

inline double_complex Li2(const double x)
{
  DOUBLE_COMPLEX result;
  li2sub_(&result, &x);
  return ToComplex(result);
}

/****************************************************************/

#define clearcache clearcache_
#define markcache markcache_
#define restorecache restorecache_

inline void setmudim(const double mudim)
{
  ffregul_.mudim = mudim;
  clearcache();
}

inline double getmudim() { return ffregul_.mudim; }


inline void setdelta(const double delta)
{
  ffregul_.delta = delta;
  clearcache();
}

inline double getdelta() { return ffregul_.delta; }


inline void setlambda(const double lambda)
{
  ffregul_.lambda = lambda;
  clearcache();
}

inline double getlambda() { return ffregul_.lambda; }


inline void setmaxdev(const double maxdev)
{
  ltvars_.maxdev = maxdev;
}

inline double getmaxdev() { return ltvars_.maxdev; }


inline void setwarndigits(const long warndigits)
{
  ltvars_.warndigits = warndigits;
}

inline long getwarndigits() { return ltvars_.warndigits; }


inline void seterrdigits(const long errdigits)
{
  ltvars_.errdigits = errdigits;
}

inline long geterrdigits() { return ltvars_.errdigits; }


inline void setversionkey(const long versionkey)
{
  ltvars_.versionkey = versionkey;
  clearcache();
}

inline long getversionkey() { return ltvars_.versionkey; }


inline void setdebugkey(const long debugkey)
{
  ltvars_.debugkey = debugkey;
}

inline long getdebugkey() { return ltvars_.debugkey; }


inline void setdebugrange(const long debugfrom, const long debugto)
{
  ltvars_.debugfrom = debugfrom;
  ltvars_.debugto = debugto;
}

  }
}

#endif

