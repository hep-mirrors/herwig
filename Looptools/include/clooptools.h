/* -*- C++ -*-
  clooptools.h
  the C++ header file with all definitions for LoopTools
  this file is part of LoopTools
  last modified 22 Jul 04 th

  Major modifiactions for Herwig 2005-06-09 D.Grellscheid
  Changes to reduce DOXYGEN warnings 2007/02/05 PR
  Output redirected to log file. 2007-07-11 dgrell
*/

#ifndef HERWIG_CLOOPTOOLS_H
#define HERWIG_CLOOPTOOLS_H

#include <complex>
typedef std::complex<double> double_complex;
#include <string>

// don't know why that line is here. The function is not declared anywhere.
// #define cachelookup cachelookup_

// =========== declarations of Fortran functions ====================
/**
 *  Definition of a struct to represent complex numbers
 */
struct dcomplex 
{
  /**
   *  The real part
   */
  double r;

  /**
   *  The complex part
   */ 
  double i;
};

extern "C" {

  void a0sub_(dcomplex *, const double *);
  void ca0sub_(dcomplex *, const dcomplex *);

  void b0sub_(dcomplex *, const double *, const double *, const double *);
  void cb0sub_(dcomplex *, const dcomplex *, 
	       const dcomplex *, const dcomplex *);
  void db0sub_(dcomplex *, const double *, const double *, const double *);
  void cdb0sub_(dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *);
  void b1sub_(dcomplex *, const double *, const double *, const double *);
  void cb1sub_(dcomplex *, const dcomplex *, 
	       const dcomplex *, const dcomplex *);
  void db1sub_(dcomplex *, const double *, const double *, const double *);
  void cdb1sub_(dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *);
  void b00sub_(dcomplex *, const double *, const double *, const double *);
  void cb00sub_(dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *);
  void db00sub_(dcomplex *, const double *, const double *, const double *);
  void cdb00sub_(dcomplex *, const dcomplex *, 
		 const dcomplex *, const dcomplex *);
  void b11sub_(dcomplex *, const double *, const double *, const double *);
  void cb11sub_(dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *);
  void db11sub_(dcomplex *, const double *, 
		const double *, const double *);
  void cdb11sub_(dcomplex *, const dcomplex *, 
		 const dcomplex *, const dcomplex *);

  void c0sub_(dcomplex *,
	      const double *, const double *, const double *,
	      const double *, const double *, const double *);
  void cc0sub_(dcomplex *,
	       const dcomplex *, const dcomplex *, const dcomplex *,
	       const dcomplex *, const dcomplex *, const dcomplex *);
  void c0isub_(dcomplex *, const int *,
	       const double *, const double *, const double *,
	       const double *, const double *, const double *);
  void cc0isub_(dcomplex *, const int *,
		const dcomplex *, const dcomplex *, const dcomplex *,
		const dcomplex *, const dcomplex *, const dcomplex *);
  int cget_(
	    const double *, const double *, const double *,
	    const double *, const double *, const double *);
  int ccget_(
	     const dcomplex *, const dcomplex *, const dcomplex *,
	     const dcomplex *, const dcomplex *, const dcomplex *);

  void d0sub_(dcomplex *,
	      const double *, const double *, const double *, const double *,
	      const double *, const double *,
	      const double *, const double *, const double *, const double *);
  void cd0sub_(dcomplex *,
	       const dcomplex *, const dcomplex *, 
	       const dcomplex *, const dcomplex *,
	       const dcomplex *, const dcomplex *,
	       const dcomplex *, const dcomplex *, 
	       const dcomplex *, const dcomplex *);
  void d0isub_(dcomplex *, const int *,
	       const double *, const double *, 
	       const double *, const double *,
	       const double *, const double *,
	       const double *, const double *, 
	       const double *, const double *);
  void cd0isub_(dcomplex *, const int *,
		const dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *,
		const dcomplex *, const dcomplex *,
		const dcomplex *, const dcomplex *, 
		const dcomplex *, const dcomplex *);
  int dget_(
	    const double *, const double *, 
	    const double *, const double *,
	    const double *, const double *,
	    const double *, const double *, 
	    const double *, const double *);
  int cdget_(
	     const dcomplex *, const dcomplex *, 
	     const dcomplex *, const dcomplex *,
	     const dcomplex *, const dcomplex *,
	     const dcomplex *, const dcomplex *, 
	     const dcomplex *, const dcomplex *);

  void ffini_();
  void ffexi_();
  

  void setmudim_(const double *);
  double getmudim_();
  void setdelta_(const double *);
  double getdelta_();
  void setlambda_(const double *);
  double getlambda_();

  void setcachelast_(const dcomplex *, const int *);
  int getcachelast_(const dcomplex *);

  extern dcomplex cbase_[], ccbase_[], dbase_[], cdbase_[];

} // extern "C"

// ============= end of Fortran function declarations ============


namespace Herwig {
  namespace Looptools {

    /**
     * aliases for 3pt-function coefficient IDs
     */
    enum CType {cc0=1,cc1=2,cc2=3,
		cc00=4,cc11=5,cc12=6,cc22=7,
		cc001=8,cc002=9,cc111=10,cc112=11,cc122=12,cc222=13};

    /**
     * aliases for 4pt-function coefficient IDs
     */
    enum DType {dd0=1,dd1=2,dd2=3,dd3=4,
		dd00=5,dd11=6,dd12=7,dd13=8,dd22=9,dd23=10,dd33=11,
		dd001=12,dd002=13,dd003=14,dd111=15,dd112=16,dd113=17,
		dd122=18,dd123=19,dd133=20,dd222=21,dd223=22,dd233=23,
		dd333=24,
		dd0000=25,dd0011=26,dd0012=27,dd0013=28,dd0022=29,dd0023=30,
		dd0033=31,dd1111=32,dd1112=33,dd1113=34,dd1122=35,dd1123=36,
		dd1133=37,dd1222=38,dd1223=39,dd1233=40,dd1333=41,dd2222=42,
		dd2223=43,dd2233=44,dd2333=45,dd3333=46};

    // ========== C++ wrappers for Fortran functions =============
    // for some reason not all functions have a wrapper

    /**
     *  Cache for \f$C\f$ functions
     */
    inline double_complex Ccache(int pos) {
      return double_complex(cbase_[pos - 1].r, cbase_[pos - 1].i);
    }

    /**
     *  Cache for \f$C\f$ functions
     */
    inline double_complex CCcache(int pos) {
      return double_complex(ccbase_[pos - 1].r, ccbase_[pos - 1].i);
    }

    /**
     *  Cache for \f$D\f$ functions
     */
    inline double_complex Dcache(int pos) {
      return double_complex(dbase_[pos - 1].r, dbase_[pos - 1].i);
    }

    /**
     *  Cache for \f$D\f$ functions
     */
    inline double_complex CDcache(int pos) {
      return double_complex(cdbase_[pos - 1].r, cdbase_[pos - 1].i);
    }


    /**
     *  \f$C\f$ functions
     */
    inline double_complex Cval(CType id, int pos) {
      return Ccache(pos + id);
    }

    /**
     *  \f$D\f$ functions
     */
    inline double_complex Dval(DType id, int pos) {
      return Dcache(pos + id);
    }

    /**
     *  The \f$A_0\f$ function
     */
    inline double_complex A0(const double m)
    {
      dcomplex result;

      a0sub_(&result, &m);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$B_0\f$ function
     */
    inline double_complex B0(const double p,
			     const double m1, const double m2)
    {
      dcomplex result;

      b0sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The derivative of the \f$B_0\f$ function
     */
    inline double_complex DB0(const double p,
			      const double m1, const double m2)
    {
      dcomplex result;

      db0sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$B_1\f$ function
     */
    inline double_complex B1(const double p,
			     const double m1, const double m2)
    {
      dcomplex result;

      b1sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The derivative of the \f$B_1\f$ function
     */
    inline double_complex DB1(const double p,
			      const double m1, const double m2)
    {
      dcomplex result;

      db1sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$B_{00}\f$ function
     */
    inline double_complex B00(const double p,
			      const double m1, const double m2)
    {
      dcomplex result;

      b00sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The derivative of the \f$B_{00}\f$ function
     */
    inline double_complex DB00(const double p,
			       const double m1, const double m2)
    {
      dcomplex result;

      db00sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$B_{11}\f$ function
     */
    inline double_complex B11(const double p,
			      const double m1, const double m2)
    {
      dcomplex result;

      b11sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The derivative of the \f$B_{11}\f$ function
     */
    inline double_complex DB11(const double p,
			       const double m1, const double m2)
    {
      dcomplex result;

      db11sub_(&result, &p, &m1, &m2);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$C_{0}\f$ function
     */
    inline double_complex C0(const double p1,
			     const double p2, const double p1p2,
			     const double m1, const double m2, const double m3)
    {
      dcomplex result;

      c0sub_(&result, &p1, &p2, &p1p2, &m1, &m2, &m3);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$C\f$ functions
     */
    inline double_complex C0i(const int id, const double p1,
			      const double p2, const double p1p2,
			      const double m1, const double m2, const double m3)
    {
      dcomplex result;

      c0isub_(&result, &id, &p1, &p2, &p1p2, &m1, &m2, &m3);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$C\f$ functions
     */
    inline int Cget(const double p1,
		    const double p2, const double p1p2,
		    const double m1, const double m2, const double m3)
    {
      return cget_(&p1, &p2, &p1p2, &m1, &m2, &m3);
    }

    /**
     *  The \f$D_0\f$ function
     */
    inline double_complex D0(const double p1,
			     const double p2, const double p3, const double p4,
			     const double p1p2, const double p2p3,
			     const double m1, const double m2, const double m3, const double m4)
    {
      dcomplex result;

      d0sub_(&result, &p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$D\f$ functions
     */
    inline double_complex D0i(const int id,
			      const double p1, const double p2, const double p3, const double p4,
			      const double p1p2, const double p2p3,
			      const double m1, const double m2, const double m3, const double m4)
    {
      dcomplex result;

      d0isub_(&result, &id, &p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
      return double_complex(result.r, result.i);
    }

    /**
     *  The \f$D\f$ functions
     */
    inline int Dget(const double p1,
		    const double p2, const double p3, const double p4,
		    const double p1p2, const double p2p3,
		    const double m1, const double m2, const double m3, const double m4)
    {
      return dget_(&p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
    }

    /**
     *  Looptools initialisation
     */
    void ffini(std::string logfilename = std::string("Looptools.log"));

    /**
     *  Looptools termination
     */
    void ffexi(std::string logfilename = std::string("Looptools.log"));

    /**
     *  Set \f$\mu\f$
     */
    inline void setmudim(const double newmudim) { setmudim_(&newmudim); }

    /**
     *  Get \f$\mu\f$
     */
    inline double getmudim() { return getmudim_(); }

    /**
     *  Set \f$\delta\f$
     */
    inline void setdelta(const double newdelta) { setdelta_(&newdelta); }

    /**
     *  Get \f$\delta\f$
     */
    inline double getdelta() { return getdelta_(); }

    /**
     *  Set \f$\lambda\f$
     */
    inline void setlambda(const double newlambda) { setlambda_(&newlambda); }

    /**
     *  Get \f$\lambda\f$
     */
    inline double getlambda() {return getlambda_(); }

    /**
     *  Set the cache
     */
    inline void setcachelast(const dcomplex *buffer, const int offset)
    {
      setcachelast_(buffer, &offset);
    }

    /**
     *  Get the cache
     */
    inline int getcachelast(const dcomplex * buffer) 
    {
      return getcachelast_(buffer);
    }

    // ========== end of C++ wrappers for Fortran functions =============

  } // namespace Looptools
} // namespace Herwig


#endif

