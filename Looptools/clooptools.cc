/* -*- C++ -*-
   clooptools.cc
   the C++ file with the definitions for fortran IO redirection
   Output redirected to log file. 2007-07-18 dgrell

   Definitions moved here from clooptools.h 2011-01-21 dgrell
*/

#include "Herwig/Looptools/clooptools.h"

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>
#include <cassert>
#include <string>

#include "ThePEG/Repository/CurrentGenerator.h"


extern "C" {

#include <unistd.h>


  extern void a0sub_(DOUBLE_COMPLEX *result, AARGS(_Fr_));
  extern void a0subc_(DOUBLE_COMPLEX *result, AARGS(_Fc_));
  extern void a00sub_(DOUBLE_COMPLEX *result, AARGS(_Fr_));
  extern void a00subc_(DOUBLE_COMPLEX *result, AARGS(_Fc_));

  extern long bget_(BARGS(_Fr_));
  extern long bgetc_(BARGS(_Fc_));

  extern void c0sub_(DOUBLE_COMPLEX *result, CARGS(_Fr_));
  extern void c0subc_(DOUBLE_COMPLEX *result, CARGS(_Fc_));
  extern long cget_(CARGS(_Fr_));
  extern long cgetc_(CARGS(_Fc_));

  extern void d0sub_(DOUBLE_COMPLEX *result, DARGS(_Fr_));
  extern void d0subc_(DOUBLE_COMPLEX *result, DARGS(_Fc_));
  extern long dget_(DARGS(_Fr_));
  extern long dgetc_(DARGS(_Fc_));

  extern void e0sub_(DOUBLE_COMPLEX *result, EARGS(_Fr_));
  extern void e0subc_(DOUBLE_COMPLEX *result, EARGS(_Fc_));
  extern long eget_(EARGS(_Fr_));
  extern long egetc_(EARGS(_Fc_));

  extern void li2sub_(DOUBLE_COMPLEX *result, const double *x);
  extern void li2csub_(DOUBLE_COMPLEX *result, CDOUBLE_COMPLEX *x);

  extern void ltini_(void);
  extern void ltexi_(void);

  extern void clearcache_(void);
  extern void markcache_(void);
  extern void restorecache_(void);

  extern struct {		/* MUST match common block ltvars in lt.h! */
    DOUBLE_COMPLEX cache[8][2];
    DOUBLE_COMPLEX savedptr[8];
    double maxdev;
    long warndigits, errdigits;
    long serial, versionkey;
    long debugkey, debugfrom, debugto;
  } ltvars_;

  extern struct {		/* MUST match common block ltcache in lt.h! */
    long cmpbits;
  } ltcache_;

  extern struct {		/* MUST match common block ltregul in ff.h! */
    double mudim, delta, lambda, minmass;
  } ltregul_;

}

/****************************************************************/

namespace {

  int start_redirection(std::string logfilename) {
    if ( ! ThePEG::CurrentGenerator::isVoid() 
	 && ThePEG::CurrentGenerator::current().useStdOut() ) return -1;
    // redirect C stdout --- unix specific solution,
    // see C FAQ: http://c-faq.com/stdio/undofreopen.html
    int    fd;
    fflush(stdout);
    fd = dup(fileno(stdout));
    freopen(logfilename.c_str(), "a", stdout);
    return fd;
  }
  
  void stop_redirection(int fd) {
    if ( ! ThePEG::CurrentGenerator::isVoid() 
	 && ThePEG::CurrentGenerator::current().useStdOut() ) return;
    fflush(stdout);
    close(fileno(stdout));
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
  }
  
} // namespace

namespace Herwig {
  namespace Looptools {

    static int initcount = 0;

    void ltini(std::string logfilename) {
      assert( initcount >= 0 );
      if ( initcount == 0 ) {
	int rd = start_redirection(logfilename);
	ltini_();
	stop_redirection(rd);
      }
      ++initcount;
    }

    void ltexi(std::string logfilename) {
      assert( initcount > 0 );
      --initcount;
      if ( initcount == 0 ) {
	int rd = start_redirection(logfilename);
	ltexi_();
	stop_redirection(rd);
      }
    }


    double_complex ToComplex(DOUBLE_COMPLEX c) {
      return double_complex(c.re, c.im);
    }

    double_complex A0(AARGS(_Cr_))
    {
      DOUBLE_COMPLEX result;
      a0sub_(&result, AARGS(_Frp_));
      return ToComplex(result);
    }

    double_complex A0C(AARGS(_Cc_))
    {
      DOUBLE_COMPLEX result;
      a0subc_(&result, AARGS(_Fcp_));
      return ToComplex(result);
    }

    double_complex A00(AARGS(_Cr_))
    {
      DOUBLE_COMPLEX result;
      a00sub_(&result, AARGS(_Frp_));
      return ToComplex(result);
    }

    double_complex A00C(AARGS(_Cc_))
    {
      DOUBLE_COMPLEX result;
      a00subc_(&result, AARGS(_Fcp_));
      return ToComplex(result);
    }

    /****************************************************************/

    long Bget(BARGS(_Cr_))
    {
      return bget_(BARGS(_Frp_));
    }

    long BgetC(BARGS(_Cc_))
    {
      return bgetc_(BARGS(_Fcp_));
    }

    DOUBLE_COMPLEX *Bcache(const long integral)
    { return &ltvars_.cache[0][integral]; }

    DOUBLE_COMPLEX *BcacheC(const long integral)
    { return &ltvars_.cache[1][integral]; }

    double_complex Bval(const int i, const long integral)
    { return ToComplex(Bcache(integral)[i]); }

    double_complex BvalC(const int i, const long integral)
    { return ToComplex(BcacheC(integral)[i]); }

    double_complex B0i(const int i, BARGS(_Cr_))
    { return Bval(i, Bget(BARGS(_Id_))); }

    double_complex B0iC(const int i, BARGS(_Cc_))
    { return BvalC(i, BgetC(BARGS(_Id_))); }

    double_complex B0(BARGS(_Cr_))
    { return B0i(bb0, BARGS(_Id_)); }
    double_complex B1(BARGS(_Cr_))
    { return B0i(bb1, BARGS(_Id_)); }
    double_complex B00(BARGS(_Cr_))
    { return B0i(bb00, BARGS(_Id_)); }
    double_complex B11(BARGS(_Cr_))
    { return B0i(bb11, BARGS(_Id_)); }
    double_complex B001(BARGS(_Cr_))
    { return B0i(bb001, BARGS(_Id_)); }
    double_complex B111(BARGS(_Cr_))
    { return B0i(bb111, BARGS(_Id_)); }
    double_complex DB0(BARGS(_Cr_))
    { return B0i(dbb0, BARGS(_Id_)); }
    double_complex DB1(BARGS(_Cr_))
    { return B0i(dbb1, BARGS(_Id_)); }
    double_complex DB00(BARGS(_Cr_))
    { return B0i(dbb00, BARGS(_Id_)); }
    double_complex DB11(BARGS(_Cr_))
    { return B0i(dbb11, BARGS(_Id_)); }

    double_complex B0C(BARGS(_Cc_))
    { return B0iC(bb0, BARGS(_Id_)); }
    double_complex B1C(BARGS(_Cc_))
    { return B0iC(bb1, BARGS(_Id_)); }
    double_complex B00C(BARGS(_Cc_))
    { return B0iC(bb00, BARGS(_Id_)); }
    double_complex B11C(BARGS(_Cc_))
    { return B0iC(bb11, BARGS(_Id_)); }
    double_complex B001C(BARGS(_Cc_))
    { return B0iC(bb001, BARGS(_Id_)); }
    double_complex B111C(BARGS(_Cc_))
    { return B0iC(bb111, BARGS(_Id_)); }
    double_complex DB0C(BARGS(_Cc_))
    { return B0iC(dbb0, BARGS(_Id_)); }
    double_complex DB1C(BARGS(_Cc_))
    { return B0iC(dbb1, BARGS(_Id_)); }
    double_complex DB00C(BARGS(_Cc_))
    { return B0iC(dbb00, BARGS(_Id_)); }
    double_complex DB11C(BARGS(_Cc_))
    { return B0iC(dbb11, BARGS(_Id_)); }

    /****************************************************************/

    double_complex C0(CARGS(_Cr_))
    {
      DOUBLE_COMPLEX result;
      c0sub_(&result, CARGS(_Frp_));
      return ToComplex(result);
    }

    double_complex C0C(CARGS(_Cc_))
    {
      DOUBLE_COMPLEX result;
      c0subc_(&result, CARGS(_Fcp_));
      return ToComplex(result);
    }

    long Cget(CARGS(_Cr_))
    {
      return cget_(CARGS(_Frp_)); 
    }

    long CgetC(CARGS(_Cc_))
    {
      return cgetc_(CARGS(_Fcp_));
    }

    DOUBLE_COMPLEX *Ccache(const long integral)
    { return &ltvars_.cache[2][integral]; }

    DOUBLE_COMPLEX *CcacheC(const long integral)
    { return &ltvars_.cache[3][integral]; }

    double_complex Cval(const int i, const long integral)
    { return ToComplex(Ccache(integral)[i]); }

    double_complex CvalC(const int i, const long integral)
    { return ToComplex(CcacheC(integral)[i]); }

    double_complex C0i(const int i, CARGS(_Cr_))
    { return Cval(i, Cget(CARGS(_Id_))); }

    double_complex C0iC(const int i, CARGS(_Cc_))
    { return CvalC(i, CgetC(CARGS(_Id_))); }

    /****************************************************************/

    double_complex D0(DARGS(_Cr_))
    {
      DOUBLE_COMPLEX result;
      d0sub_(&result, DARGS(_Frp_));
      return ToComplex(result);
    }

    double_complex D0C(DARGS(_Cc_))
    {
      DOUBLE_COMPLEX result;
      d0subc_(&result, DARGS(_Fcp_));
      return ToComplex(result);
    }

    long Dget(DARGS(_Cr_))
    {
      return dget_(DARGS(_Frp_));
    }

    long DgetC(DARGS(_Cc_))
    {
      return dgetc_(DARGS(_Fcp_));
    }

    DOUBLE_COMPLEX *Dcache(const long integral)
    { return &ltvars_.cache[4][integral]; }

    DOUBLE_COMPLEX *DcacheC(const long integral)
    { return &ltvars_.cache[5][integral]; }

    double_complex Dval(const int i, const long integral)
    { return ToComplex(Dcache(integral)[i]); }

    double_complex DvalC(const int i, const long integral)
    { return ToComplex(DcacheC(integral)[i]); }

    double_complex D0i(const int i, DARGS(_Cr_))
    { return Dval(i, Dget(DARGS(_Id_))); }

    double_complex D0iC(const int i, DARGS(_Cc_))
    { return DvalC(i, DgetC(DARGS(_Id_))); }

    /****************************************************************/

    double_complex E0(EARGS(_Cr_))
    {
      DOUBLE_COMPLEX result;
      e0sub_(&result, EARGS(_Frp_));
      return ToComplex(result);
    }

    double_complex E0C(EARGS(_Cc_))
    {
      DOUBLE_COMPLEX result;
      e0subc_(&result, EARGS(_Fcp_));
      return ToComplex(result);
    }

    long Eget(EARGS(_Cr_))
    {
      return eget_(EARGS(_Frp_));
    }

    long EgetC(EARGS(_Cc_))
    {
      return egetc_(EARGS(_Fcp_));
    }

    DOUBLE_COMPLEX *Ecache(const long integral)
    { return &ltvars_.cache[6][integral]; }

    DOUBLE_COMPLEX *EcacheC(const long integral)
    { return &ltvars_.cache[7][integral]; }

    double_complex Eval(const int i, const long integral)
    { return ToComplex(Ecache(integral)[i]); }

    double_complex EvalC(const int i, const long integral)
    { return ToComplex(EcacheC(integral)[i]); }

    double_complex E0i(const int i, EARGS(_Cr_))
    { return Eval(i, Eget(EARGS(_Id_))); }

    double_complex E0iC(const int i, EARGS(_Cc_))
    { return EvalC(i, EgetC(EARGS(_Id_))); }

    /****************************************************************/

    double_complex Li2(const double x)
    {
      DOUBLE_COMPLEX result;
      li2sub_(&result, _Frp_(x));
      return ToComplex(result);
    }

    double_complex Li2C(const double_complex x)
    {
      DOUBLE_COMPLEX result;
      li2csub_(&result, _Fcp_(x));
      return ToComplex(result);
    }

    /****************************************************************/


    void setmudim(const double mudim)
    {
      ltregul_.mudim = mudim;
      clearcache();
    }

    double getmudim() { return ltregul_.mudim; }


    void setdelta(const double delta)
    {
      ltregul_.delta = delta;
      clearcache();
    }

    double getdelta() { return ltregul_.delta; }


    void setlambda(const double lambda)
    {
      ltregul_.lambda = lambda;
      clearcache();
    }

    double getlambda() { return ltregul_.lambda; }


    void setminmass(const double minmass)
    {
      ltregul_.minmass = minmass;
      clearcache();
    }

    double getminmass() { return ltregul_.minmass; }


    void setmaxdev(const double maxdev)
    {
      ltvars_.maxdev = maxdev;
    }

    double getmaxdev() { return ltvars_.maxdev; }


    void setwarndigits(const long warndigits)
    {
      ltvars_.warndigits = warndigits;
    }

    long getwarndigits() { return ltvars_.warndigits; }


    void seterrdigits(const long errdigits)
    {
      ltvars_.errdigits = errdigits;
    }

    long geterrdigits() { return ltvars_.errdigits; }


    void setversionkey(const long versionkey)
    {
      ltvars_.versionkey = versionkey;
      clearcache();
    }

    long getversionkey() { return ltvars_.versionkey; }


    void setdebugkey(const long debugkey)
    {
      ltvars_.debugkey = debugkey;
    }

    long getdebugkey() { return ltvars_.debugkey; }


    void setdebugrange(const long debugfrom, const long debugto)
    {
      ltvars_.debugfrom = debugfrom;
      ltvars_.debugto = debugto;
    }


    void setcmpbits(const long cmpbits)
    {
      ltcache_.cmpbits = cmpbits;
    }

    long getcmpbits() { return ltcache_.cmpbits; }


    void clearcache()   { clearcache_(); }
    void markcache()    { markcache_(); }
    void restorecache() { restorecache_(); }

  } // namespace Looptools
} // namespace Herwig
