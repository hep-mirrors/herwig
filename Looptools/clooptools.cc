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

  extern memindex aget_(AARGS(_lt_Fr_));
  extern memindex agetc_(AARGS(_lt_Fc_));
  extern void aput_(COMPLEX *res, AARGS(_lt_Fr_));
  extern void aputc_(COMPLEX *res, AARGS(_lt_Fc_));
  extern void aputnocache_(COMPLEX *res, AARGS(_lt_Fr_));
  extern void aputnocachec_(COMPLEX *res, AARGS(_lt_Fc_));

  extern memindex bget_(BARGS(_lt_Fr_));
  extern memindex bgetc_(BARGS(_lt_Fc_));
  extern void bput_(COMPLEX *res, BARGS(_lt_Fr_));
  extern void bputc_(COMPLEX *res, BARGS(_lt_Fc_));
  extern void bputnocache_(COMPLEX *res, BARGS(_lt_Fr_));
  extern void bputnocachec_(COMPLEX *res, BARGS(_lt_Fc_));

  extern memindex cget_(CARGS(_lt_Fr_));
  extern memindex cgetc_(CARGS(_lt_Fc_));
  extern void cput_(COMPLEX *res, CARGS(_lt_Fr_));
  extern void cputc_(COMPLEX *res, CARGS(_lt_Fc_));
  extern void c0nocache_(COMPLEX *res, CARGS(_lt_Fr_));
  extern void c0nocachec_(COMPLEX *res, CARGS(_lt_Fc_));

  extern memindex dget_(DARGS(_lt_Fr_));
  extern memindex dgetc_(DARGS(_lt_Fc_));
  extern void dput_(COMPLEX *res, DARGS(_lt_Fr_));
  extern void dputc_(COMPLEX *res, DARGS(_lt_Fc_));
  extern void d0nocache_(COMPLEX *res, DARGS(_lt_Fr_));
  extern void d0nocachec_(COMPLEX *res, DARGS(_lt_Fc_));

  extern memindex eget_(EARGS(_lt_Fr_));
  extern memindex egetc_(EARGS(_lt_Fc_));
  extern void eput_(COMPLEX *res, EARGS(_lt_Fr_));
  extern void eputc_(COMPLEX *res, EARGS(_lt_Fc_));
  extern void e0nocache_(COMPLEX *res, EARGS(_lt_Fr_));
  extern void e0nocachec_(COMPLEX *res, EARGS(_lt_Fc_));

  extern void li2sub_(COMPLEX *res, _lt_Fr_(x));
  extern void li2csub_(COMPLEX *res, _lt_Fc_(x));

  extern void li2omxsub_(COMPLEX *res, _lt_Fr_(x));
  extern void li2omxcsub_(COMPLEX *res, _lt_Fc_(x));

  extern void ltini_(void);
  extern void ltexi_(void);

  extern void clearcache_(void);
  extern void markcache_(void);
  extern void restorecache_(void);


  extern struct {         /* MUST match common block ltvars in lt.h! */
    COMPLEX cache[10][2];
    COMPLEX savedptr[10];
    REAL maxdev;
    INTEGER epsi, warndigits, errdigits;
    INTEGER serial, versionkey;
    INTEGER debugkey, debugfrom, debugto;
  } ltvars_;
  
  extern struct {         /* MUST match common block ltcache in lt.h! */
    INTEGER cmpbits;
  } ltcache_;
  
  extern struct {         /* MUST match common block ltregul in ff.h! */
    REAL mudim, im_mudim, delta, uvdiv, lambda, minmass;
    REAL diffeps, zeroeps;
  } ltregul_;

}

#define CACHEPTR(n,i) &ltvars_.cache[n][i]

#define EPSINDEX(i) i+ltvars_.epsi

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


    ComplexType ToComplex(COMPLEX c) {
      return ComplexType(c.re, c.im);
    }

    /****************************************************************/

    memindex Aget(AARGS(_lt_Cr_))
    {
      return aget_(AARGS(_lt_Frp_));
    }

    memindex AgetC(AARGS(_lt_Cc_))
    {
      return agetc_(AARGS(_lt_Fcp_));
    }

    void Aput(ComplexType *res, AARGS(_lt_Cr_))
    {
      aput_(_lt_Fap_(res), AARGS(_lt_Frp_));
    }

    void AputC(ComplexType *res, AARGS(_lt_Cc_))
    {
      aputc_(_lt_Fap_(res), AARGS(_lt_Fcp_));
    }

    void Aputnocache(ComplexType *res, AARGS(_lt_Cr_))
    {
      aputnocache_(_lt_Fap_(res), AARGS(_lt_Frp_));
    }
    
    void AputnocacheC(ComplexType *res, AARGS(_lt_Cc_))
    {
      aputnocachec_(_lt_Fap_(res), AARGS(_lt_Fcp_));
    }
    
    COMPLEX *Acache(const memindex integral)
      { return CACHEPTR(0,integral); }
    
    COMPLEX *AcacheC(const memindex integral)
      { return CACHEPTR(1,integral); }
    
    ComplexType Aval(const int i, const memindex integral)
      { return ToComplex(Acache(integral)[i]); }
    
    ComplexType AvalC(const int i, const memindex integral)
      { return ToComplex(AcacheC(integral)[i]); }
    
    ComplexType A0i(const int i, AARGS(_lt_Cr_))
      { return Aval(EPSINDEX(i), Aget(AARGS(_lt_Id_))); }
    
    ComplexType A0iC(const int i, AARGS(_lt_Cc_))
      { return AvalC(EPSINDEX(i), AgetC(AARGS(_lt_Id_))); }
    
    ComplexType A0(AARGS(_lt_Cr_))
      { return A0i(aa0, AARGS(_lt_Id_)); }
    ComplexType A00(AARGS(_lt_Cr_))
      { return A0i(aa00, AARGS(_lt_Id_)); }
    
    ComplexType A0C(AARGS(_lt_Cc_))
      { return A0iC(aa0, AARGS(_lt_Id_)); }
    ComplexType A00C(AARGS(_lt_Cc_))
      { return A0iC(aa00, AARGS(_lt_Id_)); }


    /****************************************************************/

    memindex Bget(BARGS(_lt_Cr_))
    {
      return bget_(BARGS(_lt_Frp_));
    }

    memindex BgetC(BARGS(_lt_Cc_))
    {
      return bgetc_(BARGS(_lt_Fcp_));
    }

    void Bput(ComplexType *res, BARGS(_lt_Cr_))
    {
      bput_(_lt_Fap_(res), BARGS(_lt_Frp_));
    }

    void BputC(ComplexType *res, BARGS(_lt_Cc_))
    {
      bputc_(_lt_Fap_(res), BARGS(_lt_Fcp_));
    }

    void Bputnocache(ComplexType *res, BARGS(_lt_Cr_))
    {
      bputnocache_(_lt_Fap_(res), BARGS(_lt_Frp_));
    }
    
    void BputnocacheC(ComplexType *res, BARGS(_lt_Cc_))
    {
      bputnocachec_(_lt_Fap_(res), BARGS(_lt_Fcp_));
    }

    COMPLEX *Bcache(const memindex integral)
    { return CACHEPTR(2,integral); }

    COMPLEX *BcacheC(const memindex integral)
    { return CACHEPTR(3,integral); }

    ComplexType Bval(const int i, const memindex integral)
    { return ToComplex(Bcache(integral)[i]); }

    ComplexType BvalC(const int i, const memindex integral)
    { return ToComplex(BcacheC(integral)[i]); }

    ComplexType B0i(const int i, BARGS(_lt_Cr_))
    { return Bval(EPSINDEX(i), Bget(BARGS(_lt_Id_))); }

    ComplexType B0iC(const int i, BARGS(_lt_Cc_))
    { return BvalC(EPSINDEX(i), BgetC(BARGS(_lt_Id_))); }

    ComplexType B0(BARGS(_lt_Cr_))
    { return B0i(bb0, BARGS(_lt_Id_)); }
    ComplexType B1(BARGS(_lt_Cr_))
    { return B0i(bb1, BARGS(_lt_Id_)); }
    ComplexType B00(BARGS(_lt_Cr_))
    { return B0i(bb00, BARGS(_lt_Id_)); }
    ComplexType B11(BARGS(_lt_Cr_))
    { return B0i(bb11, BARGS(_lt_Id_)); }
    ComplexType B001(BARGS(_lt_Cr_))
    { return B0i(bb001, BARGS(_lt_Id_)); }
    ComplexType B111(BARGS(_lt_Cr_))
    { return B0i(bb111, BARGS(_lt_Id_)); }
    ComplexType DB0(BARGS(_lt_Cr_))
    { return B0i(dbb0, BARGS(_lt_Id_)); }
    ComplexType DB1(BARGS(_lt_Cr_))
    { return B0i(dbb1, BARGS(_lt_Id_)); }
    ComplexType DB00(BARGS(_lt_Cr_))
    { return B0i(dbb00, BARGS(_lt_Id_)); }
    ComplexType DB11(BARGS(_lt_Cr_))
    { return B0i(dbb11, BARGS(_lt_Id_)); }

    ComplexType B0C(BARGS(_lt_Cc_))
    { return B0iC(bb0, BARGS(_lt_Id_)); }
    ComplexType B1C(BARGS(_lt_Cc_))
    { return B0iC(bb1, BARGS(_lt_Id_)); }
    ComplexType B00C(BARGS(_lt_Cc_))
    { return B0iC(bb00, BARGS(_lt_Id_)); }
    ComplexType B11C(BARGS(_lt_Cc_))
    { return B0iC(bb11, BARGS(_lt_Id_)); }
    ComplexType B001C(BARGS(_lt_Cc_))
    { return B0iC(bb001, BARGS(_lt_Id_)); }
    ComplexType B111C(BARGS(_lt_Cc_))
    { return B0iC(bb111, BARGS(_lt_Id_)); }
    ComplexType DB0C(BARGS(_lt_Cc_))
    { return B0iC(dbb0, BARGS(_lt_Id_)); }
    ComplexType DB1C(BARGS(_lt_Cc_))
    { return B0iC(dbb1, BARGS(_lt_Id_)); }
    ComplexType DB00C(BARGS(_lt_Cc_))
    { return B0iC(dbb00, BARGS(_lt_Id_)); }
    ComplexType DB11C(BARGS(_lt_Cc_))
    { return B0iC(dbb11, BARGS(_lt_Id_)); }

    /****************************************************************/

    memindex Cget(CARGS(_lt_Cr_))
    {
      return cget_(CARGS(_lt_Frp_)); 
    }

    memindex CgetC(CARGS(_lt_Cc_))
    {
      return cgetc_(CARGS(_lt_Fcp_));
    }

    void Cput(ComplexType *res, CARGS(_lt_Cr_))
    {
      cput_(_lt_Fap_(res), CARGS(_lt_Frp_));
    }

    void CputC(ComplexType *res, CARGS(_lt_Cc_))
    {
      cputc_(_lt_Fap_(res), CARGS(_lt_Fcp_));
    }

    void C0nocache(ComplexType *res, CARGS(_lt_Cr_))
    {
      c0nocache_(_lt_Fap_(res), CARGS(_lt_Frp_));
    }
    
    void C0nocacheC(ComplexType *res, CARGS(_lt_Cc_))
    {
      c0nocachec_(_lt_Fap_(res), CARGS(_lt_Fcp_));
    }

    COMPLEX *Ccache(const memindex integral)
    { return CACHEPTR(4,integral); }

    COMPLEX *CcacheC(const memindex integral)
    { return CACHEPTR(5,integral); }

    ComplexType Cval(const int i, const memindex integral)
    { return ToComplex(Ccache(integral)[i]); }

    ComplexType CvalC(const int i, const memindex integral)
    { return ToComplex(CcacheC(integral)[i]); }

    ComplexType C0i(const int i, CARGS(_lt_Cr_))
    { return Cval(EPSINDEX(i), Cget(CARGS(_lt_Id_))); }

    ComplexType C0iC(const int i, CARGS(_lt_Cc_))
    { return CvalC(EPSINDEX(i), CgetC(CARGS(_lt_Id_))); }

    ComplexType C0(CARGS(_lt_Cr_))
    { return C0i(cc0, CARGS(_lt_Id_)); }
    ComplexType C0C(CARGS(_lt_Cc_))
    { return C0iC(cc0, CARGS(_lt_Id_)); }

    /****************************************************************/

    memindex Dget(DARGS(_lt_Cr_))
    {
      return dget_(DARGS(_lt_Frp_)); 
    }

    memindex DgetC(DARGS(_lt_Cc_))
    {
      return dgetc_(DARGS(_lt_Fcp_));
    }

    void Dput(ComplexType *res, DARGS(_lt_Cr_))
    {
      dput_(_lt_Fap_(res), DARGS(_lt_Frp_));
    }

    void DputC(ComplexType *res, DARGS(_lt_Cc_))
    {
      dputc_(_lt_Fap_(res), DARGS(_lt_Fcp_));
    }

    void D0nocache(ComplexType *res, DARGS(_lt_Cr_))
    {
      d0nocache_(_lt_Fap_(res), DARGS(_lt_Frp_));
    }
    
    void D0nocacheC(ComplexType *res, DARGS(_lt_Cc_))
    {
      d0nocachec_(_lt_Fap_(res), DARGS(_lt_Fcp_));
    }

    COMPLEX *Dcache(const memindex integral)
    { return CACHEPTR(6,integral); }

    COMPLEX *DcacheC(const memindex integral)
    { return CACHEPTR(7,integral); }

    ComplexType Dval(const int i, const memindex integral)
    { return ToComplex(Dcache(integral)[i]); }

    ComplexType DvalC(const int i, const memindex integral)
    { return ToComplex(DcacheC(integral)[i]); }

    ComplexType D0i(const int i, DARGS(_lt_Cr_))
    { return Dval(EPSINDEX(i), Dget(DARGS(_lt_Id_))); }

    ComplexType D0iC(const int i, DARGS(_lt_Cc_))
    { return DvalC(EPSINDEX(i), DgetC(DARGS(_lt_Id_))); }

    ComplexType D0(DARGS(_lt_Cr_))
    { return D0i(dd0, DARGS(_lt_Id_)); }
    ComplexType D0C(DARGS(_lt_Cc_))
    { return D0iC(dd0, DARGS(_lt_Id_)); }

    /****************************************************************/

    memindex Eget(EARGS(_lt_Cr_))
    {
      return eget_(EARGS(_lt_Frp_)); 
    }

    memindex EgetC(EARGS(_lt_Cc_))
    {
      return egetc_(EARGS(_lt_Fcp_));
    }

    void Eput(ComplexType *res, EARGS(_lt_Cr_))
    {
      eput_(_lt_Fap_(res), EARGS(_lt_Frp_));
    }

    void EputC(ComplexType *res, EARGS(_lt_Cc_))
    {
      eputc_(_lt_Fap_(res), EARGS(_lt_Fcp_));
    }

    void E0nocache(ComplexType *res, EARGS(_lt_Cr_))
    {
      e0nocache_(_lt_Fap_(res), EARGS(_lt_Frp_));
    }
    
    void E0nocacheC(ComplexType *res, EARGS(_lt_Cc_))
    {
      e0nocachec_(_lt_Fap_(res), EARGS(_lt_Fcp_));
    }

    COMPLEX *Ecache(const memindex integral)
    { return CACHEPTR(8,integral); }

    COMPLEX *EcacheC(const memindex integral)
    { return CACHEPTR(9,integral); }

    ComplexType Eval(const int i, const memindex integral)
    { return ToComplex(Ecache(integral)[i]); }

    ComplexType EvalC(const int i, const memindex integral)
    { return ToComplex(EcacheC(integral)[i]); }

    ComplexType E0i(const int i, EARGS(_lt_Cr_))
    { return Eval(EPSINDEX(i), Eget(EARGS(_lt_Id_))); }

    ComplexType E0iC(const int i, EARGS(_lt_Cc_))
    { return EvalC(EPSINDEX(i), EgetC(EARGS(_lt_Id_))); }

    ComplexType E0(EARGS(_lt_Cr_))
    { return E0i(ee0, EARGS(_lt_Id_)); }
    ComplexType E0C(EARGS(_lt_Cc_))
    { return E0iC(ee0, EARGS(_lt_Id_)); }

    /****************************************************************/

    ComplexType Li2(_lt_Cr_(x))
    {
      COMPLEX res;
      li2sub_(&res, _lt_Frp_(x));
      return ToComplex(res);
    }

    ComplexType Li2C(_lt_Cc_(x))
    {
      COMPLEX res;
      li2csub_(&res, _lt_Fcp_(x));
      return ToComplex(res);
    }

    ComplexType Li2omx(_lt_Cr_(x))
    {
      COMPLEX res;
      li2omxsub_(&res, _lt_Frp_(x));
      return ToComplex(res);
    }

    ComplexType Li2omxC(_lt_Cc_(x))
    {
      COMPLEX res;
      li2omxcsub_(&res, _lt_Fcp_(x));
      return ToComplex(res);
    }

    /****************************************************************/

    void clearcache()   { clearcache_(); }
    void markcache()    { markcache_(); }
    void restorecache() { restorecache_(); }

#define _lt_set_(v) \
  ltregul_.v = v

#define _lt_setx_(v) \
  if( fabs(ltregul_.v - v) > (ltregul_.diffeps) ) \
    clearcache(); \
  ltregul_.v = v

#define _lt_get_(v) \
  return ltregul_.v


    void setmudim(cRealType mudim) { _lt_setx_(mudim); }

    RealType getmudim() { _lt_get_(mudim); }


    void setdelta(cRealType delta) { _lt_setx_(delta); }

    RealType getdelta() { _lt_get_(delta); }


    void setuvdiv(cRealType uvdiv) { _lt_setx_(uvdiv); }

    RealType getuvdiv() { _lt_get_(uvdiv); }


    void setlambda(RealType lambda) {
      int epsi = 0;
      if( lambda <= 0 ) {
        if( lambda != 0 && lambda != -1 && lambda != -2 ) {
          fprintf(stderr, "illegal value for lambda\n");
          lambda = 0;
        }
        epsi = -(int)lambda;
      }
      _lt_setx_(lambda);
      ltvars_.epsi = epsi;
    }

    RealType getlambda() { _lt_get_(lambda); }

    int getepsi() { 
      return ltvars_.epsi;
    }


    void setminmass(cRealType minmass) { _lt_setx_(minmass); }

    RealType getminmass() { _lt_get_(minmass); }


    void setmaxdev(cRealType maxdev) {
      ltvars_.maxdev = maxdev;
    }

    RealType getmaxdev() { return ltvars_.maxdev; }


    void setwarndigits(const int warndigits) {
      ltvars_.warndigits = warndigits;
    }

    int getwarndigits() { return ltvars_.warndigits; }


    void seterrdigits(const int errdigits) {
      ltvars_.errdigits = errdigits;
    }

    int geterrdigits() { return ltvars_.errdigits; }


    void setversionkey(const int versionkey) {
      ltvars_.versionkey = versionkey;
      clearcache();
    }

    int getversionkey() { return ltvars_.versionkey; }


    void setdebugkey(const int debugkey) {
      ltvars_.debugkey = debugkey;
    }

    int getdebugkey() { return ltvars_.debugkey; }


    void setdebugrange(const int debugfrom, const int debugto)
    {
      ltvars_.debugfrom = debugfrom;
      ltvars_.debugto = debugto;
    }


    void setcmpbits(const int cmpbits)
    {
      ltcache_.cmpbits = cmpbits;
    }

    int getcmpbits() { return ltcache_.cmpbits; }


    void setdiffeps(cRealType diffeps) { _lt_set_(diffeps); }

    RealType getdiffeps() { _lt_get_(diffeps); }


    void setzeroeps(cRealType zeroeps) { _lt_set_(zeroeps); }

    RealType getzeroeps() { _lt_get_(zeroeps); }

  } // namespace Looptools
} // namespace Herwig
