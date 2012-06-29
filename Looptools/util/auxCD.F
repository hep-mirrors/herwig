* auxCD.F
* auxillary functions used by the three- and four-point integrals
* these functions are adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* they are used for double-checking the results of FF
* last modified 1 Feb 10 th

#include "defs.h"


	double complex function ln(x, s)
	implicit none
	double precision x, s

#include "ff.h"

	if( x .le. 0 ) then
#ifdef WARNINGS
	  if( s .eq. 0 ) print *, "ln: argument on cut"
#endif
	  ln = log(-x) + cI*sign(pi, s)
	else
	  ln = log(x)
	endif
	end

************************************************************************

	double complex function cln(z, s)
	implicit none
	double complex z
	double precision s

#include "ff.h"

	if( DIMAG(z) .eq. 0 .and. DBLE(z) .le. 0 ) then
#ifdef WARNINGS
	  if( s .eq. 0 ) print *, "cln: argument on cut"
#endif
	  cln = log(-z) + cI*sign(pi, s)
	else
	  cln = log(z)
	endif
	end

************************************************************************
* lnrat(x, y) = log(x - i eps) - log(y - i eps)
* original version by R.K. Ellis
* this function is hardwired for the sign of epsilon
* we must adjust the sign of x and y to get the right sign for epsilon

	double complex function lnrat(x, y)
	implicit none
	double precision x, y

#include "ff.h"

	lnrat = log(abs(x/y)) +
     &    cI*(sign(.5D0*pi, x) - sign(.5D0*pi, y))
	end

************************************************************************

	double complex function lndiv0(x, y)
	implicit none
	double precision x, y

	double precision den

	double complex lnrat
	external lnrat

	den = 1 - x/y
	if( abs(den) .lt. 1D-7 ) then
	  lndiv0 = -1 - den*(.5D0 + den/3D0)
	else
	  lndiv0 = lnrat(x, y)/den
	endif
	end

************************************************************************

	double complex function lndiv1(x, y)
	implicit none
	double precision x, y

	double precision den

	double complex lnrat
	external lnrat

	den = 1 - x/y
	if( abs(den) .lt. 1D-7 ) then
	  lndiv1 = -.5D0 - den/3D0*(1 + .75D0*den)
	else
	  lndiv1 = (lnrat(x, y)/den + 1)/den
	endif
	end

************************************************************************
* Li2omrat(x, y) = Li2(1 - (x - i eps)/(y - i eps)) for real x and y
* hence arguments are typically negative invariants
* original version by R.K. Ellis

	double complex function Li2omrat(x, y)
	implicit none
	double precision x, y

#include "ff.h"

	double precision omarg

	double complex spence, lnrat
	external spence, lnrat

	omarg = x/y
	if( omarg .lt. 0 ) then
	  Li2omrat = pi6 - spence(DCMPLX(omarg), 0D0) -
     &      log(1 - omarg)*lnrat(x, y)
	else
	  Li2omrat = spence(DCMPLX(1 - omarg), 0D0)
	endif
	end

************************************************************************
* Li2omx2 = Li2(1 - (z1 + i eps1) (z2 + i eps2)) for complex z1, z2
* for z1 z2 < 1: +Li2(1 - z1 z2)
* for z1 z2 > 1: -Li2(1 - 1/(z1 z2)) - 1/2 (ln(z1) + ln(z2))^2
* original version by R.K. Ellis

	double complex function Li2omx2(z1, s1, z2, s2)
	implicit none
	double complex z1, z2
	double precision s1, s2

#include "ff.h"

	double complex z12, l12
	double precision s12

	double complex cln, spence
	external cln, spence

	z12 = z1*z2
	if( abs(z12) .lt. eps ) then
	  Li2omx2 = 0
	else if( abs(z12 - 1) .eq. acc ) then
	  Li2omx2 = pi6
	else
	  l12 = cln(z1, s1) + cln(z2, s2)
	  s12 = sign(1D0, DBLE(z2))*s1 + sign(1D0, DBLE(z1))*s2
	  if( abs(z12) .le. 1 ) then
	    Li2omx2 = pi6 - spence(z12, s12) -
     &        l12*cln(1 - z12, -s12)
	  else
	    z12 = 1/z12
	    Li2omx2 = -pi6 + spence(z12, s12) -
     &        l12*(cln(1 - z12, -s12) + .5D0*l12)
	  endif
	endif
	end

************************************************************************
* Li2omx3 = Li2(1 - (z1 + i eps1) (z2 + i eps2)) for complex z1, z2
* for z1 z2 < 1: +Li2(1 - z1 z2)
* for z1 z2 > 1: -Li2(1 - 1/(z1 z2)) - 1/2 (ln(z1) + ln(z2))^2
* original version by R.K. Ellis

	double complex function Li2omx3(z1, s1, z2, s2, z3, s3)
	implicit none
	double complex z1, z2, z3
	double precision s1, s2, s3

#include "ff.h"

	double complex z123, l123
	double precision s123

	double complex cln, spence
	external cln, spence

	z123 = z1*z2*z3
	if( abs(DIMAG(z123)) .lt. eps )
     &    s123 = sign(1D0,
     &      DBLE(z2*z3)*s1 + DBLE(z1*z3)*s2 + DBLE(z1*z2)*s3)

	if( abs(z123) .le. 1 ) then
	  Li2omx3 = pi6 - spence(z123, s123)
	  if( abs(z123) .gt. eps .and. abs(z123 - 1) .gt. acc )
     &      Li2omx3 = Li2omx3 - cln(1 - z123, 0D0)*
     &        (cln(z1, s1) + cln(z2, s2) + cln(z3, s3))
	else
	  z123 = 1/z123
	  l123 = cln(z1, s1) + cln(z2, s2) + cln(z3, s3)
	  Li2omx3 = -pi6 + spence(z123, s123) -
     &      l123*(.5D0*l123 - cln(1 - z123, 0D0))
	endif
	end

************************************************************************
* Li2omrat2 = Li2(1 - (n1 - i eps) (n2 - i eps)/(d1 - i eps)/(d2 - i eps))
* for real n1, n2, d1, d2
* original version by R.K. Ellis

	double complex function Li2omrat2(n1, d1, n2, d2)
	implicit none
	double precision n1, d1, n2, d2

#include "ff.h"

	double precision r12
	double complex l12

	double complex lnrat, spence
	external lnrat, spence

	r12 = n1*n2/(d1*d2)
	if( r12 .lt. 1 ) then
	  Li2omrat2 = pi6 - spence(DCMPLX(r12), 0D0)
	  if( abs(r12*(1 - r12)) .gt. acc )
     &      Li2omrat2 = Li2omrat2 -
     &        (lnrat(n1, d1) + lnrat(n2, d2))*log(1 - r12)
	else
	  r12 = 1/r12
	  l12 = lnrat(n1, d1) + lnrat(n2, d2)
	  Li2omrat2 = -pi6 + spence(DCMPLX(r12), 0D0) -
     &      l12*(.5D0*l12 + log(1 - r12))
	endif
	end

************************************************************************
* original version by R.K. Ellis

	double complex function Li2rat(r1, s1, r2, s2)
	implicit none
	double complex r1, r2
	double precision s1, s2

#include "ff.h"

	double complex r12, l12

	double complex Li2omx2, cln, spence
	external Li2omx2, cln, spence

	if( abs(DIMAG(r1)) + abs(DIMAG(r2)) .lt. eps ) then
	  Li2rat = Li2omx2(r1, s1, r2, s2)
	  return
	endif

	r12 = r1*r2
	if( abs(r12) .lt. 1 ) then
	  Li2rat = pi6 - spence(r12, 0D0)
	  if( abs(r12*(1 - r12)) .gt. acc )
     &      Li2rat = Li2rat - (cln(r1, s1) + cln(r2, s2))*log(1 - r12)
	else
	  r12 = 1/r12
	  l12 = cln(r1, s1) + cln(r2, s2)
	  Li2rat = -pi6 + spence(r12, 0D0) -
     &      l12*(.5D0*l12 - log(1 - r12))
	endif
	end

************************************************************************

	double complex function spence(z, s)
	implicit none
	double complex z
	double precision s

#include "ff.h"

	double complex z1
	double precision az1

	double complex Li2series, cln
	external Li2series, cln

	z1 = 1 - z
	az1 = abs(z1)

#ifdef WARNINGS
	if( s .eq. 0 .and. 
     &    DIMAG(z) .eq. 0 .and. abs(DBLE(z1)) .lt. acc )
     &    print *, "spence: argument on cut"
#endif

	if( az1 .lt. acc ) then
	  spence = pi6
	else if( DBLE(z) .lt. .5D0 ) then
	  if( abs(z) .lt. 1 ) then
	    spence = Li2series(z, s)
	  else
	    spence = -pi6 -
     &        .5D0*cln(-z, -s)**2 - Li2series(1/z, -s)
	  endif
	else
	  if( az1 .lt. 1 ) then
	    spence = pi6 -
     &        cln(z, s)*cln(z1, -s) - Li2series(z1, -s)
	  else
	    spence = 2*pi6 +
     &        .5D0*cln(-z1, -s)**2 - cln(z, s)*cln(z1, -s) +
     &        Li2series(1/z1, s)
	  endif
	endif
	end

************************************************************************

	double complex function Li2series(z, s)
	implicit none
	double complex z
	double precision s

	double complex xm, x2, new
	integer j

	double complex cln
	external cln

* these are the even-n Bernoulli numbers, already divided by (n + 1)!
* as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 50, 2}]
	double precision b(25)
	data b /
     &    0.02777777777777777777777777777777777777777778774D0, 
     &    -0.000277777777777777777777777777777777777777777778D0, 
     &    4.72411186696900982615268329554043839758125472D-6, 
     &    -9.18577307466196355085243974132863021751910641D-8, 
     &    1.89788699889709990720091730192740293750394761D-9, 
     &    -4.06476164514422552680590938629196667454705711D-11, 
     &    8.92169102045645255521798731675274885151428361D-13, 
     &    -1.993929586072107568723644347793789705630694749D-14, 
     &    4.51898002961991819165047655285559322839681901D-16, 
     &    -1.035651761218124701448341154221865666596091238D-17, 
     &    2.39521862102618674574028374300098038167894899D-19, 
     &    -5.58178587432500933628307450562541990556705462D-21, 
     &    1.309150755418321285812307399186592301749849833D-22, 
     &    -3.087419802426740293242279764866462431595565203D-24, 
     &    7.31597565270220342035790560925214859103339899D-26, 
     &    -1.740845657234000740989055147759702545340841422D-27, 
     &    4.15763564461389971961789962077522667348825413D-29, 
     &    -9.96214848828462210319400670245583884985485196D-31, 
     &    2.394034424896165300521167987893749562934279156D-32, 
     &    -5.76834735536739008429179316187765424407233225D-34, 
     &    1.393179479647007977827886603911548331732410612D-35, 
     &    -3.372121965485089470468473635254930958979742891D-37, 
     &    8.17820877756210262176477721487283426787618937D-39, 
     &    -1.987010831152385925564820669234786567541858996D-40, 
     &    4.83577851804055089628705937311537820769430091D-42 /

	xm = -cln(1 - z, -s)
	x2 = xm**2
	Li2series = xm - x2/4D0
	do j = 1, 25
	  xm = xm*x2
	  new = Li2series + xm*b(j)
	  if( new .eq. Li2series ) return
	  Li2series = new
	enddo
#ifdef WARNINGS
	print *, "Li2series: bad convergence"
#endif
	end

************************************************************************

	integer function eta(z1, s1, z2, s2, s12)
	implicit none
	double complex z1, z2
	double precision s1, s2, s12

	double precision im1, im2, im12

	im1 = DIMAG(z1)
	if( im1 .eq. 0 ) im1 = s1
	im2 = DIMAG(z2)
	if( im2 .eq. 0 ) im2 = s2
	im12 = DIMAG(z1*z2)
	if( im12 .eq. 0 ) im12 = s12

	if( im1 .lt. 0 .and. im2 .lt. 0 .and. im12 .gt. 0 ) then
	  eta = 1
	else
     &  if( im1 .gt. 0 .and. im2 .gt. 0 .and. im12 .lt. 0 ) then
	  eta = -1
	else
	  eta = 0
#ifdef WARNINGS
	  if( .not. (im2 .eq. 0 .and. DBLE(z2) .gt. 0 .or.
     &               im1 .eq. 0 .and. DBLE(z1) .gt. 0) .and.
     &        (im1 .eq. 0 .and. DBLE(z1) .lt. 0 .or.
     &         im2 .eq. 0 .and. DBLE(z2) .lt. 0 .or.
     &         im12 .eq. 0 .and. DBLE(z1*z2) .lt. 0) )
     &      print *, "eta not defined"
#endif
	endif
	end

