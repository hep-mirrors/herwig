* D0.F
* the scalar four-point function
* this file is part of LoopTools
* last modified 14 Dec 10 th

#include "defs.h"

#define legs 4
#define M(i) xpi(i)
#define P(i) xpi(i+legs)


	double complex function D0(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4

#include "lt.h"

	double complex res(0:1)
	double precision xpi(13)
	integer key, ier

	external D0softDR, D0collDR, D0soft, D0coll

	M(1) = m1
	M(2) = m2
	M(3) = m3
	M(4) = m4
	P(1) = p1 
	P(2) = p2
	P(3) = p3
	P(4) = p4
	P(5) = p1p2
	P(6) = p2p3
	P(7) = 0
	P(8) = 0
	P(9) = 0

	if( lambda .le. 0 ) then
	  call DDispatch(D0, xpi, D0softDR, D0collDR)
	  return
	endif

	ier = 0
	key = ibits(versionkey, KeyD0, 2)

	if( key .ne. 1 ) then
	  call ffxd0(res(0), xpi, ier)
	  if( ier .gt. warndigits ) then
	    ier = 0
	    call ffxd0r(res(0), xpi, ier)
	    if( ier .gt. warndigits ) key = ior(key, 2)
	    if( ier .ge. errdigits ) key = ior(key, 3)
	  endif
	endif

	if( key .ne. 0 ) then
	  call DDispatch(res(1), xpi, D0soft, D0coll)
	  if( key .gt. 1 .and.
     &        abs(res(0) - res(1)) .gt. maxdev*abs(res(0)) ) then
	    print *, "Discrepancy in D0:"
	    print *, "  p1   =", p1
	    print *, "  p2   =", p2
	    print *, "  p3   =", p3
	    print *, "  p4   =", p4
	    print *, "  p1p2 =", p1p2
	    print *, "  p2p3 =", p2p3
	    print *, "  m1   =", m1
	    print *, "  m2   =", m2
	    print *, "  m3   =", m3
	    print *, "  m4   =", m4
	    print *, "D0 a   =", res(0)
	    print *, "D0 b   =", res(1)
	    if( ier .gt. errdigits ) res(0) = res(1)
	  endif
	endif

	D0 = res(iand(key, 1))
	end

************************************************************************
* adapter code for C++

	subroutine d0sub(res, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double complex res
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4

	double complex D0
	external D0

	res = D0(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
	end

************************************************************************

	subroutine DDispatch(res, xpi, soft, coll)
	implicit none
	double complex res
	double precision xpi(13)
	external soft, coll

#include "lt.h"

	integer i, z, s, perm

	integer xpiperm(12), mperm(0:7)
	data xpiperm /
     &    O'1234561234', O'1635421243', O'5264131324',
     &    O'2341652341', O'6351242431', O'2546132314',
     &    O'3412563412', O'5462313142', O'3615243421',
     &    O'4123654123', O'4526314132', O'6153424213' /
	data mperm / O'1234561234', O'1234561234',
     &    O'5264131324', O'1234561234', O'4321561432',
     &    O'1635421243', O'5361421342', O'1234561234' /

* 0 1 1xxx  O'1234561234'
* 1 2 12xx  O'1234561234'
* 2 2 13xx  O'5264131324'
* 3 3 123x  O'1234561234'
* 4 2 14xx  O'4321561432'
* 5 3 124x  O'1635421243'
* 6 3 134x  O'5361421342'
* 7 4 xxxx  O'1234561234'

#define pj(p,j) ibits(p,3*(10-j),3)
#define mj(p,j) ibits(p,3*(4-j),3)

#define Px(j) P(pj(perm,j))
#define Mx(j) M(mj(perm,j))

	z = 0
	s = 0
	do i = 1, 12
	  perm = xpiperm(i)
	  if( abs(Mx(1)) .lt. eps ) then
	    if( abs(Px(1)) + abs(Mx(2)) .lt. eps ) then
	      call coll(res, xpi, perm)
	      if( res .ne. perm ) return
	    endif
	    if( s .eq. 0 .and.
     &        abs(Px(1) - Mx(2)) +
     &        abs(Px(4) - Mx(4)) .lt. acc ) s = perm
	    if( z .eq. 0 ) z = perm
	  endif
	enddo

	if( s .ne. 0 ) then
	  call soft(res, xpi, s)
	  return
	endif

	if( z .eq. 0 ) then
	  call D0m4(res, xpi)
	  return
	endif

	perm = z
	z = 0
	if( abs(Mx(2)) .lt. eps ) z = 1
	if( abs(Mx(3)) .lt. eps ) z = z + 2
	if( abs(Mx(4)) .lt. eps ) z = z + 4
	s = mperm(z)
	if( s .ne. O'1234561234' ) perm =
     &    pj(perm, pj(s, 1))*8**9 +
     &    pj(perm, pj(s, 2))*8**8 +
     &    pj(perm, pj(s, 3))*8**7 +
     &    pj(perm, pj(s, 4))*8**6 +
     &    pj(perm, pj(s, 5))*8**5 +
     &    pj(perm, pj(s, 6))*8**4 +
     &    mj(perm, mj(s, 1))*8**3 +
     &    mj(perm, mj(s, 2))*8**2 +
     &    mj(perm, mj(s, 3))*8**1 +
     &    mj(perm, mj(s, 4))*8**0

	goto (2, 2, 3, 2, 3, 3, 4) z

	call D0m3(res, xpi, perm)
	return

2	call D0m2(res, xpi, perm)
	return

3	call D0m1(res, xpi, perm)
	return

4	call D0m0(res, xpi)
	end

************************************************************************

	subroutine D0soft(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m3, p1, p2, p3, p4, p1p2, p2p3
	double precision r1, r3, r4
	double complex xs, x2, x3, y, c, fac
	double complex lxs, lx2, lx3, l1x2, l1x3, ly, lm

	double complex ln, spence, bdK
	external ln, spence, bdK

	m3 = Mx(3)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)
	p4 = Px(4)
	p1p2 = Px(5)
	p2p3 = Px(6)

	r1 = sqrt(p1)
	r4 = sqrt(p4)
	fac = .5D0/(r1*r4*(p1p2 - m3))
	xs = bdK(p2p3, r1, r4)
	lxs = -1
	if( xs .ne. 1 ) then
	  lxs = log(xs)
	  fac = 2*xs/((1 - xs)*(1 + xs))*fac
	endif

* massless case
	if( abs(m3) .lt. eps ) then
	  if( abs(p1 - p2) + abs(p3 - p4) .lt. acc ) then
	    res = -2*ln(-lambda/p1p2, 1D0)*lxs*fac
	    return
	  endif
	  y = (r1*(p3 - p4 + cI*eps))/(r4*(p2 - p1 + cI*eps))
	  ly = log(y)
	  c = ln(lambda/(r1*r4), 0D0) +
     &      ln((p2 - p1)/p1p2, p1 - p2) +
     &      ln((p3 - p4)/p1p2, p4 - p3)
	  if( xs .eq. 1 ) then
	    res = fac*(c - 2 - (1 + y)/(1 - y)*ly)
	  else
	    res = fac*(pi6 -
     &        spence(xs/y, 0D0) -
     &        (lxs + log(1/y))*log(1 - xs/y) -
     &        spence(xs*y, 0D0) -
     &        (lxs + ly)*(log(1 - xs*y) + .5D0*(lxs - ly)) +
     &        spence(xs**2, 0D0) +
     &        lxs*(2*log((1 - xs)*(1 + xs)) - c))
	  endif
	  return
	endif

* massive case
	r3 = sqrt(m3)
	x2 = bdK(p2, r1, r3)
	x3 = bdK(p3, r4, r3)
	lx2 = log(x2)
	lx3 = log(x3)
	l1x3 = log(1/x3)
	lm = 2*ln(r3*sqrt(lambda)/(m3 - p1p2), 1D0)
	if( xs .eq. 1 ) then
	  c = -2
	  if( abs(x2 - x3) .gt. acc ) then
	    c = (1 + x2/x3)/(1 - x2/x3)*(lx2 + l1x3) +
     &        (1 + x2*x3)/(1 - x2*x3)*(lx2 + lx3) + 2
	  else if( abs(x2 - 1) .gt. acc ) then
	    c = -2*(x2**2 + 1)/((x2 - 1)*(x2 + 1))*lx2
	  endif
	  res = fac*(lm - c)
	else
	  l1x2 = log(1/x2)
	  res = fac*( .5D0*pi**2 +
     &      lxs*(2*log((1 - xs)*(1 + xs)) - lm) +
     &      spence(xs**2, 0D0) + lx2**2 + lx3**2 -
     &      spence(xs/(x2*x3), 0D0) -
     &      (lxs + l1x2 + l1x3)*log(1 - xs/(x2*x3)) -
     &      spence(xs*x2/x3, 0D0) -
     &      (lxs + lx2 + l1x3)*log(1 - xs*x2/x3) -
     &      spence(xs/x2*x3, 0D0) -
     &      (lxs + l1x2 + lx3)*log(1 - xs/x2*x3) -
     &      spence(xs*x2*x3, 0D0) -
     &      (lxs + lx2 + lx3)*log(1 - xs*x2*x3) )
	endif
	end

************************************************************************

	double complex function bdK(x, m1, m2)
* this is actually -K from the Beenakker/Denner paper for D0soft
	implicit none
	double precision x, m1, m2

#include "lt.h"

	double precision d
	double complex t

	d = x - (m1 - m2)**2
	if( abs(d) .lt. acc ) then
	  bdK = 1
	else
	  t = 4*m1*m2/(d + cI*eps)
	  bdK = -t/(sqrt(1 - t) + 1)**2
	endif
	end

************************************************************************

	subroutine D0coll(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	logical ini
	data ini /.FALSE./

	Px(1) = max(minmass, 1D-14)
	res = perm

	if( ini ) return
	print *, "collinear-divergent D0, using mass cutoff ", Px(1)
	ini = .TRUE.
	end

************************************************************************
* IR-divergent D0 in dim reg
* from W. Beenakker and A. Denner, NPB 338 (1990) 349

	subroutine D0softDR(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m2, m3, m4, p2, p3, t, p2p3, q2, q3
	double precision r1, r3, r4, m24, sy
	double complex c, fac, xs, x2, x3, lxs, lx2, lx3, lm, y

	double complex bdK, ln, cln, lnrat, Li2omx2, Li2omx3
	external bdK, ln, cln, lnrat, Li2omx2, Li2omx3

c	PRINT *, "D0softDR"

	if( lambda .eq. -2 ) then
	  res = 0
	  return
	endif

	m2 = Mx(2)
	m3 = Mx(3)
	m4 = Mx(4)
	t = m3 - Px(5)
	p2p3 = Px(6)

	p2 = Px(2)
	q2 = m2 - p2
	p3 = Px(3)
	q3 = m4 - p3

	r1 = sqrt(m2)
	r4 = sqrt(m4)

	fac = .5D0/(r1*r4*t)
	xs = bdK(p2p3, r1, r4)
	lxs = -1
	if( xs .ne. 1 ) then
	  lxs = log(xs)
	  fac = 2*xs/((1 - xs)*(1 + xs))*fac
	endif

	if( abs(m3) .lt. eps ) then
	  if( abs(q2) + abs(q3) .lt. acc ) then
* qlbox14: D0(m2, m2, m4, m4; p1p2, p2p3; 0, m2, 0, m4)
c	PRINT *, "D0softDR: qlbox14"
	    res = 2*fac*lxs
	    if( lambda .ne. -1 ) res = res*lnrat(mudim, t)
	    return
	  endif

* qlbox15: D0(m2, p2, p3, m4; p1p2, p2p3; 0, m2, 0, m4)
* Beenakker-Denner Eq. (2.11)
c	PRINT *, "D0softDR: qlbox15"
	  if( lambda .eq. -1 ) then
	    res = fac*lxs
	    return
	  endif

	  if( abs(q2*q3) .lt. acc ) then
	    m24 = m2
	    if( abs(q2) .lt. acc ) m24 = m4
	    res = fac*( lxs*(lxs + log(mudim/m24) +
     &        2*lnrat(q2 + q3, t)) +
     &        Li2omx2(xs, 1D0, xs, 1D0) )
	    return
	  endif

	  y = r1*q3/(r4*q2)
	  sy = sign(.5D0, r1*q3) - sign(.5D0, r4*q2)

	  if( xs .eq. 1 ) then
	    res = fac*( -log(mudim/(r1*r4)) +
     &        lnrat(q2, t) + lnrat(q3, t) + 2 +
     &        (1 + y)/(1 - y)*ln(y, sy) )
	  else
	    res = fac*( -.5D0*ln(y, sy)**2 +
     &        lxs*(.5D0*lxs + lnrat(q2, t) + lnrat(q3, t) +
     &          log(mudim/(r1*r4))) +
     &        Li2omx2(xs, 1D0, xs, 1D0) -
     &        Li2omx2(xs, 1D0, y, sy) -
     &        Li2omx2(xs, 1D0, 1/y, -sy) )
	  endif
	  return
	endif

* qlbox16: D0(m2, p2, p3, m4; p1p2, p2p3; 0, m2, m3, m4)
* Beenakker-Denner Eq. (2.9)
c	PRINT *, "D0softDR: qlbox16"

	if( lambda .eq. -1 ) then
	  res = fac*lxs
	  return
	endif

	r3 = sqrt(m3)
	x2 = bdK(p2, r1, r3)
	x3 = bdK(p3, r4, r3)
	lx2 = log(x2)
	lx3 = log(x3)

	lm = 2*lnrat(sqrt(m3*mudim), t)

	if( xs .eq. 1 ) then
	  c = -2
	  if( abs(x2 - x3) .gt. acc ) then
	    c = (1 + x2/x3)/(1 - x2/x3)*(lx2 + log(1/x3)) +
     &        (1 + x2*x3)/(1 - x2*x3)*(lx2 + lx3) + 2
	  else if( abs(x2 - 1) .gt. acc ) then
	    c = -2*(x2**2 + 1)/((x2 - 1)*(x2 + 1))*lx2
	  endif
	  res = fac*(c - lm)
	else
	  res = fac*(lm*lxs - lx2**2 - lx3**2 +
     &      Li2omx2(xs, 1D0, xs, 1D0) -
     &      Li2omx3(xs, 1D0, x2, 1D0, x3, 1D0) -
     &      Li2omx3(xs, 1D0, 1/x2, -1D0, 1/x3, -1D0) -
     &      Li2omx3(xs, 1D0, x2, 1D0, 1/x3, -1D0) -
     &      Li2omx3(xs, 1D0, 1/x2, -1D0, x3, 1D0))
	endif
	end

************************************************************************

	subroutine D0collDR(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	integer z, s

* see D0table.f
	integer pperm(0:127)
	data pperm /
     &    -898440548,  -898440548,  -865381277,  -898440548,
     &    -898440548, -1972182372,  -865381277, -1972182372,
     &    -898440548,  -898440548,  -865381277, -1972182372,
     &    -898440548, -1972182372,  -865381277,  1249043100,
     &    -898440548,  -898440548, -1939123101, -1939123101,
     &    -898440548, -1972182372, -1939123101,  1368319628,
     &    -898440548,  -898440548, -1939123101,  1282102371,
     &    -898440548, -1972182372, -1939123101,   175301276,
     &    -898440548, -1921432484,  -865381277, -1921432484,
     &    -898440548, -1972182372,  -865381277, -1972182372,
     &    -898440548, -1921432484,  -865381277,  1299792988,
     &    -898440548, -1921432484,  -865381277,  1249043100,
     &    -898440548, -1921432484, -1939123101, -1939123101,
     &    -898440548, -1972182372, -1939123101,  1368319628,
     &    -898440548, -1921432484, -1939123101,  1282102371,
     &    -898440548, -1972182372, -1939123101,   175301276,
     &    -898440548,  -898440548, -1904860509, -1904860509,
     &    -898440548, -1972182372, -1904860509, -1972182372,
     &    -898440548,  -898440548, -1904860509,  1316364963,
     &    -898440548, -1972182372, -1904860509,  1249043100,
     &    -898440548,  -898440548, -1904860509, -1939123101,
     &    -898440548, -1972182372, -1904860509,  1368319628,
     &    -898440548,  -898440548, -1904860509,  1282102371,
     &    -898440548, -1972182372, -1904860509,   175301276,
     &    -898440548, -1921432484, -1904860509,  1907239051,
     &    -898440548, -1972182372, -1904860509,  1907239051,
     &    -898440548, -1921432484, -1904860509,   833497227,
     &    -898440548, -1972182372, -1904860509,   833497227,
     &    -898440548, -1921432484, -1904860509,  1907239051,
     &    -898440548, -1972182372, -1904860509,   850499737,
     &    -898440548, -1921432484, -1904860509,   833497227,
     &    -898440548, -1972182372, -1904860509,   175301276 /

	z = 0
	if( abs(Mx(3)) .lt. eps ) z = 1
	if( abs(Mx(4)) .lt. eps ) z = z + 2
	if( abs(Px(2)) .lt. eps ) z = z + 4
	if( abs(Px(3)) .lt. eps ) z = z + 8
	if( abs(Px(4)) .lt. eps ) z = z + 16
	if( abs(Px(5)) .lt. eps ) z = z + 32
	if( abs(Px(6)) .lt. eps ) z = z + 64
	s = pperm(z)
	if( iand(s, O'7777777777') .ne. O'1234561234' ) perm =
     &    pj(perm, pj(s, 1))*8**9 +
     &    pj(perm, pj(s, 2))*8**8 +
     &    pj(perm, pj(s, 3))*8**7 +
     &    pj(perm, pj(s, 4))*8**6 +
     &    pj(perm, pj(s, 5))*8**5 +
     &    pj(perm, pj(s, 6))*8**4 +
     &    mj(perm, mj(s, 1))*8**3 +
     &    mj(perm, mj(s, 2))*8**2 +
     &    mj(perm, mj(s, 3))*8**1 +
     &    mj(perm, mj(s, 4))*8**0

	goto (22,22,22,23, 22,22,22,23, 10,11,12,13)
     &    ibits(s, 30, 2) + ibits(z, 0, 2)*4 - 3

	call D0m2p3(res, xpi, perm)
	return

23	call D0m1p3(res, xpi, perm)
	return

22	call D0m1p2(res, xpi, perm)
	return

13	call D0m0p3(res, xpi, perm)
	return

12	call D0m0p2(res, xpi, perm)
	return

11	call D0m0p1(res, xpi, perm)
	return

10	call D0m0p0(res, xpi, perm)
	end

************************************************************************
* qlbox1: D0(0, 0, 0, 0; p1p2, p2p3; 0, 0, 0, 0)
* Bern, Dixon, Kosower, NPB 412 (1994) 751 [hep-ph/9306240], Eq. (I.11)

	subroutine D0m0p0(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, fac

	double complex lnrat
	external lnrat

c	PRINT *, "D0m0p0: qlbox1"
	s = -Px(5)
	t = -Px(6)
	fac = 1/(s*t)
	if( lambda .eq. -2 )then
	  res = 4*fac
	else if( lambda .eq. -1 ) then
	  res = 2*fac*(-lnrat(t, mudim) - lnrat(s, mudim))
	else
	  res = fac*(lnrat(t, mudim)**2 + lnrat(s, mudim)**2 -
     &      lnrat(t, s)**2 - pi**2)
	endif
	end

************************************************************************
* qlbox2: D0(0, 0, 0, p4; p1p2, p2p3; 0, 0, 0, 0)
* One-mass integral as given in
* Ellis, Giele, Zanderighi, Eq. (A22).

	subroutine D0m0p1(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, m4, fac
	double complex l1, l2

	double complex lnrat, Li2omrat
	external lnrat, Li2omrat

c	PRINT *, "D0m0p1: qlbox2"
	s = -Px(5)
	t = -Px(6)
	m4 = -Px(4)
	fac = 1/(s*t)
	if( lambda .eq. -2 ) then
	  res = 2*fac
	else if( lambda .eq. -1 ) then
	  res = 2*fac*(lnrat(m4, mudim) -
     &      lnrat(t, mudim) - lnrat(s, mudim))
	else
	  l1 = sqrt(lnrat(t, mudim)**2 + lnrat(m4, t)**2 +
     &              lnrat(s, mudim)**2 + lnrat(m4, s)**2)
          l2 = sqrt(lnrat(m4, mudim)**2 + lnrat(t, s)**2)
	  res = fac*((l1 - l2)*(l1 + l2) +
     &      2*(Li2omrat(t, m4) + Li2omrat(s, m4) - pi6))
	endif
	end

************************************************************************

	subroutine D0m0p2(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, q2, q3, q4, fac, r
	double complex ls, lt, lq2, lq3, lq4

	double complex lnrat, lndiv0, lndiv1, Li2omrat, Li2omrat2
	external lnrat, lndiv0, lndiv1, Li2omrat, Li2omrat2

c	PRINT *, "D0m0p2"
	s = -Px(5)
	t = -Px(6)
	fac = 1/(s*t)
	q4 = -Px(4)
	q3 = -Px(3)

	if( abs(q3) .lt. eps ) then
* qlbox3: D0(0, p2, 0, p4; p1p2, p2p3; 0, 0, 0, 0)
* Bern, Dixon, Kosower, NPB 412 (1994) 751 [hep-ph/9306240], Eq. (I.13)
c	PRINT *, "D0m0p2: qlbox3"

	  if( lambda .eq. -2 ) then
	    res = 0
	    return
	  endif

	  q2 = -Px(2)
	  r = 1 - q2*q4*fac

* Use expansion only in cases where signs (s,t,m2,m4) are not
* ++-- or --++
	  if( abs(r) .lt. 1D-6 .and.
     &      (fac .lt. 0 .or. q2*q4 .lt. 0) ) then
* expanded case
	    if( lambda .eq. -1 ) then
	      res = -(2 + r)*fac
	    else
	      res = fac*(2 - .5D0*r +
     &          (2 + r)*(lnrat(s, mudim) + lnrat(t, q4)) +
     &          2*(lndiv0(q4, t) + lndiv0(q4, s)) +
     &          r*(lndiv1(q4, t) + lndiv1(q4, s)))
	    endif
	  else
* general case
	    fac = 1/(s*t - q2*q4)
	    if( lambda .eq. -1 ) then
	      res = 2*fac*(lnrat(q2, s) + lnrat(q4, t))
	    else
	      ls = lnrat(s, mudim)
	      lt = lnrat(t, mudim)
	      lq2 = lnrat(q2, mudim)
	      lq4 = lnrat(q4, mudim)
	      res = fac*(
     &          (ls - lq2)*(ls + lq2) +
     &          (lt - lq4)*(lt + lq4) - lnrat(s, t)**2 +
     &          2*(Li2omrat2(q2, s, q4, t) -
     &            Li2omrat(q2, s) - Li2omrat(q2, t) -
     &            Li2omrat(q4, s) - Li2omrat(q4, t)) )
	    endif
	  endif
	  return
	endif

* qlbox4: D0(0, 0, p3, p4; p1p2, p2p3; 0, 0, 0, 0)
* Bern, Dixon, Kosower, NPB 412 (1994) 751 [hep-ph/9306240], Eq. (I.14)
c	PRINT *, "D0m0p2: qlbox4"

	if( lambda .eq. -2 ) then
	  res = fac
	else if( lambda .eq. -1 ) then
	  res = -fac*(lnrat(s, q3) + lnrat(t, q4) + lnrat(t, mudim))
	else
	  ls = lnrat(s, mudim)
	  lt = lnrat(t, mudim)
	  lq3 = lnrat(q3, mudim)
	  lq4 = lnrat(q4, mudim)
	  res = fac*(
     &      .5D0*((ls - lq3)*(ls + lq3) +
     &            (lt - lq4)*(lt + lq4) + lt**2) +
     &      lnrat(s, q3)*lnrat(s, q4) - lnrat(s, t)**2 -
     &      2*(Li2omrat(q3, t) + Li2omrat(q4, t)) )
	endif
	end

************************************************************************
* qlbox5: D0(0, p2, p3, p4; p1p2, p2p3; 0, 0, 0, 0)
* Bern, Dixon, Kosower, NPB 412 (1994) 751 [hep-ph/9306240], Eq. (I.15)
* or from hep-ph/0508308v3 Eq. (A27)
* (v3 corrects previous versions)

	subroutine D0m0p3(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, q2, q3, q4, fac, r

	double complex lnrat, lndiv0, lndiv1, Li2omrat, Li2omrat2
	external lnrat, lndiv0, lndiv1, Li2omrat, Li2omrat2

c	PRINT *, "D0m0p3: qlbox5"

	if( lambda .eq. -2 ) then
	  res = 0
	  return
	endif

	s = -Px(5)
	t = -Px(6)
	fac = 1/(s*t)
	q2 = -Px(2)
	q3 = -Px(3)
	q4 = -Px(4)

	r = 1 - q2*q4*fac

* Use expansion only in cases where signs of (s,t,q2,q4) are
* not ++-- or --++
	if( abs(r) .lt. 1D-6 .and.
     &    (fac .lt. 0 .or. q2*q4 .lt. 0) ) then
* expanded case
	  if( lambda .eq. -1 ) then
	    res = -.5D0*(2 + r)*fac
	  else
	    res = lndiv0(q4, t)
	    res = fac*(
     &        .5D0*(2 + r)*(2 + (1 + q4/t)*res -
     &          lnrat(mudim, s) - lnrat(q3, t)) +
     &        r*(lndiv1(q4, t) - res - 1) )
	  endif
	else
* general case
	  fac = 1/(s*t - q2*q4)
	  if( lambda .eq. -1 ) then
	    res = fac*(lnrat(q2, t) + lnrat(q4, s))
	  else
	    res = fac*(
     &        (lnrat(q3, t) + lnrat(mudim, t))*lnrat(q2, t) +
     &        (lnrat(q3, s) + lnrat(mudim, s))*lnrat(q4, s) -
     &        .5D0*(lnrat(t, q2)**2 + lnrat(s, q4)**2) -
     &        lnrat(s, t)**2 -
     &        2*(Li2omrat(q2, s) + Li2omrat(q4, t) -
     &          Li2omrat2(q2, s, q4, t)) )
	  endif
	endif
	end

************************************************************************

	subroutine D0m1p2(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m4, s, t, q3, q4, fac
	double complex lm, ls, lt, lq
	integer ir

	double complex lnrat, Li2omrat, Li2omrat2
	external lnrat, Li2omrat, Li2omrat2

c	PRINT *, "D0m1p2"
	m4 = Mx(4)
	s = -Px(5)
	t = m4 - Px(6)
	fac = 1/(s*t)
	q3 = m4 - Px(3)
	q4 = m4 - Px(4)

	ir = 0
	if( abs(q3) .lt. acc ) ir = 1
	if( abs(q4) .lt. acc ) then
	  ir = ir + 1
	  q4 = q3
	endif

	if( lambda .eq. -2 ) then
	  res = .5D0*(2 + ir)*fac
	  return
	endif

	goto (1, 2) ir

* qlbox8: D0(0, 0, p3, p4; p1p2, p2p3; 0, 0, 0, m4)
c	PRINT *, "D0m1p2: qlbox8"
	lm = lnrat(s, mudim)
	if( lambda .eq. -1 ) then
	  res = fac*(lnrat(q3, t) + lnrat(q4, t) - lm)
	else
	  ls = lnrat(s, m4)
	  res = fac*(-2*(Li2omrat(q3, t) + Li2omrat(q4, t)) -
     &      Li2omrat2(q3, s, q4, m4) - pi6 +
     &      .5D0*(lm - ls)*(lm + ls) + 2*lm*lnrat(t, m4) -
     &      lnrat(q3, mudim)*lnrat(q3, m4) -
     &      lnrat(q4, mudim)*lnrat(q4, m4))
	endif
	return

1	continue
* qlbox7: D0(0, 0, m4, p4; p1p2, p2p3; 0, 0, 0, m4)
c	PRINT *, "D0m1p2: qlbox7"
	ls = lnrat(s, m4)
	lt = lnrat(t, m4)
	lm = lnrat(mudim, m4)
	lq = lnrat(q4, m4)
	if( lambda .eq. -1 ) then
	  res = fac*(1.5D0*lm - 2*lt - ls + lq)
	else
	  res = fac*(2*ls*lt - lq**2 - 5*pi12 +
     &      lm*(.75D0*lm - 2*lt - ls + lq) -
     &      2*Li2omrat(q4, t))
	endif
	return

2	continue
* qlbox6: D0(0, 0, m4, m4; p1p2, p2p3; 0, 0, 0, m4)
c	PRINT *, "D0m1p2: qlbox6"
	ls = lnrat(s, m4)
	lt = lnrat(t, m4)
	lm = lnrat(mudim, m4)
	if( lambda .eq. -1 ) then
	  res = fac*(2*(lm - lt) - ls)
	else
	  res = fac*((lm - ls)*(lm - 2*lt) - .5D0*pi**2)
	endif
	end

************************************************************************

	subroutine D0m1p3(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, m4, q2, q3, q4, m4mu, fac

	double complex lnrat, Li2omrat, Li2omrat2
	external lnrat, Li2omrat, Li2omrat2

c	PRINT *, "D0m1p3"

	if( lambda .eq. -2 ) then
	  res = 0
	  return
	endif

	q2 = -Px(2)
	s = -Px(5)
	m4 = Mx(4)
	q3 = m4 - Px(3)
	q4 = m4 - Px(4)
	t = m4 - Px(6)

	if( abs(t) .lt. acc ) then
	  t = q4
	  q4 = 0
	  s = q2
	  q2 = -Px(5)
	endif

	m4mu = sqrt(m4*mudim)

* qlbox9: D0(0, p2, p3, m4; p1p2, p2p3; 0, 0, 0, m4)
	if( abs(q4) .lt. acc ) then
c	PRINT *, "D0m1p3: qlbox9"
	  fac = 1/(s*t)
	  if( lambda .eq. -1 ) then
	    res = -fac*(lnrat(t, m4mu) + lnrat(s, q2))
	  else
	    res = fac*(Li2omrat2(q3, q2, t, m4) + 2*Li2omrat(s, q2) +
     &        lnrat(t, m4mu) + lnrat(s, q2) + pi12)
	  endif
	  return
	endif

* qlbox10: D0(0, p2, p3, p4; p1p2, p2p3; 0, 0, 0, m4)
c	PRINT *, "D0m1p3: qlbox10"
	fac = 1/(s*t - q2*q4)
	res = fac*(lnrat(q2, mudim) + lnrat(q4, mudim) -
     &    lnrat(s, mudim) - lnrat(t, mudim))
	if( lambda .ne. -1 ) then
	  res = 2*res*lnrat(m4mu, t) +
     &      fac*(Li2omrat2(q3, q2, t, m4) - Li2omrat2(q3, s, q4, m4) +
     &        2*(Li2omrat2(q2, s, q4, t) +
     &          Li2omrat(q2, s) - Li2omrat(t, q4)))
	endif
	end

************************************************************************

	subroutine D0m2p3(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision s, t, m3, m4, q3, q4, p3, fac, m3mu, m4mu
	double precision p34, c, s3t, s4s
	double complex ls, lt, lq3, lq4, d
	double complex x43(4), r3t, r4s, r43p, r43m
	double complex logs, dilogs
	integer ir, case

	double complex lnrat, cln, Li2rat, Li2omrat, Li2omrat2
	external lnrat, cln, Li2rat, Li2omrat, Li2omrat2

c	PRINT *, "D0m2p3"
	m3 = Mx(3)
	s = m3 - Px(5)
	q3 = m3 - Px(2)
	m4 = Mx(4)
	t = m4 - Px(6)
	q4 = m4 - Px(4)

	if( abs(s*t) .lt. eps ) then
	  s = q3
	  q3 = m3 - Px(5)
	  t = q4
	  q4 = m4 - Px(6)
	endif

	fac = 1/(s*t - q3*q4)

	ir = 0
	if( abs(q3) .lt. acc ) ir = 1
	if( abs(q4) .lt. acc ) then
	  ir = ir + 1
	  q4 = q3
	  s = t
	  t = m3 - Px(5)
	  m4 = m3
	  m3 = Mx(4)
	endif

	if( lambda .eq. -2 ) then
	  res = .5D0*fac*ir
	  return
	endif

	if( lambda .eq. -1 ) goto (10, 11, 12) ir + 1

	p3 = Px(3)
	if( abs(p3) .lt. eps ) then
	  case = 1
	  logs = lnrat(m3, m4)**2
	else
	  p34 = p3 + m3 - m4
	  c = -4*p3*m3
	  d = sqrt(DCMPLX(p34**2 + c))

	  x43(1) = -p34 - d
	  x43(2) =  p34 - d
	  if( abs(x43(1)) .lt. abs(x43(2)) ) then
	    x43(1) = c/x43(2)
	  else
	    x43(2) = c/x43(1)
	  endif

	  p34 = -p3 + m3 - m4
	  c = -4*p3*m4

	  x43(3) = -p34 - d
	  x43(4) =  p34 - d
	  if( abs(x43(3)) .lt. abs(x43(4)) ) then
	    x43(3) = c/x43(4)
	  else
	    x43(4) = c/x43(3)
	  endif

	  if( abs(DIMAG(d)) .lt. eps ) then
	    case = 2
	    logs = lnrat(x43(1), x43(3))**2 +
     &             lnrat(x43(2), x43(4))**2
	  else
	    case = 3
	    r43p = x43(1)/x43(3)
	    r43m = x43(2)/x43(4)
	    logs = cln(r43p, 0D0)**2 + cln(r43m, 0D0)**2
	  endif
	endif

	goto (1, 2) ir

* qlbox13: D0(0, p2, p3, p4; p1p2, p2p3; 0, 0, m3, m4)
c	PRINT *, "D0m2p3: qlbox13"

	ls = lnrat(s, mudim)
	lt = lnrat(t, mudim)
	lq3 = lnrat(q3, mudim)
	lq4 = lnrat(q4, mudim)

	if( case .eq. 1 ) then
	  dilogs = Li2omrat2(q3, t, -1D0, -1D0) +
     &             Li2omrat2(q3, t, m4, m3) +
     &             Li2omrat2(q4, s, m3, m4) +
     &             Li2omrat2(q4, s, -1D0, -1D0)
	else if( case .eq. 2 ) then
	  dilogs = Li2omrat2(q3, t, x43(4), x43(2)) +
     &             Li2omrat2(q3, t, x43(3), x43(1)) +
     &             Li2omrat2(q4, s, x43(1), x43(3)) +
     &             Li2omrat2(q4, s, x43(2), x43(4))
	else
	  r3t = q3/t
	  s3t = sign(.5D0, q3) - sign(.5D0, t)
	  r4s = q4/s
	  s4s = sign(.5D0, q4) - sign(.5D0, s)
	  dilogs = Li2rat(r3t,s3t, 1/r43m,0D0) +
     &             Li2rat(r3t,s3t, 1/r43p,0D0) +
     &             Li2rat(r4s,s4s, r43p,0D0) +
     &             Li2rat(r4s,s4s, r43m,0D0)
	endif

	res = -fac*(dilogs + .5D0*logs + lq3**2 + lq4**2 +
     &    2*(Li2omrat(q3, s) + Li2omrat(q4, t) -
     &      Li2omrat2(q3, s, q4, t) - ls*lt) +
     &    (lt - lq3)*log(m3/mudim) + (ls - lq4)*log(m4/mudim))
	return

10	res = fac*(lnrat(q3, mudim) + lnrat(q4, mudim) -
     &    lnrat(s, mudim) - lnrat(t, mudim))
	return

1	continue
* qlbox12: D0(0, m3, p3, p4; p1p2, p2p3; 0, 0, m3, m4)
c	PRINT *, "D0m2p3: qlbox12"

	m3mu = sqrt(m3*mudim)
	ls = lnrat(s, m3mu)
	lt = lnrat(t, m3mu)
	lq4 = lnrat(q4, m3mu)

	if( case .eq. 1 ) then
	  dilogs = 0
	else if( case .eq. 2 ) then
	  dilogs = Li2omrat2(q4, s, x43(1), x43(3)) +
     &             Li2omrat2(q4, s, x43(2), x43(4))
	else
	  r4s = q4/s
	  s4s = sign(.5D0, q4) - sign(.5D0, s)
	  dilogs = Li2rat(r4s,s4s, r43p,0D0) +
     &             Li2rat(r4s,s4s, r43m,0D0)
	endif

	res = -fac*(dilogs + .5D0*logs + pi12 +
     &    2*(Li2omrat(q4, t) - ls*lt) +
     &    lq4**2 + (ls - lq4)*log(m4/m3))
	return

11	m3mu = sqrt(m3*mudim)
	res = fac*(lnrat(q4, m3mu) - lnrat(s, m3mu) - lnrat(t, m3mu))
	return

2	continue
* qlbox11: D0(0, m3, p3, m4; p1p2, p2p3; 0, 0, m3, m4)
c	PRINT *, "D0m2p3: qlbox11"

	m3mu = sqrt(m3*mudim)
	m4mu = sqrt(m4*mudim)
	ls = lnrat(s, m3mu)
	lt = lnrat(t, m4mu)

	res = fac*(.25D0*log(m3/m4)**2 - .5D0*(logs + pi**2) +
     &    2*lnrat(s, sqrt(m3*mudim))*lnrat(t, sqrt(m4*mudim)))
	return

12	res = -fac*(lnrat(s, sqrt(m3*mudim)) +
     &    lnrat(t, sqrt(m4*mudim)))
	end

************************************************************************
* this routine is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* it is used for double-checking the results of FF
* M. Rauch: implemented the log branch cuts for k13 < 2
* (from Denner, Nierste, Scharf; Nucl Phys B367 (1991) 637)

#define AddEps(k) k*DCMPLX(1D0, -sign(eps, k))
#define k2r(k) (.5D0*k*(1 + sqrt(DCMPLX((1 - 2/k)*(1 + 2/k)))))

	subroutine D0m4(res, xpi)
	implicit none
	double complex res
	double precision xpi(13)

#include "lt.h"

	double precision tmp, ir1324, gamma, s1, s2
	double precision kij(6), irij(6), ix(2,4)
	double complex rij(6), x(2,4), l(2,4), q13, q24
	double complex a, b, c, d, disc, ki, etas
	integer j

	double precision k12, k13, k14, k23, k24, k34
	double precision ir12, ir13, ir14, ir23, ir24, ir34
	double complex r12, r14, r13, r23, r24, r34
	equivalence (kij(1), k12), (rij(1), r12), (irij(1), ir12)
	equivalence (kij(2), k23), (rij(2), r23), (irij(2), ir23)
	equivalence (kij(3), k34), (rij(3), r34), (irij(3), ir34)
	equivalence (kij(4), k14), (rij(4), r14), (irij(4), ir14)
	equivalence (kij(5), k13), (rij(5), r13), (irij(5), ir13)
	equivalence (kij(6), k24), (rij(6), r24), (irij(6), ir24)

	double complex cln, xspence, xeta, xetatilde
	integer eta
	external cln, xspence, xeta, xetatilde, eta

	k12 = (M(1) + M(2) - P(1))/sqrt(M(1)*M(2))
	k23 = (M(2) + M(3) - P(2))/sqrt(M(2)*M(3))
	k34 = (M(3) + M(4) - P(3))/sqrt(M(3)*M(4))
	k14 = (M(1) + M(4) - P(4))/sqrt(M(1)*M(4))
	k13 = (M(1) + M(3) - P(5))/sqrt(M(1)*M(3))
	k24 = (M(2) + M(4) - P(6))/sqrt(M(2)*M(4))

* test if r_13 can be made real by a permutation
* if one of the r_ij is real r_13 must be made real => case 1
	if( abs(k13) .ge. 2 ) then
*	  nothing to do
* otherwise try all permutations
	else if( abs(k12) .ge. 2 ) then
* 2 <-> 3
	  tmp = k12
	  k12 = k13
	  k13 = tmp
	  tmp = k24
	  k24 = k34
	  k34 = tmp
	else if( abs(k14) .ge. 2 ) then
* 3 <-> 4
	  tmp = k13
	  k13 = k14
	  k14 = tmp
	  tmp = k23
	  k23 = k24
	  k24 = tmp
	else if( abs(k23) .ge. 2 ) then
* 1 <-> 2
	  tmp = k13
	  k13 = k23
	  k23 = tmp
	  tmp = k14
	  k14 = k24
	  k24 = tmp
	else if( abs(k24) .ge. 2 ) then
*  1 -> 4, 2 -> 1, 3 -> 2, 4 -> 3
	  tmp = k12
	  k12 = k23
          k23 = k34
	  k34 = k14
	  k14 = tmp
	  tmp = k13
	  k13 = k24
	  k24 = tmp
	else if( abs(k34) .ge. 2 ) then
* 1 <-> 4
	  tmp = k12
	  k12 = k24
	  k24 = tmp
	  tmp = k13
	  k13 = k34
	  k34 = tmp
* 	else
* nothing found => all r_ij on the complex unit circle => case 2
	endif

	r12 = k2r(k12)
	r23 = k2r(k23)
	r34 = k2r(k34)
	r14 = k2r(k14)
	r13 = 1/k2r(k13)
	r24 = 1/k2r(k24)

	do j = 1, 6
	  if( DIMAG(rij(j)) .eq. 0 ) then
	    ki = kij(j) - cI*eps
	    irij(j) = sign(1D0, abs(rij(j)) - 1)*
     &        DIMAG(k2r(ki))
	  else
	    irij(j) = 0
	  endif
	enddo

	ir1324 = sign(1D0, DBLE(r24))*ir13 -
     &           sign(1D0, DBLE(r13))*ir24

	a = k34/r24 - k23 + (k12 - k14/r24)*r13
	b = (1/r13 - r13)*(1/r24 - r24) + k12*k34 - k14*k23
	c = k34*r24 - k23 + (k12 - k14*r24)/r13
	d = k23 + (r24*k14 - k12)*r13 - r24*k34
	disc = sqrt(b**2 - 4*a*(c + cI*eps*d))
	ix(1,4) = DIMAG(.5D0/a*(b - disc))
	ix(2,4) = DIMAG(.5D0/a*(b + disc))

	disc = sqrt(b**2 - 4*a*c)
	x(1,4) = .5D0/a*(b - disc)
	x(2,4) = .5D0/a*(b + disc)
	if( abs(x(1,4)) .gt. abs(x(2,4)) ) then
	  x(2,4) = c/(a*x(1,4))
	else    
	  x(1,4) = c/(a*x(2,4))
	endif

	x(1,1) = x(1,4)/r24
	x(2,1) = x(2,4)/r24
	x(1,2) = x(1,4)*r13/r24
	x(2,2) = x(2,4)*r13/r24
	x(1,3) = x(1,4)*r13
	x(2,3) = x(2,4)*r13

	s1 = sign(1D0, DBLE(x(1,4)))
	s2 = sign(1D0, DBLE(x(2,4)))
	ix(1,1) = ix(1,4)*DBLE(x(1,1))*s1
	ix(2,1) = ix(2,4)*DBLE(x(2,1))*s2
	ix(1,2) = ix(1,4)*DBLE(x(1,2))*s1
	ix(2,2) = ix(2,4)*DBLE(x(2,2))*s2
	ix(1,3) = ix(1,4)*DBLE(x(1,3))*s1
	ix(2,3) = ix(2,4)*DBLE(x(2,3))*s2
 
	res = 0
	do j = 1, 4
	  res = res + Sgn(j)*(
     &      xspence(x(1,j), ix(1,j), rij(j), irij(j)) +
     &      xspence(x(1,j), ix(1,j), 1/rij(j), -irij(j)) )
	enddo

	gamma = sign(1D0, DBLE(a*(x(2,4) - x(1,4))))
	l(1,4) = c2ipi*eta(r13, ir13, 1/r24, -ir24, ir1324)
	l(2,4) = l(1,4)

	if( DIMAG(r13) .eq. 0 ) then
	  r12 = k12 - r24*k14
	  r23 = k23 - r24*k34
	  r34 = k34 - r13*k14
	  r14 = k23 - r13*k12
	  q13 = k13 - 2*r13
	  q24 = k24 - 2*r24

	  c = gamma*sign(1D0, DIMAG(r24) + ir24)
	  l(1,1) = cln(-x(1,1), -ix(1,1)) +
     &      cln(r14 - q13/x(1,1), -1D0) +
     &      cln((r12 - q24*x(1,4))/d, c)
	  l(2,1) = cln(-x(2,1), -ix(2,1)) +
     &      cln(r14 - q13/x(2,1), -1D0) +
     &      cln((r12 - q24*x(2,4))/d, -c)

	  c = gamma*sign(1D0, DBLE(r13)*(DIMAG(r24) + ir24))
	  l(1,2) = cln(-x(1,2), -ix(1,2)) +
     &      cln(r14 - q13/x(1,1), -1D0) +
     &      cln((r23 - q24*x(1,3))/d, c)
	  l(2,2) = cln(-x(2,2), -ix(2,2)) +
     &      cln(r14 - q13/x(2,1), -1D0) +
     &      cln((r23 - q24*x(2,3))/d, -c)

	  l(1,3) = cln(-x(1,3), -ix(1,3)) +
     &      cln(r34 - q13/x(1,4), -1D0) +
     &      cln((r23 - q24*x(1,3))/d, c)
	  l(2,3) = cln(-x(2,3), -ix(2,3)) +
     &      cln(r34 - q13/x(2,4), -1D0) +
     &      cln((r23 - q24*x(2,3))/d, -c)

	  etas =
     &      xetatilde(x(1,4), ix(1,4), r13, ir13, l(1,3)) +
     &      xetatilde(x(1,4), ix(1,4), 1/r24, -ir24, l(1,1)) -
     &      xetatilde(x(1,4), ix(1,4), r13/r24, ir1324, l(1,2)) +
     &      xetatilde(x(1,4), ix(1,4), -r13/r24, -ir1324, l(1,4))
	else
	  do j = 1, 3
	    l(1,j) = log(-x(1,j)) +
     &        cln(kij(j) - 1/x(1,j) - x(1,j), -x(1,j)*b*gamma)
	    l(2,j) = log(-x(2,j)) +
     &        cln(kij(j) - 1/x(2,j) - x(2,j), -x(2,j)*b*gamma)
	  enddo

	  etas =
     &      xeta(x(1,4), ix(1,4), r13, ir13, ix(1,3), l(1,3)) +
     &      xeta(x(1,4), ix(1,4), 1/r24, -ir24, ix(1,1), l(1,1)) -
     &      xeta(x(1,4), ix(1,4), r13/r24, ir1324, ix(1,2), l(1,2)) +
     &      xeta(x(1,4), ix(1,4), -r13/r24, -ir1324, ix(1,4), l(1,4))*
     &        (1 - sign(1D0, DBLE(b))*gamma)
 	endif

	res = (res - c2ipi*etas + (l(2,2) - l(1,2))*l(1,4))/
     &    (sqrt(M(1)*M(2)*M(3)*M(4))*disc)
	end

************************************************************************

	subroutine D0m3(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m2, m3, m4, p1, p2, p3, p4, p1p2, p2p3
	double precision m, k12, k13, k14, k23, k24, k34
	double precision ir12, ir14, ir24, ix1(2), ix4(2)
	double complex r12, r14, r24, q12, q24
	double complex x1(2), x4(2), l4(2)
	double complex a, b, c, d

	double complex cln, xspence, xetatilde
	external cln, xspence, xetatilde

	m2 = Mx(2)
	m3 = Mx(3)
	m4 = Mx(4)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)
	p4 = Px(4)
	p1p2 = Px(5)
	p2p3 = Px(6)

	m = sqrt(m3*m4)
	k23 = (m4 - p4)/m
	k12 = (m4 + m3 - p3)/m
	r12 = k2r(k12)
	ir12 = 0
	if( k12 .lt. -2 ) ir12 = sign(10D0, 1 - abs(r12))

	m = sqrt(m2*m3)
	k34 = (m2 - p1)/m
	k14 = (m2 + m3 - p2)/m
	r14 = k2r(k14)
	ir14 = 0
	if( k14 .lt. -2 ) ir14 = sign(10D0, 1 - abs(r14))

	k13 = (m3 - p1p2)/m3

	m = sqrt(m2*m4)
	k24 = (m2 + m4 - p2p3)/m
	r24 = k2r(k24)
	ir24 = 0
	if( k24 .lt. -2 ) ir24 = sign(10D0, 1 - abs(r24))

	q24 = r24 - 1/r24
	q12 = k12 - r24*k14

	a = k34/r24 - k23
	b = k12*k34 - k13*q24 - k14*k23
	c = k13*q12 + r24*k34 - k23
	d = sqrt(DCMPLX((k12*k34 - k13*k24 - k14*k23)**2 -
     &    4*(k13*(k13 - k23*(k12 - k14*k24)) +
     &       k23*(k23 - k24*k34) + k34*(k34 - k13*k14))))
	x4(1) = .5D0/a*(b - d)
	x4(2) = .5D0/a*(b + d)
	if( abs(x4(1)) .gt. abs(x4(2)) ) then
	  x4(2) = c/(a*x4(1))
	else
	  x4(1) = c/(a*x4(2))
	endif

	d = -k34*r24 + k23
	ix4(1) = sign(1D0, DBLE(d))
	ix4(2) = -ix4(1)

	x1(1) = x4(1)/r24
	x1(2) = x4(2)/r24
	ix1(1) = sign(1D0, ix4(1)*DBLE(r24))
	ix1(2) = -ix1(1)

	c = cln(DCMPLX(k13), -1D0)
	l4(1) = c + cln((q12 + q24*x4(1))/d, DBLE(q24*ix4(1)/d))
	l4(2) = c + cln((q12 + q24*x4(2))/d, DBLE(q24*ix4(2)/d))

	res = (
     &    xspence(x4, ix4, r14, ir14) +
     &    xspence(x4, ix4, 1/r14, -ir14) -
     &    xspence(x4, ix4, DCMPLX(k34/k13), -k13) -
     &    xspence(x1, ix1, r12, ir12) -
     &    xspence(x1, ix1, 1/r12, -ir12) +
     &    xspence(x1, ix1, DCMPLX(k23/k13), -k13) -
     &    c2ipi*xetatilde(x4, ix4, 1/r24, -ir24, l4)
     &  )/(m3*m*a*(x4(2) - x4(1)))
	end

************************************************************************

	subroutine D0m2(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m3, m4, p1, p2, p3, p4, p1p2, p2p3
	double precision m, k12, k13, k14, k23, k24, k34
	double complex k12c, k13c, k23c, k24c, k34c
	double complex r14, x4(2)
	double complex a, b, c, disc

	double complex xspence
	external xspence

	double precision imzero(2)
	data imzero /0D0, 0D0/

	m3 = Mx(3)
	m4 = Mx(4)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)
	p4 = Px(4)
	p1p2 = Px(5)
	p2p3 = Px(6)

	k12 = (m3 - p2)/m3
	k12c = AddEps(k12)

	k13 = (m3 - p1p2)/m3
	k13c = AddEps(k13)

	k23 = -p1/m3
	k23c = AddEps(k23)

	m = sqrt(m3*m4)
	k24 = (m4 - p2p3)/m
	k24c = AddEps(k24)/k12c
	k34 = (m4 - p4)/m
	k34c = AddEps(k34)/k13c
	k14 = (m3 + m4 - p3)/m
	r14 = k2r(k14)
	r14 = r14*DCMPLX(1D0, sign(eps, DBLE(1/r14 - r14)))

	a = k34*k24 - k23
	b = k13*k24 + k12*k34 - k14*k23
	c = k13*k12 - k23*(1 - cI*eps)
	disc = sqrt(b**2 - 4*a*c)
	x4(1) = .5D0/a*(b - disc)
	x4(2) = .5D0/a*(b + disc)
	if( abs(x4(1)) .gt. abs(x4(2)) ) then
	  x4(2) = c/(a*x4(1))
	else
	  x4(1) = c/(a*x4(2))
	endif

	res = (
     &    xspence(x4, imzero, r14, 0D0) +
     &    xspence(x4, imzero, 1/r14, 0D0) -
     &    xspence(x4, imzero, k34c, 0D0) -
     &    xspence(x4, imzero, k24c, 0D0) +
     &    (log(x4(2)) - log(x4(1)))*
     &      (log(k12c) + log(k13c) - log(k23c))
     &  )/(m3*m*a*(x4(2) - x4(1)))
	end

************************************************************************

	subroutine D0m1(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(13)
	integer perm

#include "lt.h"

	double precision m4, k12, k13, k14, k23, k24, k34
	double complex k12c, k13c, k14c, k23c, k24c, k34c
	double precision a, b
	double complex c, disc, x4(2)

	double complex xspence
	external xspence

	double precision imzero(2)
	data imzero /0D0, 0D0/

	m4 = Mx(4)
	k12 = (m4 - Px(3))/m4
	k12c = AddEps(k12)
	k13 = (m4 - Px(4))/m4
	k13c = AddEps(k13)
	k14 = (m4 - Px(6))/m4
	k14c = AddEps(k14)
	k23 = -Px(5)/m4
	k23c = AddEps(k23)
	k24 = -Px(2)/m4
	k24c = AddEps(k24)/k12c
	k34 = -Px(1)/m4
	k34c = AddEps(k34)/k13c

	a = k34*k24
	b = k13*k24 + k12*k34 - k14*k23
	c = k13*k12 - k23*(1 - cI*eps)
	disc = sqrt(b*b - 4*a*c)
	x4(1) = .5D0/a*(b - disc)
	x4(2) = .5D0/a*(b + disc)
	if( abs(x4(1)) .gt. abs(x4(2)) ) then
	  x4(2) = c/(a*x4(1))
	else
	  x4(1) = c/(a*x4(2))
	endif

	res = (
     &    xspence(x4, imzero, k14c, 0D0) -
     &    xspence(x4, imzero, k34c, 0D0) -
     &    xspence(x4, imzero, k24c, 0D0) +
     &    (log(x4(2)) - log(x4(1)))*
     &      (log(k12c) + log(k13c) - log(k23c))
     &  )/(m4**2*a*(x4(2) - x4(1)))
	end

************************************************************************

	subroutine D0m0(res, xpi)
	implicit none
	double complex res
	double precision xpi(13)

#include "lt.h"

	double precision m2, k12, k13, k14, k23, k24, k34
	double complex k12c, k13c, k14c, k23c, k24c, k34c
	double precision a, b
	double complex c, disc, x4(2)

	double complex xspence
	external xspence

	double precision imzero(2)
	data imzero /0D0, 0D0/

	m2 = abs(P(6))
	k12 = -P(1)/m2
	k12c = AddEps(k12)
	k13 = -P(5)/m2
	k13c = AddEps(k13)
	k14 = -P(4)/m2
	k14c = AddEps(k14)
	k23 = -P(2)/m2
	k23c = AddEps(k23)
	k24 = -P(6)/m2
	k24c = AddEps(k24)/k12c
	k34 = -P(3)/m2
	k34c = AddEps(k34)/k13c

	a = k34*k24
	b = k13*k24 + k12*k34 - k14*k23
	c = k13*k12 + cI*eps*k23
	disc = sqrt(b*b - 4*a*c)
	x4(1) = .5D0/a*(b - disc)
	x4(2) = .5D0/a*(b + disc)
	if( abs(x4(1)) .gt. abs(x4(2)) ) then
	  x4(2) = c/(a*x4(1))
	else
	  x4(1) = c/(a*x4(2))
	endif

	res = (
     &    (log(x4(2)) - log(x4(1)))*
     &      (-.5D0*(log(x4(2)) + log(x4(1))) +
     &        log(k12c) + log(k13c) - log(k23c) - log(k14c)) -
     &    xspence(x4, imzero, k34c, 0D0) -
     &    xspence(x4, imzero, k24c, 0D0)
     &  )/(m2**2*a*(x4(2) - x4(1)))
	end

************************************************************************

	double complex function xspence(z1, im1, z2, im2)
	implicit none
	double complex z1(2), z2
	double precision im1(2), im2

#include "lt.h"

	double complex cspence
	external cspence

	xspence = cspence(z1(2), im1(2), z2, im2) -
     &    cspence(z1(1), im1(1), z2, im2)
	end

************************************************************************

	double complex function cspence(z1, im1, z2, im2)
	implicit none
	double complex z1, z2
	double precision im1, im2

#include "lt.h"

	double complex cln, spence
	integer eta
	external cln, spence, eta

	double complex z12
	double precision im12
	integer etas

	z12 = z1*z2
	im12 = im2*sign(1D0, DBLE(z1))
	if( DBLE(z12) .gt. .5D0 ) then
	  cspence = spence(1 - z12, 0D0)
	  etas = eta(z1, im1, z2, im2, im12)
	  if( etas .ne. 0 ) cspence = cspence +
     &      etas*cln(1 - z12, -im12)*c2ipi
	else if( abs(z12) .lt. 1D-4 ) then
	  cspence = pi6
	  if( abs(z12) .gt. 1D-14 ) cspence = cspence -
     &      spence(z12, 0D0) +
     &      (cln(z1, im1) + cln(z2, im2))*z12*
     &        (1 + z12*(.5D0 + z12*(1/3D0 + z12/4D0)))
	else
	  cspence = pi6 - spence(z12, 0D0) -
     &      (cln(z1, im1) + cln(z2, im2))*cln(1 - z12, 0D0)
	endif
	end

************************************************************************

	double complex function xeta(z1, im1, z2, im2, im12, l1)
	implicit none
	double complex z1(2), z2, l1(2)
	double precision im1(2), im2, im12

#include "lt.h"

	integer eta
	external eta

	xeta = l1(2)*eta(z1(2), im1(2), z2, im2, im12) -
     &    l1(1)*eta(z1(1), im1(1), z2, im2, im12)
	end

************************************************************************

	double complex function xetatilde(z1, im1, z2, im2, l1)
	implicit none
	double complex z1(2), z2, l1(2)
	double precision im1(2), im2

#include "lt.h"

	integer etatilde
	external etatilde

	xetatilde = l1(2)*etatilde(z1(2), im1(2), z2, im2) -
     &    l1(1)*etatilde(z1(1), im1(1), z2, im2)
	end

************************************************************************

	integer function etatilde(c1, im1x, c2, im2x)
	implicit none
	double complex c1, c2
	double precision im1x, im2x

	double precision im1, im2

	integer eta
	external eta

	im1 = DIMAG(c1)
	if( im1 .eq. 0 ) im1 = im1x
	im2 = DIMAG(c2)
	if( im2 .ne. 0 ) then
	  etatilde = eta(c1, im1x, c2, 0D0, 0D0)
	else if( DBLE(c2) .gt. 0 ) then
	  etatilde = 0
	else if( im1 .gt. 0 .and. im2x .gt. 0 ) then
	  etatilde = -1
	else if( im1 .lt. 0 .and. im2x .lt. 0 ) then
	  etatilde = 1
	else
	  etatilde = 0
#ifdef WARNINGS
	  if( im1 .eq. 0 .and. DBLE(c1) .lt. 0 .or.
     &        im2x .eq. 0 .and. DBLE(c1*c2) .lt. 0 )
     &      print *, "etatilde not defined"
#endif
	endif
	end

