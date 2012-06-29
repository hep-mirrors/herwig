* C0.F
* the scalar three-point function
* this file is part of LoopTools
* last modified 21 Dec 10 th

#include "defs.h"

#define legs 3
#define M(i) xpi(i)
#define P(i) xpi(i+legs)


	double complex function C0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double precision p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	external C0soft, C0coll, C0softDR, C0collDR

	double complex res(0:1)
	double precision xpi(6)
	integer key, ier

	M(1) = m1
	M(2) = m2
	M(3) = m3
	P(1) = p1
	P(2) = p2
	P(3) = p1p2

	if( lambda .le. 0 ) then
	  call CDispatch(C0, xpi, C0softDR, C0collDR)
	  return
	endif

	ier = 0
	key = ibits(versionkey, KeyC0, 2)

	if( key .ne. 1 ) then
	  call ffxc0(res(0), xpi, ier)
	  if( ier .gt. warndigits ) then
	    ier = 0
	    call ffxc0r(res(0), xpi, ier)
	    if( ier .gt. warndigits ) key = ior(key, 2)
	  endif
	endif

	if( key .ne. 0 ) then
	  call CDispatch(res(1), xpi, C0soft, C0coll)
	  if( key .gt. 1 .and.
     &        abs(res(0) - res(1)) .gt. maxdev*abs(res(0)) ) then
	    print *, "Discrepancy in C0:"
	    print *, "  p1   =", p1
	    print *, "  p2   =", p2
	    print *, "  p1p2 =", p1p2
	    print *, "  m1   =", m1
	    print *, "  m2   =", m2
	    print *, "  m3   =", m3
	    print *, "C0 a   =", res(0)
	    print *, "C0 b   =", res(1)
	    if( ier .gt. errdigits ) res(0) = res(1)
	  endif
	endif

	C0 = res(iand(key, 1))
	end

************************************************************************
* adapter code for C++

	subroutine c0sub(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	double complex res
	double precision p1, p2, p1p2, m1, m2, m3

	double complex C0
	external C0

	res = C0(p1, p2, p1p2, m1, m2, m3)
	end

************************************************************************

	subroutine CDispatch(res, xpi, soft, coll)
	implicit none
	double complex res
	double precision xpi(6)
	external soft, coll

#include "lt.h"

	integer i, z, c, perm

	double complex C0p3, C0p2, C0p1, C0p0
	external C0p3, C0p2, C0p1, C0p0

	integer p123, p231, p312
	parameter (p123 = O'123', p231 = O'231', p312 = O'312')

	integer xpiperm(3)
	data xpiperm /p123, p231, p312/

#define Px(j) P(ibits(perm,3*(3-j),3))
#define Mx(j) M(ibits(perm,3*(3-j),3))

555	z = 0
	c = 0
	do i = 1, 3
	  perm = xpiperm(i)
	  if( abs(Mx(1)) .lt. eps ) then
	    if( abs(Px(1) - Mx(2)) +
     &          abs(Px(3) - Mx(3)) .lt. acc ) then
	      if( DEBUGLEVEL .gt. 0 )
     &          print '(A,O3)', "soft C0, perm = ", perm
	      call soft(res, xpi, perm)
	      return
	    endif
	    if( abs(Px(1)) + abs(Mx(2)) .lt. eps ) c = perm
	  endif
	  if( abs(P(i)) .lt. eps ) z = z + 1
	enddo

	if( c .ne. 0 ) then
	  if( DEBUGLEVEL .gt. 0 )
     &      print '(A,O3)', "collinear C0, perm = ", perm
	  call coll(res, xpi, c)
	  if( res .eq. c ) goto 555
	  return
	endif

	goto (1, 2, 3) z

	res = C0p3(xpi, p123) + C0p3(xpi, p231) + C0p3(xpi, p312)
	return

1	res = C0p2(xpi, p123) + C0p2(xpi, p231) + C0p2(xpi, p312)
	return

2	res = C0p1(xpi, p123) + C0p1(xpi, p231) + C0p1(xpi, p312)
	return

3	res = C0p0(xpi)
	end

************************************************************************
* the following routines are adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* they are used for double-checking the results of FF

* all mom-squares != 0

	double complex function C0p3(xpi, perm)
	implicit none
	double precision xpi(6)
	integer perm

#include "lt.h"

	double precision m1, m2, m3, p1, p2, p3
	double precision m12, m13, m23, a2, n
	double complex a, b, c
	double complex y1, y2, y3, y4, x1, x2, x3, x4, z3, z4
	integer z3z4, x1z3, x3z3, x2z4, x4z4

	double complex spence
	integer eta
	external spence, eta

	m1 = Mx(1)
	m2 = Mx(2)
	m3 = Mx(3)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)

	m12 = m1 - m2
	m13 = m1 - m3
	m23 = m2 - m3

	a2 = (p1 - p2 - p3)**2 - 4*p2*p3
	a = sqrt(DCMPLX(a2))
	n = .5D0/p1
	c = (p1*(p1 - p2 - p3 - m13 - m23) - m12*(p2 - p3))/a

	y1 = n*(c + (p1 - m12))
	y4 = n*(c - (p1 - m12))
	if( abs(y1) .lt. abs(y4) )
     &    y1 = ((p1 - m12)*(p2 + m23)*p3 -
     &      (p1 - p2 - m13)*(m12*p2 + m23*p1))/(a2*p1*y4)

	y2 = n*(c - (p1 + m12))
	y4 = n*(c + (p1 + m12))
	if( abs(y2) .lt. abs(y4) )
     &    y2 = ((p1 + m12)*(p3 + m13)*p2 +
     &      (p1 - p3 - m23)*(m12*p3 - m13*p1))/(a2*p1*y4)

	b = sqrt(DCMPLX((p1 - m12)**2 - 4*p1*m2))
	y3 = n*(c + b)
	y4 = n*(c - b)
	n = (p1*((p1 - p2 - p3)*m3 + p2*p3 + m13*m23) +
     &       p2*((p2 - p3 - p1)*m1 + m12*m13) +
     &       p3*((p3 - p1 - p2)*m2 - m12*m23))/(a2*p1)
	if( abs(y3) .lt. abs(y4) ) then
	  y3 = n/y4
	else
	  y4 = n/y3
	endif
	c = a*b
	if( c .ne. 0 ) then
	  y3 = y3 + cI*eps*abs(c*y3)/c
	  y4 = y4 - cI*eps*abs(c*y4)/c
	else
	  y3 = y3*(1 + cI*eps)
	  y4 = y4*(1 - cI*eps)
	endif

	C0p3 = spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -
     &           spence(y1/y3, 0D0) - spence(y1/y4, 0D0)

	if( DIMAG(a) .ne. 0 ) then
	  c = abs(b)/b*eps*cI
	  x1 = -.5D0/p1*( p1 - m12 + b) - c
	  x2 = -.5D0/p1*( p1 - m12 - b) - c
	  x3 = -.5D0/p1*(-p1 - m12 + b) - c
	  x4 = -.5D0/p1*(-p1 - m12 - b) - c
	  z3 = 1/y3
	  z4 = 1/y4
	  z3z4 = eta(z3, 0D0, z4, 0D0, 0D0)
	  x1z3 = eta(x1, 0D0, z3, 0D0, 0D0)
	  x3z3 = eta(x3, 0D0, z3, 0D0, 0D0)
	  x2z4 = eta(x2, 0D0, z4, 0D0, 0D0)
	  x4z4 = eta(x4, 0D0, z4, 0D0, 0D0)
	  c = log(y1)*(eta(x1, 0D0, x2, 0D0, 0D0) +
     &                 z3z4 - x1z3 - x2z4) -
     &        log(y2)*(eta(x3, 0D0, x4, 0D0, 0D0) +
     &                 z3z4 - x3z3 - x4z4) +
     &        log(y3)*(x1z3 - x3z3) +
     &        log(y4)*(x2z4 - x4z4)
	  if( DIMAG(a) .gt. 0 .and. p1 .lt. 0 ) c = c - log(y1/y2)
	  C0p3 = C0p3 + c2ipi*c
	endif

	C0p3 = C0p3/a
	end

************************************************************************
* one mom-square zero

	double complex function C0p2(xpi, perm)
	implicit none
	double precision xpi(6)
	integer perm

#include "lt.h"

	double precision m1, m2, m3, p1, p2, p3
	double precision m12, m23, m13, a, c, y1, y2
	double complex b, y3, y4

	double complex spence
	external spence

	if( abs(Px(1)) .lt. eps ) then
	  C0p2 = 0
	  return
	endif

	m1 = Mx(1)
	m2 = Mx(2)
	m3 = Mx(3)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)

	m12 = m1 - m2
	m23 = m2 - m3
	m13 = m1 - m3

	if( abs(p3) .lt. eps ) then
	  a = p1 - p2
	  y1 = -2*p1*(m13 - a)
	  y2 = -2*p1*m13
	else
	  a = p3 - p1
	  y1 = -2*p1*m23
	  y2 = -2*p1*(m23 + a)
	endif

	c = p1*(p1 - p2 - p3 - m13 - m23) - m12*(p2 - p3)
	b = a*sqrt(DCMPLX((p1 - m12)**2 - 4*p1*m2))
	y3 = c + b
	y4 = c - b
	c = 4*p1*(
     &    p1*((p1 - p2 - p3)*m3 + p2*p3 + m13*m23) +
     &    p2*((p2 - p3 - p1)*m1 + m12*m13) +
     &    p3*((p3 - p1 - p2)*m2 - m12*m23) )
	if( abs(y3) .lt. abs(y4) ) then
	  y3 = c/y4
	else
	  y4 = c/y3
	endif
	c = a/p1
	if( c .ne. 0 ) then
	  y3 = y3 + cI*eps*abs(c*y3)/c
	  y4 = y4 - cI*eps*abs(c*y4)/c
	else
	  y3 = y3*(1 + cI*eps)
	  y4 = y4*(1 - cI*eps)
	endif

	C0p2 = (spence(y2/y3, 0D0) + spence(y2/y4, 0D0) -
     &            spence(y1/y3, 0D0) - spence(y1/y4, 0D0))/a
	end

************************************************************************
* two mom-squares zero

	double complex function C0p1(xpi, perm)
	implicit none
	double precision xpi(6)
	integer perm

#include "lt.h"

	double precision m1, m2, m3, p1, p2, p3
	double precision m12, m23, m13, c, y1, y2
	double complex b, y3, y4

	double complex spence
	external spence

	if( abs(Px(1)) .lt. eps ) then
	  C0p1 = 0
	  return
	endif

	m1 = Mx(1)
	m2 = Mx(2)
	m3 = Mx(3)
	p1 = Px(1)
	p2 = Px(2)
	p3 = Px(3)

	m12 = m1 - m2
	m23 = m2 - m3
	m13 = m1 - m3

	C0p1 = 0

	if( abs(m13) .gt. acc ) then
	  y1 = m23 - p1
	  y2 = m23
	  c = m23 + p1*m3/m13
	  y3 = c - cI*eps*sign(1D0, p1/m13)*abs(c)
	  C0p1 = spence(y1/y3, 0D0) - spence(y2/y3, 0D0)
	endif

	y1 = -2*p1*m23
	y2 = -2*p1*(m23 - p1)

	c = p1*(p1 - m13 - m23)
	b = p1*sqrt(DCMPLX((p1 - m12)**2 - 4*p1*m2))
	y3 = c - b
	y4 = c + b
	c = 4*p1**2*(p1*m3 + m13*m23)
	if( abs(y3) .lt. abs(y4) ) then
	  y3 = c/y4
	else
	  y4 = c/y3
	endif
	y3 = y3 - cI*eps*abs(y3)
	y4 = y4 + cI*eps*abs(y4)

	C0p1 = (C0p1 +
     &    spence(y1/y3, 0D0) + spence(y1/y4, 0D0) -
     &    spence(y2/y3, 0D0) - spence(y2/y4, 0D0))/p1
	end

************************************************************************

	double complex function C0p0(xpi)
	implicit none
	double precision xpi(6)

#include "lt.h"

	double precision m1, m2, m3
	double precision m12, m23, m13

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)

	m12 = m1 - m2
	m23 = m2 - m3
	m13 = m1 - m3

	if( abs(m23) .lt. acc ) then
	  if( abs(m13) .lt. acc ) then
	    C0p0 = -.5D0/m1
	  else
	    C0p0 = (m13 - m1*log(m1/m3))/m13**2
	  endif
	else
	  if( abs(m12) .lt. acc ) then
	    C0p0 = (-m23 + m3*log(m2/m3))/m23**2
	  else if( abs(m13) .lt. acc ) then
	    C0p0 = (m23 - m2*log(m2/m3))/m23**2
	  else
	    C0p0 = m3/(m13*m23)*log(m1/m3) - m2/(m12*m23)*log(m1/m2)
	  endif
	endif
	end

************************************************************************

	subroutine C0soft(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(6)
	integer perm

#include "lt.h"

	double complex ln, spence
	external ln, spence

	double precision s, m1, m2
	double precision a, h1, h2, h3, ps
	double complex ls

	logical ini
	data ini /.FALSE./

	s = Px(2)
	m1 = Px(1)
	m2 = Px(3)

	a = sqrt(4*m1*m2)
	if( abs(a) .lt. eps ) then
	  ps = max(minmass, 1D-14)
	  if( abs(m1) .lt. eps ) m1 = ps
	  if( abs(m2) .lt. eps ) m2 = ps
	  if( .not. ini ) then
	    print *, "collinear-divergent C0, using mass cutoff ", ps
	    ini = .TRUE.
	  endif
	endif

	if( abs(s) .lt. acc ) then
	  if( abs(m1 - m2) .lt. acc ) then
	    res = -.5D0*log(m1/lambda)/m1
	  else
	    res = -.25D0*log(m2*m1/lambda**2)*
     &        log(m1/m2)/(m1 - m2)
	  endif
	  return
	endif

	ps = s - m1 - m2
	a = (ps - a)*(ps + a)
	if( a .lt. 0 ) then
	  print *, "C0soft: complex square root not implemented"
	  a = 0
	endif
	a = sqrt(a)
	if( ps .le. 0 ) then
	  h1 = .5D0*(a - ps)
	else
	  h1 = -2*m1*m2/(a + ps)
	endif
	ps = s - m1 + m2
	if( ps .le. 0 ) then
	  h2 = .5D0*(a - ps)
	else
	  h2 = -2*s*m2/(a + ps)
	endif
	ps = s + m1 - m2
	if( ps .le. 0 ) then
	  h3 = .5D0*(a - ps)
	else
	  h3 = -2*m1*s/(a + ps)
	endif

	ls = ln(-a/s, -1D0)
	res = (-pi6 +
     &    spence(DCMPLX(h2/a), -1D0) + spence(DCMPLX(h3/a), -1D0) -
     &    .5D0*(ln(-h2/s, -1D0)**2 + ln(-h3/s, -1D0)**2) +
     &    .25D0*(ln(-m1/s, -1D0)**2 + ln(-m2/s, -1D0)**2) -
     &    ls*(ln(-h1/s, -1D0) - ls) +
     &    ln(-lambda/s, -1D0)*ln(h1/sqrt(m1*m2), 1D0))/a
	end

************************************************************************

	subroutine C0coll(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(6)
	integer perm

#include "lt.h"

	logical ini
	data ini /.FALSE./

	Px(1) = max(minmass, 1D-14)
	res = perm

	if( ini ) return
	print *, "collinear-divergent C0, using mass cutoff ", Px(1)
	ini = .TRUE.
	end

************************************************************************

	subroutine C0softDR(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(6)
	integer perm

#include "lt.h"

	double precision s, m1, m2
	double precision m, dm, r
	double complex root, fac, ls, lm, mK, lmK

	double complex Li2omx2, lnrat, cln, spence
	external Li2omx2, lnrat, cln, spence

	s = Px(2)
	m1 = Px(1)
	m2 = Px(3)

	m = sqrt(m1*m2)
	if( abs(m) .lt. eps ) then
	  if( abs(m1) .lt. eps ) then
	    m1 = m2
	    if( abs(m1) .lt. eps ) then
	      if( abs(s) .lt. eps ) then
	        print *, "C0softDR: all scales zero"
	        res = nan
	        return
	      endif
* qltri1
	      if( DEBUGLEVEL .gt. 1 ) print *, "qltri1"
	      if( lambda .eq. -2 ) then
	        res = 1/s
	      else if( lambda .eq. -1 ) then
	        res = lnrat(mudim, -s)/s
	      else
	        res = .5D0*lnrat(mudim, -s)**2/s
	      endif
	      return
	    endif
	  endif
	  if( abs(s - m1) .lt. acc ) then
* qltri5
	    if( DEBUGLEVEL .gt. 1 ) print *, "qltri5"
	    if( lambda .eq. -2 ) then
	      res = 0
	    else if( lambda .eq. -1 ) then
	      res = -.5D0/m1
	    else
	      res = (-.5D0*lnrat(mudim, m1) + 1)/m1
	    endif
	    return
	  endif

* qltri4
	  if( DEBUGLEVEL .gt. 1 ) print *, "qltri4"
	  if( lambda .eq. -2 ) then
	    res = .5D0/(s - m1)
	  else if( lambda .eq. -1 ) then
	    res = (.5D0*lnrat(mudim, m1) + lnrat(m1, m1 - s))/(s - m1)
	  else
	    ls = lnrat(m1, m1 - s)
	    lm = lnrat(mudim, m1)
	    res = (lm*(.25D0*lm + ls) + .5D0*ls**2 + pi12 -
     &        spence(DCMPLX(s/(s - m1)), 0D0))/(s - m1)
	  endif
	  return
	endif

	if( lambda .eq. -2 ) then
	  res = 0
	  return
	endif

* qltri6
	if( DEBUGLEVEL .gt. 1 ) print *, "qltri6"
	dm = sqrt(m1) - sqrt(m2)
	r = s - dm**2
	root = sqrt(DCMPLX((r - 4*m)/r))
	mK = -4*m/(r*(1 + root)**2)

	if( abs(mK - 1) .lt. acc ) then
	  if( lambda .eq. -1 ) then
	    res = .5D0/m
	  else
	    res = 0
	    if( abs(m1 - m2) .gt. acc )
     &        res = 2 + .5D0*(sqrt(m1) + sqrt(m2))/dm*log(m2/m1)
	    res = .5D0/m*(log(mudim/m) - res)
	  endif
	  return
	endif

	lmK = cln(mK, 1D0)
	fac = 1/(r*root)

	if( lambda .eq. -1 ) then
	  res = fac*lmK
	else
	  res = fac*( lmK*(.5D0*lmK + log(mudim/m)) -
     &      .125D0*log(m1/m2)**2 +
     &      Li2omx2(mK, 1D0, mK, 1D0) -
     &      Li2omx2(mK, 1D0, DCMPLX(sqrt(m1/m2)), 0D0) -
     &      Li2omx2(mK, 1D0, DCMPLX(sqrt(m2/m1)), 0D0) )
	endif
	end

************************************************************************

	subroutine C0collDR(res, xpi, perm)
	implicit none
	double complex res
	double precision xpi(6)
	integer perm

#include "lt.h"

	double precision s1, s2, m
	double precision m1, m2, r
	double complex l1, l2, lm

	double complex Li2omrat, lnrat
	external Li2omrat, lnrat

	if( lambda .eq. -2 ) then
	  res = 0
	  return
	endif

	m = Mx(3)
	s1 = Px(2)
	s2 = Px(3)

	if( abs(m) .lt. eps ) then
* qltri2
	  if( DEBUGLEVEL .gt. 1 ) print *, "qltri2"
	  r = .5D0*(s2 - s1)/s1
	  if( abs(r) .lt. acc ) then
	    if( lambda .eq. -1 ) then
	      res = (1 - r*mudim/s1)/s1
	    else
	      res = (lnrat(mudim, -s1)*(1 - r) - r)/s1
	    endif
	    return
	  endif
	  l1 = lnrat(mudim, -s1)
	  l2 = lnrat(mudim, -s2)
	  res = (l1 - l2)/(s1 - s2)
	  if( lambda .ne. -1 ) res = .5D0*(l1 + l2)*res
	  return
	endif

* qltri3
	if( DEBUGLEVEL .gt. 1 ) print *, "qltri3"
	m1 = m - s1
	m2 = m - s2
	l1 = lnrat(m1, m)
	l2 = lnrat(m2, m)
	lm = lnrat(mudim, m)

	r = .5D0*(s1 - s2)/m1
	if( abs(r) .lt. acc ) then
	  if( lambda .eq. -1 ) then
	    res = (1 - r)/m1
	  else
	    m = m/s1
	    res = (lm - (m + 1)*(l2 + r) -
     &        r*((m*(m - 2) - 1)*l2 + lm))/m1
	  endif
	  return
	endif

	res = l2 - l1
	if( lambda .ne. -1 )
     &    res = (lm - l1 - l2)*res +
     &      Li2omrat(m1, m) - Li2omrat(m2, m)
	res = res/(s1 - s2)
	end

