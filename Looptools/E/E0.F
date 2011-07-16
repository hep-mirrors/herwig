* E0.F
* the scalar four-point function
* this file is part of LoopTools
* written by M. Rauch
* last modified 8 Dec 10 th

#include "defs.h"


#ifdef COMPLEXPARA
#undef E0b
#define E0b E0C
#else
	double complex function E0(p1, p2, p3, p4, p5, 
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

#include "lt.h"

	double complex E0b
	external E0b

	double complex res(0:1), d0(5)
	double precision xpi(20)
	integer key, ier

	key = ibits(versionkey, KeyE0, 2)

	if( key .ne. 0 ) then
	  ier = 0
	  xpi(1) = m1
	  xpi(2) = m2
	  xpi(3) = m3
	  xpi(4) = m4
	  xpi(5) = m5
	  xpi(6) = p1
	  xpi(7) = p2
	  xpi(8) = p3
	  xpi(9) = p4
	  xpi(10) = p5
	  xpi(11) = p1p2
	  xpi(12) = p2p3
	  xpi(13) = p3p4
	  xpi(14) = p4p5
	  xpi(15) = p5p1

	  call ffxe0(res(1), d0, xpi, ier)
          if( ier .gt. warndigits ) key = ior(key, 2)
	endif

	if( key .ne. 1 ) then
          res(0) = E0b(p1, p2, p3, p4, p5,
     &      p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
          if( key .gt. 1 .and.
     &        abs(res(0) - res(1)) .gt. maxdev*abs(res(0)) ) then
            print *, "Discrepancy in E0:"
            print *, "  p1   =", p1
            print *, "  p2   =", p2
            print *, "  p3   =", p3
            print *, "  p4   =", p4
            print *, "  p5   =", p5
            print *, "  p1p2 =", p1p2
            print *, "  p2p3 =", p2p3
            print *, "  p3p4 =", p3p4
            print *, "  p4p5 =", p4p5
            print *, "  p5p1 =", p5p1
            print *, "  m1   =", m1
            print *, "  m2   =", m2
            print *, "  m3   =", m3
            print *, "  m4   =", m4
            print *, "  m5   =", m5
            print *, "E0 a   =", res(0)
            print *, "E0 b   =", res(1)
          endif
        endif

	E0 = res(iand(key, 1))
	end
#endif

************************************************************************
* adapter code for C++

	subroutine XE0sub(res, p1, p2, p3, p4, p5,
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	double complex res
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

	double complex XE0
	external XE0

	res = XE0(p1, p2, p3, p4, p5,
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	end

************************************************************************

	double complex function E0b(p1, p2, p3, p4, p5, 
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

#include "lt.h"

	QVAR Y(5,5), Yi(5,5), eta(5), detY
	integer i, j

	QVAR Yflat(25), Yiflat(25)
	equivalence (Y, Yflat)
	equivalence (Yi, Yiflat)

	double complex XD0
	external XD0

	Y(1,1) = 2*QPREC(m1)
	Y(2,2) = 2*QPREC(m2)
	Y(3,3) = 2*QPREC(m3)
	Y(4,4) = 2*QPREC(m4)
	Y(5,5) = 2*QPREC(m5)
	Y(1,2) = QPREC(m1) + QPREC(m2) - QPREC(p1)
	Y(2,1) = Y(1,2)
	Y(1,3) = QPREC(m1) + QPREC(m3) - QPREC(p1p2)
	Y(3,1) = Y(1,3)
	Y(1,4) = QPREC(m1) + QPREC(m4) - QPREC(p4p5)
	Y(4,1) = Y(1,4)
	Y(1,5) = QPREC(m1) + QPREC(m5) - QPREC(p5)
	Y(5,1) = Y(1,5)
	Y(2,3) = QPREC(m2) + QPREC(m3) - QPREC(p2)
	Y(3,2) = Y(2,3)
	Y(2,4) = QPREC(m2) + QPREC(m4) - QPREC(p2p3)
	Y(4,2) = Y(2,4)
	Y(2,5) = QPREC(m2) + QPREC(m5) - QPREC(p5p1)
	Y(5,2) = Y(2,5)
	Y(3,4) = QPREC(m3) + QPREC(m4) - QPREC(p3)
	Y(4,3) = Y(3,4)
	Y(3,5) = QPREC(m3) + QPREC(m5) - QPREC(p3p4)
	Y(5,3) = Y(3,5)
	Y(4,5) = QPREC(m4) + QPREC(m5) - QPREC(p4)
	Y(5,4) = Y(4,5)

	do i = 1, 5
	  do j = 1, 25
	    Yiflat(j) = Yflat(j)
	  enddo
	  do j = 1, 5
	    Yi(j,i) = 1
	  enddo
	  call XDet(5, Yi,5, eta(i))
	enddo

	call XDet(5, Y,5, detY)

	E0b = -(
     &    eta(1)*XD0(p2, p3, p4, p5p1, p2p3, p3p4, m2, m3, m4, m5) +
     &    eta(2)*XD0(p1p2, p3, p4, p5, p4p5, p3p4, m1, m3, m4, m5) +
     &    eta(3)*XD0(p1, p2p3, p4, p5, p4p5, p5p1, m1, m2, m4, m5) +
     &    eta(4)*XD0(p1, p2, p3p4, p5, p1p2, p5p1, m1, m2, m3, m5) +
     &    eta(5)*XD0(p1, p2, p3, p4p5, p1p2, p2p3, m1, m2, m3, m4)
     &  )/detY
	end

