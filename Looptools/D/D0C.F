* D0C.F
* the scalar four-point function with complex masses
* this file is part of LoopTools
* last modified 29 Jan 09 th

#include "defs.h"


	double complex function D0C(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double complex p1, p2, p3, p4, p1p2, p2p3
	double complex m1, m2, m3, m4

#include "lt.h"

	double complex cpi(10), res(0:1)
	integer key, ier, ier1

	cpi(1) = m1
	cpi(2) = m2
	cpi(3) = m3
	cpi(4) = m4
	cpi(5) = p1
	cpi(6) = p2
	cpi(7) = p3
	cpi(8) = p4
	cpi(9) = p1p2
	cpi(10) = p2p3

	ier = 0
	key = ibits(versionkey, KeyD0C, 2)

	if( key .ne. 1 ) call ffd0c(res(0), cpi, 0, ier)

	if( key .ne. 0 ) then
	  ier1 = 0
	  call ffd0c(res(1), cpi, 1, ier1)
	  if( key .gt. 1 .and.
     &        abs(res(0) - res(1)) .gt. maxdev*abs(res(0)) ) then
	    print *, "Discrepancy in D0C:"
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
	    print *, "D0C a  =", res(0)
	    print *, "D0C b  =", res(1)
	    if( ier1 .lt. ier .and. ier .gt. errdigits )
     &        res(0) = res(1)
	  endif
	endif

	D0C = res(iand(key, 1))
	end

************************************************************************
* adapter code for C++

	subroutine d0subc(res, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double complex res
	double complex p1, p2, p3, p4, p1p2, p2p3
	double complex m1, m2, m3, m4

	double complex D0C
	external D0C

	res = D0C(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
	end

