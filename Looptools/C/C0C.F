* C0C.F
* the scalar three-point function for complex parameters
* this file is part of LoopTools
* last modified 24 Jan 06 th

#include "defs.h"


	double complex function C0C(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double complex p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	double complex cpi(6)
	integer ier

	cpi(1) = m1
	cpi(2) = m2
	cpi(3) = m3
	cpi(4) = p1
	cpi(5) = p2
	cpi(6) = p1p2
	ier = 0
	call ffcc0(C0C, cpi, ier)

	if( ier .gt. warndigits ) then
	  ier = 0
	  call ffcc0r(C0C, cpi, ier)
	  if( ier .gt. warndigits ) then
	    print *, "C0C claims ", ier, " lost digits"
	    print *, "  p1   =", p1
	    print *, "  p2   =", p2
	    print *, "  p1p2 =", p1p2
	    print *, "  m1   =", m1
	    print *, "  m2   =", m2
	    print *, "  m3   =", m3
	  endif
	endif
	end

************************************************************************
* adapter code for C++

	subroutine c0subc(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	double complex res
	double complex p1, p2, p1p2, m1, m2, m3

	double complex C0C
	external C0C

	res = C0C(p1, p2, p1p2, m1, m2, m3)
	end

