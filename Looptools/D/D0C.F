* D0C.F
* the scalar four-point function for complex parameters
* this file is part of LoopTools
* last modified 24 Jan 06 th

#include "defs.h"


	double complex function D0C(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	double complex p1, p2, p3, p4, p1p2, p2p3
	double complex m1, m2, m3, m4

#include "lt.h"

	double complex cpi(13)
	integer ier

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
	cpi(11) = 0
	cpi(12) = 0
	cpi(13) = 0

	ier = 0
	call ffcd0(D0C, cpi, ier)

	if( ier .gt. warndigits ) then
	  print *, "D0C lost ", ier, " digits"
	  print *, "p1   =", p1
	  print *, "p2   =", p2
	  print *, "p3   =", p3
	  print *, "p4   =", p4
	  print *, "p1p2 =", p1p2
	  print *, "p2p3 =", p2p3
	  print *, "m1   =", m1
	  print *, "m2   =", m2
	  print *, "m3   =", m3
	  print *, "m4   =", m4
	endif
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

