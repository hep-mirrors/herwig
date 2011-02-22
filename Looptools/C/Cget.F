* Cget.F
* the three-point tensor coefficients
* this file is part of LoopTools
* improvements by M. Rauch
* last modified 28 Sep 10 th

#include "defs.h"


	integer function XCget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	integer cachelookup
	external cachelookup, XCcoeff

	integer dummy
	DVAR dvardummy

	DVAR para(Pcc)

	para(1) = p1
	para(2) = p2
	para(3) = p1p2
	para(4) = m1
	if( abs(para(4)) .lt. minmass ) para(4) = 0
	para(5) = m2
	if( abs(para(5)) .lt. minmass ) para(5) = 0
	para(6) = m3
	if( abs(para(6)) .lt. minmass ) para(6) = 0
	XCget = cachelookup(para, Cval(1,0), XCcoeff, RC*Pcc, Ncc)
	end

************************************************************************

	double complex function XC0i(i, p1, p2, p1p2, m1, m2, m3)
	implicit none
	integer i
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	integer XCget
	external XCget

	integer b

	b = XCget(p1, p2, p1p2, m1, m2, m3)
	XC0i = Cval(i,b)
	end

************************************************************************

	subroutine XCcoeff(para, C, ldpara)
	implicit none
	integer ldpara
	DVAR para(ldpara,Pcc)
	double complex C(Ncc)

#include "lt.h"

	integer XBget
	double complex XC0
	external XBget, XC0

	DVAR p1, p2, p1p2, m1, m2, m3
	DVAR f1, f2
	QVAR G(2,2)
	double complex bsum, b1sum, b00sum, b11sum, in(2)
	integer B12, B23, B13
	logical dump

#ifdef SOLVE_EIGEN
	QVAR Ginv(2,2)
#define SOLVE_SETUP XInverse(2, G,2, Ginv,2)
#define SOLVE(b) XSolve(2, G,2, Ginv,2, b)
#else
	integer perm(2)
#define IN(i) in(perm(i))
#define SOLVE_SETUP XDecomp(2, G,2, perm)
#define SOLVE(b) XSolve(2, G,2, b)
#endif

#ifdef COMPLEXPARA
	if( abs(DIMAG(para(1,1))) +
     &      abs(DIMAG(para(1,2))) +
     &      abs(DIMAG(para(1,3))) .gt. 0 )
     &    print *, "Warning: complex momenta not implemented"
	if( abs(DIMAG(para(1,4))) +
     &      abs(DIMAG(para(1,5))) +
     &      abs(DIMAG(para(1,6))) .eq. 0 ) then
	  call Ccoeff(para, C, 2)
	  return
	endif
#endif

	p1 = para(1,1)
	p2 = para(1,2)
	p1p2 = para(1,3)
	m1 = para(1,4)
	m2 = para(1,5)
	m3 = para(1,6)

	B12 = XBget(p1, m1, m2)
	B23 = XBget(p2, m2, m3)
	B13 = XBget(p1p2, m1, m3)

	serial = serial + 1
	dump = ibits(debugkey, DebugC, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(3, para, ldpara, "Ccoeff")

	f1 = QPREC(m2) - QPREC(m1) - QPREC(p1)
	f2 = QPREC(m3) - QPREC(m1) - QPREC(p1p2)

	G(1,1) = 2*QPREC(p1)
	G(2,2) = 2*QPREC(p1p2)
	G(1,2) = QPREC(p1) + QPREC(p1p2) - QPREC(p2)
	G(2,1) = G(1,2)
	call SOLVE_SETUP

	bsum = Bval(bb0,B23) + Bval(bb1,B23)
	b1sum = Bval(bb1,B23) + Bval(bb11,B23)
	b00sum = Bval(bb00,B23) + Bval(bb001,B23)
	b11sum = Bval(bb11,B23) + Bval(bb111,B23)

	C(cc0) = XC0(p1, p2, p1p2, m1, m2, m3)

	IN(1) = f1*C(cc0) - Bval(bb0,B23) + Bval(bb0,B13)
	IN(2) = f2*C(cc0) - Bval(bb0,B23) + Bval(bb0,B12)
	call SOLVE(in)
	C(cc1) = in(1)
	C(cc2) = in(2)

	C(cc00) = .5D0*(m1*C(cc0) -
     &    .5D0*(f1*C(cc1) + f2*C(cc2) - Bval(bb0,B23) - 1))

	IN(1) = f1*C(cc1) + bsum - 2*C(cc00)
	IN(2) = f2*C(cc1) + bsum + Bval(bb1,B12)
	call SOLVE(in)
	C(cc11) = in(1)
	C(cc12) = in(2)

	IN(1) = f1*C(cc2) - Bval(bb1,B23) + Bval(bb1,B13)
	IN(2) = f2*C(cc2) - Bval(bb1,B23) - 2*C(cc00)
	call SOLVE(in)
	C(cc12) = .5D0*(C(cc12) + in(1))
	C(cc22) = in(2)

	C(cc001) = 1/3D0*(m1*C(cc1) -
     &    .5D0*(f1*C(cc11) + f2*C(cc12) + bsum + 1/3D0))
	C(cc002) = 1/3D0*(m1*C(cc2) -
     &    .5D0*(f1*C(cc12) + f2*C(cc22) - Bval(bb1,B23) + 1/3D0))

	bsum = bsum + b1sum
	IN(1) = f1*C(cc11) - bsum - 4*C(cc001)
	IN(2) = f2*C(cc11) - bsum + Bval(bb11,B12)
	call SOLVE(in)
	C(cc111) = in(1)
	C(cc112) = in(2)

	IN(1) = f1*C(cc22) - Bval(bb11,B23) + Bval(bb11,B13)
	IN(2) = f2*C(cc22) - Bval(bb11,B23) - 4*C(cc002)
	call SOLVE(in)
	C(cc122) = in(1)
	C(cc222) = in(2)

	C(cc0000) = 1/4D0*(m1*C(cc00) -
     &    .5D0*(f1*C(cc001) + f2*C(cc002) - Bval(bb00,B23) -
     &          (m1 + m2 + m3 - .25D0*(p1 + p2 + p1p2))/6D0))

	IN(1) = f1*C(cc001) + b00sum - 2*C(cc0000)
	IN(2) = f2*C(cc001) + b00sum + Bval(bb001,B12)
	call SOLVE(in)
	C(cc0011) = in(1)
	C(cc0012) = in(2)

	IN(1) = f1*C(cc002) - Bval(bb001,B23) + Bval(bb001,B13)
	IN(2) = f2*C(cc002) - Bval(bb001,B23) - 2*C(cc0000)
	call SOLVE(in)
	C(cc0012) = .5D0*(C(cc0012) + in(1))
	C(cc0022) = in(2)

	bsum = bsum + b1sum + b11sum
	IN(1) = f1*C(cc111) + bsum - 6*C(cc0011)
	IN(2) = f2*C(cc111) + bsum + Bval(bb111,B12)
	call SOLVE(in)
	C(cc1111) = in(1)
	C(cc1112) = in(2)

	IN(1) = f1*C(cc222) - Bval(bb111,B23) + Bval(bb111,B13)
	IN(2) = f2*C(cc222) - Bval(bb111,B23) - 6*C(cc0022)
	call SOLVE(in)
	C(cc1222) = in(1)
	C(cc2222) = in(2)

	IN(1) = f1*C(cc122) + b11sum - 2*C(cc0022)
	IN(2) = f2*C(cc122) + b11sum - 4*C(cc0012)
	call SOLVE(in)
	C(cc1122) = in(1)
	C(cc1222) = .5D0*(C(cc1222) + in(2))

	if( dump ) call XDumpCoeff(3, C)
	end

