* BcoeffC.F
* invoke the two-point tensor coefficients
* this file is part of LoopTools
* last modified 24 Aug 09 th

#include "defs.h"


	subroutine BcoeffC(para, B, ldpara)
	implicit none
	integer ldpara
	double complex para(ldpara,Pbb)
	double complex B(Nbb)

#include "lt.h"

	integer ier(Nbb), i
	logical ini

	character*5 name(Nbb)
	data name /"bb0", "bb1", "bb11", "bb00", "bb001", "bb111",
     &    "dbb0", "dbb1", "dbb00", "dbb11"/

	if( lambda .lt. 0 ) then
	  do i = 1, Nbb
	    B(i) = 0
	  enddo
	  if( lambda .eq. -1 ) then
	    B(bb0) = 1
	    B(bb1) = -.5D0
	    B(bb00) = -(para(1,1) - 3*(para(1,2) + para(1,3)))/12D0
	    B(bb11) = 1/3D0
	    B(bb001) = -(para(1,1) - 2*para(1,2) - 4*para(1,3))/24D0
	    B(bb111) = -.25D0
	    B(dbb00) = -1/12D0
	  endif
	  return
	endif

	call BcoeffaC(para, B, ldpara, ier)

	ini = .TRUE.
	do i = 1, Nbb
	  if( ier(i) .gt. warndigits ) then
	    if( ini ) then
	      print *, "Loss of digits in BgetC for:"
	      call DumpParaC(2, para, ldpara, " ")
	      ini = .FALSE.
	    endif
	    print *, name(i), " claims ", ier(i), "lost digits"
	  endif
	enddo
	end

