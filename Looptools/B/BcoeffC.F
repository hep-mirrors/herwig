* BcoeffC.F
* invoke the two-point tensor coefficients
* this file is part of LoopTools
* last modified 16 Nov 06 th

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

	call BcoeffaC(para, B, ldpara, ier)

	ini = .TRUE.
	do i = 1, Nbb
	  if( ier(i) .gt. warndigits ) then
	    if( ini ) then
	      print *, "Loss of digits in CBget for:"
	      call DumpParaC(2, para, ldpara, " ")
	      ini = .FALSE.
	    endif
	    print *, name(i), " claims ", ier(i), "lost digits"
	  endif
	enddo
	end

