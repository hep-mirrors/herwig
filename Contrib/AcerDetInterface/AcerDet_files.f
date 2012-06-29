      subroutine acerdet_files
      implicit none
      integer nwpawc
      real hmemor
      PARAMETER (NWPAWC = 5 000 000)
      COMMON /PAWC/ HMEMOR(NWPAWC)
      INCLUDE "acdnout.inc"
C--unit numbers for the inpt and output files
      NINP  = 16
      NOUT  = 11
C--open the files
      OPEN(NINP   ,file='acerdet.dat',status='old')
      OPEN(NOUT   ,file='acerdet.out')
C--hbook 
      CALL HLIMIT(NWPAWC)
C------initialize histo output
      CALL HOUTPU(NOUT)
      end
