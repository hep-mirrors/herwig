      PROGRAM TEST
      DIMENSION XPDFGM(-6:6)
      ISET=2
      IP2=0
      X=0.0860466
      Q2=262.249
      P2=6.67495e-08
      CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
      DO I=-5,5
         PRINT *, XPDFGM(I)
      ENDDO
      END
