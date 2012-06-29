      DOUBLE PRECISION FUNCTION	mrstalphas(T,alambda,flavor,
     &     qsct,qsdt,iord)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI/3.14159/
      DATA TOL/.0005/
      ITH=0
      TT=T
      qsdtt=qsdt/4.
      qsctt=qsct/4.
      AL=ALAMBDA
      AL2=AL*AL
      FLAV=4.
      QS=AL2*dEXP(T)
      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=dlog(qs/al2)
          tt=t
      endif

      IF(QS.gt.QSCTT) GO	TO 12  
      IF(QS.lt.QSDTT) GO	TO 312  
   11 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)1,1,2
c     IF(IORD)2,2,2	!TAKE CARE !!
    1 CONTINUE
      MRSTALPHAS=4.*PI/B0/T
      RETURN
    2 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)
    5 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)3,3,4
    3 CONTINUE
      MRSTALPHAS=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (13,14,15) ITH
    4 CONTINUE
      AS=AS2
      GO TO 5
   12 ITH=1
      T=dLOG(QSCTT/AL2)
      GO TO 11
   13 ALFQC4=MRSTALPHAS
      FLAV=5.   
      ITH=2
      GO TO 11
   14 ALFQC5=MRSTALPHAS
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=MRSTALPHAS
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      MRSTALPHAS=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)31,31,32
c     IF(IORD)32,32,32	!TAKE CARE !!
   31 CONTINUE
      MRSTALPHAS=4.*PI/B0/T
      RETURN
   32 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)
   35 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)33,33,34
   33 CONTINUE
      MRSTALPHAS=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
   34 CONTINUE
      AS=AS2
      GO TO 35
  312 ITH=1
      T=dLOG(QSDTT/AL2)
      GO TO 311
  313 ALFQC4=MRSTALPHAS
      FLAV=3.   
      ITH=2
      GO TO 311
  314 ALFQC3=MRSTALPHAS
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=MRSTALPHAS
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      MRSTALPHAS=1./ALFINV
      RETURN
      END
