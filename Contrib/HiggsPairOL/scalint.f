************************************************************************
        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
************************************************************************
        IMPLICIT REAL*8 (A-Z)
	REAL*8 M(4),P(4,4),K(4,4)
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 C04,D04,CSPEN,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
        COMPLEX*16 Z2,B,SC,TC,WP,WM,BS,XS
	INTEGER GEN,I,J

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D04 = D04 + (-1D0)**(I+J)*(
     *		CSPEN(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPEN(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )

	   D04 = D04 + (3D0-2D0*I)*(
     *		 ETAS(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA(RS(1,3),1D0/RS(2,4))
     *		   *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 + (3D0-2D0*I)*(
     *		+ETA(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA(-XX(I),R(1,3)/R(2,4))
     *		  +ETA(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA(R(1,3),1D0/R(2,4))
     *		   *ETA(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       -CSPEN(1D0+SS(1)*X(I,1))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPEN(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPEN(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPEN(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END

************************************************************************
        FUNCTION C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
      COMPLEX*16 C03,CSPEN,ETA,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM
C     REAL*8 KAPPA
      COMPLEX*16 KAPPA
C     KAPPA(A,B,C) = DSQRT(A**2+B**2+C**2-2*(A*B+A*C+B*C))
C     KAPPA(A,B,C) = DSQRT(DABS(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      KAPPA(A,B,C) = CDSQRT(DCMPLX(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      EPS = 1.D-8
      IM = DCMPLX(0.D0,1.D0)
      IEPS = DCMPLX(0.D0,1.D-17)
      PI = 4*DATAN(1.D0)
      XX = 0.D0
C     IF(P1.LT.0.D0.OR.P2.LT.0.D0.OR.P3.LT.0.D0) XX=1.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
      ALPHA = KAPPA(Q10,Q21,Q20)
      ALP(0) = KAPPA(Q21,SM1,SM2)*(1+IEPS*Q21)
      ALP(1) = KAPPA(Q20,SM2,SM0)*(1+IEPS*Q20)
      ALP(2) = KAPPA(Q10,SM0,SM1)*(1+IEPS*Q10)
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
      CDUM=0.D0
      DO I=0,2
       DO J=1,2
        CDUM = CDUM + CSPEN((Y0(I)-1)/Y(I,J)) - CSPEN(Y0(I)/Y(I,J))
        CX = ETA(1-X(I,J),1/Y(I,J))
        IF(CX.NE.0.D0)THEN
         CDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = ETA(-X(I,J),1/Y(I,J))
        IF(CY.NE.0.D0)THEN
         CDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = ETA(-X(I,1),-X(I,2))
       IF(CX.NE.0.D0)THEN
        CDUM = CDUM - CX*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       CY = ETA(Y(I,1),Y(I,2))
       IF(CY.NE.0.D0)THEN
        CDUM = CDUM + CY*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2*PI*IM*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      C03 = CDUM/ALPHA
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPEN,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)/                                                   
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=854513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPEN=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPEN=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPEN=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPEN=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
        END                          
***********************************************************************
        FUNCTION ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END   
***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END 

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                            


      SUBROUTINE FORMFAC(AMQ,  S, T, U, M1, M2, C0AB,C0AC,C0AD,C0BC,C0BD
     &,C0CD,D0ABC,D0BAC,D0ACB)
C--FORM FACTORS FOR LO MATRIX ELEMENTS FOR GG -> HH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMPLEX*16 CQ2,CA5
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
  
      SS=S/AMQ**2
      TT=T/AMQ**2
      UU=U/AMQ**2
      R1=M1**2/AMQ**2
      R2=M2**2/AMQ**2
      RZ=MZ**2/AMQ**2
      T1=TT-R1
      U1=UU-R1
    
      A1 =  - 4*C0AB*SS
 
      A2 = 0
 
      H1 =  - 4*(C0AB*SS - 4*C0AB - 2)
 
      H2 = 0
 
      Z1 = ( - 4*(RZ - SS)*(2*C0AB + 1)*(U1 + T1 + SS))/RZ
   
      Z2 = 0
 
      AA1 = 0
      AA2 = 0
      HH1 = 0
      HH2 = 0
      AH1 = 0
      AH2 = 0

      AA1 = (16*C0AB*SS + 2*C0AC*T1*(U1 + T1 + SS + 2*R1) + 2*C0AD*(
     . - U1*T1 - U1*SS - T1**2 - 2*T1*SS - 2*T1*R1 - SS**2
     . - 2*SS*R1) + 2*C0BC*U1*(U1 + T1 + SS + 2*R1) + 2*C0BD
     .*( - U1**2 - U1*T1 - 2*U1*SS - 2*U1*R1 - T1*SS - SS**2 -
     .2*SS*R1) + 4*D0BAC*SS*( - U1 - T1 - 2*R1) + 2*D0ACB*(
     . - U1**2*T1 - U1*T1**2 - U1*T1*SS - 2*U1*T1*R1 + U1*SS
     .*R1 - 2*U1*SS + T1*SS*R1 - 2*T1*SS + SS**2*R1 + 2*SS*
     .R1**2 - 4*SS*R1) + 4*D0ABC*SS*( - U1 - T1 - 2*R1) + 8*
     .SS)/SS
 
      AA2 = (2*C0AB*SS*(U1**2 + 4*U1*R1 + T1**2 + 4*T1*R1 + 2*SS*R1 +
     .4*R1**2) + 2*C0AC*T1*(U1*R1 + T1**2 + 3*T1*R1 + SS*R1
     . + 2*R1**2) + 2*C0CD*( - U1**3 - U1**2*T1 - 2*U1**2*R1
     . - U1*T1**2 + 2*U1*SS*R1 - T1**3 - 2*T1**2*R1 + 2*T1*
     .SS*R1 + 4*SS*R1**2) + 2*C0AD*( - U1**2*T1 - U1**2*SS
     . - 3*U1*T1*R1 - 3*U1*SS*R1 - T1**2*R1 - 2*T1*SS*R1 - 2
     .*T1*R1**2 - SS**2*R1 - 2*SS*R1**2) + 2*C0BC*U1*(U1**2
     . + 3*U1*R1 + T1*R1 + SS*R1 + 2*R1**2) + 2*C0BD*( - U1
     .**2*R1 - U1*T1**2 - 3*U1*T1*R1 - 2*U1*SS*R1 - 2*U1*R1
     .**2 - T1**2*SS - 3*T1*SS*R1 - SS**2*R1 - 2*SS*R1**2))
      AA2 = (AA2
     . + 2*D0BAC*( - 2*U1**2*T1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*
     .U1*T1*R1 - U1*SS*R1**2 + 2*U1*SS*R1 - T1**3*SS - 4*T1
     .**2*SS*R1 - T1*SS**2*R1 - 5*T1*SS*R1**2 + 2*T1*SS*R1
     . - SS**2*R1**2 - 2*SS*R1**3 + 4*SS*R1**2) + 4*D0ACB*(
     . - U1**2*T1 - U1*T1**2 - 2*U1*T1*R1 + U1*SS*R1 + T1*SS
     .*R1 + 2*SS*R1**2) + 2*D0ABC*( - U1**3*SS - 2*U1**2*T1
     . - 4*U1**2*SS*R1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*T1*
     .R1 - U1*SS**2*R1 - 5*U1*SS*R1**2 + 2*U1*SS*R1 - T1*SS*
     .R1**2 + 2*T1*SS*R1 - SS**2*R1**2 - 2*SS*R1**3 + 4*SS*
     .R1**2))/(U1*T1 - SS*R1)
 
      HH1 = (16*C0AB*SS + 2*C0AC*T1*(U1 + T1 + SS + 2*R1 - 8) + 2*C0AD
     .*( - U1*T1 - U1*SS - T1**2 - 2*T1*SS - 2*T1*R1 + 8*T1 -
     .SS**2 - 2*SS*R1 + 8*SS) + 2*C0BC*U1*(U1 + T1 + SS + 2*
     .R1 - 8) + 2*C0BD*( - U1**2 - U1*T1 - 2*U1*SS - 2*U1*R1
     . + 8*U1 - T1*SS - SS**2 - 2*SS*R1 + 8*SS) + 4*D0BAC*SS
     .*( - U1 - T1 - 2*SS - 2*R1 + 8) + 2*D0ACB*( - U1**2*T1 -
     .U1*T1**2 - U1*T1*SS - 2*U1*T1*R1 + 8*U1*T1 + U1*SS*R1
     . - 2*U1*SS + T1*SS*R1 - 2*T1*SS + SS**2*R1 - 4*SS**2
     . + 2*SS*R1**2 - 12*SS*R1 + 16*SS) + 4*D0ABC*SS*( - U1
     . - T1 - 2*SS - 2*R1 + 8) + 8*SS)/SS
 
      HH2 = (2*C0AB*SS*(U1**2 + 4*U1*R1 - 8*U1 + T1**2 + 4*T1*R1 - 8*
     .T1 + 2*SS*R1 + 4*R1**2 - 16*R1) + 2*C0AC*T1*(U1*R1 +
     .T1**2 + 3*T1*R1 - 8*T1 + SS*R1 + 2*R1**2 - 8*R1) + 2*
     .C0CD*( - U1**3 - U1**2*T1 - 2*U1**2*R1 + 8*U1**2 - U1*T1
     .**2 + 2*U1*SS*R1 - T1**3 - 2*T1**2*R1 + 8*T1**2 + 2
     .*T1*SS*R1 + 4*SS*R1**2 - 16*SS*R1) + 2*C0AD*( - U1
     .**2*T1 - U1**2*SS - 3*U1*T1*R1 + 8*U1*T1 - 3*U1*SS*R1
     . + 8*U1*SS - T1**2*R1 - 2*T1*SS*R1 - 2*T1*R1**2 + 8*T1
     .*R1 - SS**2*R1 - 2*SS*R1**2 + 8*SS*R1) + 2*C0BC*U1*(U1
     .**2 + 3*U1*R1 - 8*U1 + T1*R1 + SS*R1 + 2*R1**2 - 8*R1))
      HH2 = (HH2
     . + 2*C0BD*( - U1**2*R1 - U1*T1**2 - 3*U1*T1*R1 + 8*U1*T1
     . - 2*U1*SS*R1 - 2*U1*R1**2 + 8*U1*R1 - T1**2*SS - 3*T1
     .*SS*R1 + 8*T1*SS - SS**2*R1 - 2*SS*R1**2 + 8*SS*R1) +
     .2*D0BAC*( - 2*U1**2*T1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*
     .T1*R1 + 16*U1*T1 - U1*SS*R1**2 + 2*U1*SS*R1 - T1
     .**3*SS - 4*T1**2*SS*R1 + 8*T1**2*SS - T1*SS**2*
     .R1 - 5*T1*SS*R1**2 + 18*T1*SS*R1 - SS**2*R1**2
     . - 2*SS*R1**3 + 12*SS*R1**2 - 16*SS*R1))
      HH2 = ( HH2 + 4*
     .D0ACB*( - U1**2*T1 - U1*T1**2 - 2*U1*T1*R1 + 8*U1*T1 + U1
     .*SS*R1 + T1*SS*R1 + 2*SS*R1**2 - 8*SS*R1) + 2*
     .D0ABC*( - U1**3*SS - 2*U1**2*T1 - 4*U1**2*SS*R1 + 8*U1**2
     .*SS - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*T1*R1 + 16*
     .U1*T1 - U1*SS**2*R1 - 5*U1*SS*R1**2 + 18*U1*SS*R1
     . - T1*SS*R1**2 + 2*T1*SS*R1 - SS**2*R1**2 - 2*SS*
     .R1**3 + 12*SS*R1**2 - 16*SS*R1))/(U1*T1 - SS*R1)
 
      AH1 = (2*C0AC*T1*(U1 + T1 + SS) + 2*C0AD*( - U1*T1 - U1*SS - T1
     .**2 - 2*T1*SS - SS**2) + 2*C0BC*U1*(U1 + T1 + SS) + 2*
     .C0BD*( - U1**2 - U1*T1 - 2*U1*SS - T1*SS - SS**2) + 4*
     .D0BAC*SS*( - U1 - T1 - 2*SS) + 2*D0ACB*( - U1**2*T1 - U1*
     .T1**2 - U1*T1*SS + U1*SS*R1 - 2*U1*SS + T1*SS*R1 - 2*
     .T1*SS + SS**2*R1 - 4*SS**2) + 4*D0ABC*SS*( - U1 - T1
     . - 2*SS))/SS
 
      AH2 = (2*C0AB*SS*( - U1**2 - 2*U1*R1 + T1**2 + 2*T1*R1) + 2*C0AC
     .*T1*( - U1*R1 + T1**2 + T1*R1 - SS*R1) + 2*C0CD*(U1**3 +
     .U1**2*T1 - U1*T1**2 - 4*U1*SS*R1 - T1**3 + 4*T1*SS*R1)
     . + 2*C0AD*(U1**2*T1 + U1**2*SS + U1*T1*R1 + U1*SS*R1 - T1
     .**2*R1 - 2*T1*SS*R1 - SS**2*R1) + 2*C0BC*U1*( - U1**2
     . - U1*R1 + T1*R1 + SS*R1) + 2*C0BD*(U1**2*R1 - U1*T1**
     .2 - U1*T1*R1 + 2*U1*SS*R1 - T1**2*SS - T1*SS*R1 + SS**
     .2*R1) + 2*D0BAC*(2*U1**2*T1 - 2*U1*T1**2 + U1*T1*SS*R1
     . + U1*SS*R1**2 - 2*U1*SS*R1 - T1**3*SS - 2*T1**2*SS*R1
     . + T1*SS**2*R1 - T1*SS*R1**2 + 2*T1*SS*R1 + SS**2*R1**2))
      AH2 = (AH2
     . + 4*D0ACB*(U1**2*T1 - U1*T1**2 - U1*SS*R1 + T1*SS*
     .R1) + 2*D0ABC*(U1**3*SS + 2*U1**2*T1 + 2*U1**2*SS*R1
     . - 2*U1*T1**2 - U1*T1*SS*R1 - U1*SS**2*R1 + U1*SS*R1**
     .2 - 2*U1*SS*R1 - T1*SS*R1**2 + 2*T1*SS*R1 - SS**2*R1**
     .2))/(U1*T1 - SS*R1)

      A1 = A1*AMQ
      H1 = H1*AMQ
      Z1 = -Z1*AMQ**2
      A2 = A2*AMQ
      H2 = H2*AMQ
      Z2 = -Z2*AMQ**2
c      WRITE(*,*) H1 
c      WRITE(*,*) HH1
c      WRITE (*,*) HH2
      RETURN
      END
