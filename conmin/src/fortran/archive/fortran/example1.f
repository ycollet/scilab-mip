      PROGRAM EXAMPL1
      implicit double precision(a-h,o-z)
      DIMENSION S(6),G1(11),G2(11),B(11,11),C(11),MS1(22)
      DIMENSION VLB(6),VUB(6),SCAL(6),DF(6),A(6,11),
     .              ISC(11),IC(11)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,
     .               ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,
     .               NFDG,NSCAL,LINOBJ,ITMAX,ITRM,ICNDIR,IGOTO,NAC,
     .               INFO,INFOG,ITER
      COMMON /VARABLE/ AOBJ,X(6),G(11)
      COMMON /ANDATA/ LOOPCNT
      NAMELIST /CONPAR/ INFOG,INFO,NFDG,IPRINT,NDV,ITMAX,NCON,NSIDE,
     .                  ICNDIR,NSCAL,FDCH,FDCHM,CT,CTMIN,CTLMIN,THETA,
     .                  PHI,DELFUN,DABFUN,LINOBJ,ITRM,X,VLB,VUB,
     .                  N1,N2,N3,N4,N5,ALPHAX,ABOBJ1,CTL,ISC,SCAL
C
C     THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM ONE OF THE CONMIN
C     MANUAL.
C
C
C  INITIALIZE
C
      INFOG=0
      INFO=0
      NFDG=0
      IPRINT=2
      NDV=4
      ITMAX=40
      NCON=3
      NSIDE=0
      ICNDIR=0
      NSCAL=0
      FDCH=0.0
      FDCHM=0.0
      CT=0.0
      CTMIN=0.0
      CTL=0.0
      CTLMIN=0.0
      THETA=0.0
      PHI=0.0
      DELFUN=0.0
      DABFUN=0.0
      LINOBJ=0.0
      ITRM=0
      N1=6
      N2=11
      N3=11
      N4=11
      N5=22
      ALPHAX=0.0
      ABOBJ1=0.0
      CTL=0.0
      DO 5 I=1,NDV
        X(I)=1.0
        VLB(I)=-99999.
        VUB(I)= 99999.
    5 CONTINUE
C
      DO 6 J=1,NCON
        ISC(J)=0
    6 CONTINUE
C
C     READ THE PARAMETERS FOR CONMIN
C
CCC   READ(5,CONPAR)                  USE DEFAULT VALUES
      WRITE(6,CONPAR)
      NLIM=ITMAX*(NDV+5)
C
C     NON-ITERATIVE PART OF ANALYSIS
C
      IGOTO = 0
C
C     ITERATIVE PART OF ANALYSIS
C
      DO 1000 I = 1,NLIM
        LOOPCNT=I
C
C       CALL THE OPTIMIZATION ROUTINE CONMIN
C
        CALL CONMIN(X,VLB,VUB,G,SCAL,DF,A,S,G1,G2,B,C,ISC,IC,MS1,N1,N2,
     .              N3,N4,N5)
C
        IF(IGOTO.EQ.0) LOOPCNT=-999
C
C       ANALYSIS MODULE
C
        CALL ANALYS
        OBJ=AOBJ
        IF (IGOTO.EQ.0) GO TO 1100
 1000 CONTINUE
C
C
 1100 CONTINUE
      STOP
      END



      SUBROUTINE ANALYS
      implicit double precision(a-h,o-z)
      COMMON /VARABLE/ AOBJ,X(6),G(11)
C
C   ROUTINE TO CALCULATE OBJECTIVE FUNCTION AND
C   CONSTRAINTS
C
C
C  OBJECTIVE FUNCTION
C
      AOBJ = X(1)**2 - 5.*X(1) + X(2)**2 - 5.*X(2) + 2.*X(3)**2
     .       - 21.*X(3) + X(4)**2 + 7.0*X(4) + 50.
C
C
C   CONSTRAINT VALUES
C
      G(1) = X(1)**2 + X(1) + X(2)**2 - X(2) + X(3)**2 + X(3)
     .       + X(4)**2 - X(4) - 8.0
C
      G(2) = X(1)**2 - X(1) + 2. * X(2)**2 + X(3)**2 + 2.*X(4)**2
     .       - X(4) - 10.0
C
      G(3) = 2.*X(1)**2 + 2.*X(1) + X(2)**2 - X(2) + X(3)**2 - X(4) -5.0
C
      RETURN
      END
