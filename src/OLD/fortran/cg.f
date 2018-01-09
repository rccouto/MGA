      SUBROUTINE ZXCGR  (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFN,IER
      DOUBLE PRECISION   ACC,DFPRED,X(N),G(N),F,W(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG,
     1                   IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NCALLS,NFBEG,
     2                   NFOPT
      DOUBLE PRECISION   BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA,
     1                   GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH,
     2                   STMIN,SUM,WORK
      DATA               MAXLIN/5/,MXFCON/2/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  THE WORKING SPACE ARRAY IS SPLIT
C                                    INTO SIX VECTORS OF LENGTH N. THE
C                                    FIRST PART IS USED FOR THE SEARCH
C                                    DIRECTION OF AN ITERATION. THE
C                                    SECOND AND THIRD PARTS CONTAIN THE
C                                    INFORMATION THAT IS REQUIRED BY
C                                    THE CONJUGACY CONDITIONS OF THE
C                                    RESTART PROCEDURE. THE FOURTH PART
C                                    CONTAINS THE GRADIENT AT THE START
C                                    OF AN ITERATION. THE FIFTH PART
C                                    CONTAINS THE PARAMETERS THAT GIVE
C                                    THE LEAST CALCULATED VALUE OF F.
C                                    THE SIXTH PART CONTAINS THE
C                                    GRADIENT VECTOR WHERE F IS LEAST.
      IRSDX = N
      IRSDG = IRSDX+N
      IGINIT = IRSDG+N
      IXOPT = IGINIT+N
      IGOPT = IXOPT+N
C                                  SET SOME PARAMETERS TO BEGIN THE
C                                    CALCULATION. ITERC AND
C                                    NCALLS COUNT THE NUMBER OF
C                                    ITERATIONS AND CALLS OF FUNCT.
C                                    ITERFM IS THE NUMBER OF THE MOST
C                                    RECENT ITERATION THAT DECREASES F.
      ITERC = 0
      NCALLS = 0
      ITERFM = ITERC
C                                  CALL SUBROUTINE FUNCT. LET THE
C                                    INITIAL SEARCH DIRECTION BE MINUS
C                                    THE GRADIENT VECTOR. USUALLY THE
C                                    PARAMETER ITERRS GIVES THE
C                                    ITERATION NUMBER OF THE MOST
C                                    RECENT RESTART, BUT IT IS SET TO
C                                    ZERO WHEN THE STEEPEST DESCENT
C                                    DIRECTION IS USED.
    5 NCALLS = NCALLS+1
      CALL FUNCT (N,X,F,G)
      IF (NCALLS.GE.2) GO TO 20
   10 DO 15 I=1,N
   15 W(I) = -G(I)
      ITERRS = 0
      IF (ITERC.GT.0) GO TO 80
C                                  SET SUM TO G SQUARED. GMIN AND GNEW
C                                    ARE THE OLD AND THE NEW
C                                    DIRECTIONAL DERIVATIVES ALONG THE
C                                    CURRENT SEARCH DIRECTION. LET FCH
C                                    BE THE DIFFERENCE BETWEEN F AND
C                                    THE PREVIOUS BEST VALUE OF THE
C                                    OBJECTIVE FUNCTION.
   20 GNEW = 0.0D0
      SUM = 0.0D0
      DO 25 I=1,N
         GNEW = GNEW+W(I)*G(I)
   25 SUM = SUM+G(I)**2
      IF (NCALLS.EQ.1) GO TO 35
      FCH = F-FMIN
C                                  STORE THE VALUES OF X, F AND G, IF
C                                    THEY ARE THE BEST THAT HAVE BEEN
C                                    CALCULATED SO FAR, AND NOTE G
C                                    SQUARED AND THE VALUE OF NCALLS.
C                                    TEST FOR CONVERGENCE.
      IF (FCH) 35,30,50
   30 IF (GNEW/GMIN.LT.-1.0D0) GO TO 45
   35 FMIN = F
      GSQRD = SUM
      NFOPT = NCALLS
      DO 40 I=1,N
         W(IXOPT+I) = X(I)
   40 W(IGOPT+I) = G(I)
   45 IF (SUM.LE.ACC) GO TO 9005
C                                  TEST IF THE VALUE OF MAXFN ALLOWS
C                                    ANOTHER CALL OF FUNCT.
   50 IF (NCALLS.NE.MAXFN) GO TO 55
      IER = 131
      GO TO 9000
   55 IF (NCALLS.GT.1) GO TO 100
C                                  SET DFPR TO THE ESTIMATE OF THE
C                                    REDUCTION IN F GIVEN IN THE
C                                    ARGUMENT LIST, IN ORDER THAT THE
C                                    INITIAL CHANGE TO THE PARAMETERS
C                                    IS OF A SUITABLE SIZE. THE VALUE
C                                    OF STMIN IS USUALLY THE
C                                    STEP-LENGTH OF THE MOST RECENT
C                                    LINE SEARCH THAT GIVES THE LEAST
C                                    CALCULATED VALUE OF F.
      DFPR = DFPRED
      STMIN = DFPRED/GSQRD
C                                  BEGIN THE ITERATION
   80 ITERC = ITERC+1
C                                  STORE THE INITIAL FUNCTION VALUE AND
C                                    GRADIENT, CALCULATE THE INITIAL
C                                    DIRECTIONAL DERIVATIVE, AND BRANCH
C                                    IF ITS VALUE IS NOT NEGATIVE. SET
C                                    SBOUND TO MINUS ONE TO INDICATE
C                                    THAT A BOUND ON THE STEP IS NOT
C                                    KNOWN YET, AND SET NFBEG TO THE
C                                    CURRENT VALUE OF NCALLS. THE
C                                    PARAMETER IRETRY SHOWS THE NUMBER
C                                    OF ATTEMPTS AT SATISFYING THE BETA
C                                    CONDITION.
      FINIT = F
      GINIT = 0.0D0
      DO 85 I=1,N
         W(IGINIT+I) = G(I)
   85 GINIT = GINIT+W(I)*G(I)
      IF (GINIT.GE.0.0D0) GO TO 165
      GMIN = GINIT
      SBOUND = -1.0D0
      NFBEG = NCALLS
      IRETRY = -1
C                                  SET STEPCH SO THAT THE INITIAL
C                                    STEP-LENGTH IS CONSISTENT WITH THE
C                                    PREDICTED REDUCTION IN F, SUBJECT
C                                    TO THE CONDITION THAT IT DOES NOT
C                                    EXCEED THE STEP-LENGTH OF THE
C                                    PREVIOUS ITERATION. LET STMIN BE
C                                    THE STEP TO THE LEAST CALCULATED
C                                    VALUE OF F.
      STEPCH = DMIN1(STMIN,DABS(DFPR/GINIT))
      STMIN = 0.0D0
C                                  CALL SUBROUTINE FUNCT AT THE VALUE
C                                    OF X THAT IS DEFINED BY THE NEW
C                                    CHANGE TO THE STEP-LENGTH, AND LET
C                                    THE NEW STEP-LENGTH BE STEP. THE
C                                    VARIABLE WORK IS USED AS WORK
C                                    SPACE.
   90 STEP = STMIN+STEPCH
      WORK = 0.0D0
      DO 95 I=1,N
         X(I) = W(IXOPT+I)+STEPCH*W(I)
   95 WORK = DMAX1(WORK,DABS(X(I)-W(IXOPT+I)))
      IF (WORK.GT.0.0D0) GO TO 5
C                                  TERMINATE THE LINE SEARCH IF STEPCH
C                                    IS EFFECTIVELY ZERO.
      IF (NCALLS.GT.NFBEG+1) GO TO 115
      IF (DABS(GMIN/GINIT)-0.2D0) 170,170,115
C                                  LET SPLN BE THE QUADRATIC SPLINE
C                                    THAT INTERPOLATES THE CALCULATED
C                                    FUNCTION VALUES AND DIRECTIONAL
C                                    DERIVATIVES AT THE POINTS STMIN
C                                    AND STEP OF THE LINE SEARCH, WHERE
C                                    THE KNOT OF THE SPLINE IS AT
C                                    0.5*(STMIN+STEP). REVISE STMIN,
C                                    GMIN AND SBOUND, AND SET DDSPLN TO
C                                    THE SECOND DERIVATIVE OF SPLN AT
C                                    THE NEW STMIN. HOWEVER, IF FCH IS
C                                    ZERO, IT IS ASSUMED THAT THE
C                                    MAXIMUM ACCURACY IS ALMOST
C                                    ACHIEVED, SO DDSPLN IS CALCULATED
C                                    USING ONLY THE CHANGE IN THE
C                                    GRADIENT.
  100 WORK = (FCH+FCH)/STEPCH-GNEW-GMIN
      DDSPLN = (GNEW-GMIN)/STEPCH
      IF (NCALLS.GT.NFOPT) SBOUND = STEP
      IF (NCALLS.GT.NFOPT) GO TO 105
      IF (GMIN*GNEW.LE.0.0D0) SBOUND = STMIN
      STMIN = STEP
      GMIN = GNEW
      STEPCH = -STEPCH
  105 IF (FCH.NE.0.0D0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH
C
C                                  TEST FOR CONVERGENCE OF THE LINE
C                                    SEARCH, BUT FORCE AT LEAST TWO
C                                    STEPS TO BE TAKEN IN ORDER NOT TO
C                                    LOSE QUADRATIC TERMINATION.
      IF (GMIN.EQ.0.0D0) GO TO 170
      IF (NCALLS.LE.NFBEG+1) GO TO 120
      IF (DABS(GMIN/GINIT).LE.0.2D0) GO TO 170
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MAXLIN.
  110 IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120
  115 IER = 129
      GO TO 170
C                                  SET STEPCH TO THE GREATEST CHANGE TO
C                                    THE CURRENT VALUE OF STMIN THAT IS
C                                    ALLOWED BY THE BOUND ON THE LINE
C                                    SEARCH. SET GSPLN TO THE GRADIENT
C                                    OF THE QUADRATIC SPLINE AT
C                                    (STMIN+STEPCH). HENCE CALCULATE
C                                    THE VALUE OF STEPCH THAT MINIMIZES
C                                    THE SPLINE FUNCTION, AND THEN
C                                    OBTAIN THE NEW FUNCTION AND
C                                    GRADIENT VECTOR, FOR THIS VALUE OF
C                                    THE CHANGE TO THE STEP-LENGTH.
  120 STEPCH = 0.5D0*(SBOUND-STMIN)
      IF (SBOUND.LT.-0.5D0) STEPCH = 9.0D0*STMIN
      GSPLN = GMIN+STEPCH*DDSPLN
      IF (GMIN*GSPLN.LT.0.0D0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)
      GO TO 90
C                                  CALCULATE THE VALUE OF BETA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION.
  125 SUM = 0.0D0
      DO 130 I=1,N
  130 SUM = SUM+G(I)*W(IGINIT+I)
      BETA = (GSQRD-SUM)/(GMIN-GINIT)
C                                  TEST THAT THE NEW SEARCH DIRECTION
C                                    CAN BE MADE DOWNHILL. IF IT
C                                    CANNOT, THEN MAKE ONE ATTEMPT TO
C                                    IMPROVE THE ACCURACY OF THE LINE
C                                    SEARCH.
      IF (DABS(BETA*GMIN).LE.0.2D0*GSQRD) GO TO 135
      IRETRY = IRETRY+1
      IF (IRETRY.LE.0) GO TO 110
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MXFCON.
C                                    SET DFPR TO THE PREDICTED
C                                    REDUCTION IN F ON THE NEXT
C                                    ITERATION.
  135 IF (F.LT.FINIT) ITERFM = ITERC
      IF (ITERC.LT.ITERFM+MXFCON) GO TO 140
      IER = 132
      GO TO 9000
  140 DFPR = STMIN*GINIT
C                                  BRANCH IF A RESTART PROCEDURE IS
C                                    REQUIRED DUE TO THE ITERATION
C                                    NUMBER OR DUE TO THE SCALAR
C                                    PRODUCT OF CONSECUTIVE GRADIENTS.
      IF (IRETRY.GT.0) GO TO 10
      IF (ITERRS.EQ.0) GO TO 155
      IF (ITERC-ITERRS.GE.N) GO TO 155
      IF (DABS(SUM).GE.0.2D0*GSQRD) GO TO 155
C                                  CALCULATE THE VALUE OF GAMA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION, AND SET SUM TO A SCALAR
C                                    PRODUCT FOR THE TEST BELOW. THE
C                                    VALUE OF GAMDEN IS SET BY THE
C                                    RESTART PROCEDURE.
      GAMA = 0.0D0
      SUM = 0.0D0
      DO 145 I=1,N
         GAMA = GAMA+G(I)*W(IRSDG+I)
  145 SUM = SUM+G(I)*W(IRSDX+I)
      GAMA = GAMA/GAMDEN
C                                  RESTART IF THE NEW SEARCH DIRECTION
C                                    IS NOT SUFFICIENTLY DOWNHILL.
C
      IF (DABS(BETA*GMIN+GAMA*SUM).GE.0.2D0*GSQRD) GO TO 155
C
C                                  CALCULATE THE NEW SEARCH DIRECTION.
      DO 150 I=1,N
  150 W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)
      GO TO 80
C                                  APPLY THE RESTART PROCEDURE.
  155 GAMDEN = GMIN-GINIT
      DO 160 I=1,N
         W(IRSDX+I) = W(I)
         W(IRSDG+I) = G(I)-W(IGINIT+I)
  160 W(I) = -G(I)+BETA*W(I)
      ITERRS = ITERC
      GO TO 80
C                                  SET IER TO INDICATE THAT THE SEARCH
C                                    DIRECTION IS UPHILL.
  165 IER = 130
C                                  ENSURE THAT F, X AND G ARE OPTIMAL.
  170 IF (NCALLS.EQ.NFOPT) GO TO 180
      F = FMIN
      DO 175 I=1,N
         X(I) = W(IXOPT+I)
  175 G(I) = W(IGOPT+I)
  180 IF (IER.EQ.0) GO TO 125
 9000 CONTINUE
      CALL UERTST (IER,6HZXCGR )
 9005 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  UERTST
C  ZXMJN
C  USPKD
C  UGETIO
C
      SUBROUTINE UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER            NAME(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQ,IEQDF,IOUNIT,LEVEL,LEVOLD,NAMEQ(6),
     *                   NAMSET(6),NAMUPK(6),NIN,NMTB
      DATA               NAMSET/1HU,1HE,1HR,1HS,1HE,1HT/
      DATA               NAMEQ/6*1H /
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
C                                  UNPACK NAME INTO NAMUPK
C                                  FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
C                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
C                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
C
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAMUPK
C                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE ZXMJN  (A,N,Z,SIG,W,IR,MK,EPS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR,MK
      DOUBLE PRECISION   A(1),Z(N),SIG,W(N),EPS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,JJ,IJ,JP1,I,II,MM
      DOUBLE PRECISION   ZERO,ONE,FOUR,TI,V,TIM,AL,R,B,GM,Y
      DATA               ZERO/0.0D0/,ONE/1.0D0/,FOUR/4.0D0/
C                                  UPDATE FACTORS GIVEN IN A
C                                    SIG*Z*Z-TRANSPOSE IS ADDED
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.GT.1) GO TO 5
C                                  N .EQ. 1
      A(1) = A(1)+SIG*Z(1)*Z(1)
      IR = 1
      IF (A(1).GT.ZERO) GO TO 9005
      A(1) = ZERO
      IR = 0
      GO TO 9005
C                                  N .GT. 1
    5 IF (SIG.GT.ZERO) GO TO 65
      IF (SIG.EQ.ZERO.OR.IR.EQ.0) GO TO 9005
      TI = ONE/SIG
      JJ = 0
      IF (MK.EQ.0) GO TO 15
C                                  L*W = Z ON INPUT
      DO 10 J=1,N
         JJ = JJ+J
         IF (A(JJ).NE.ZERO) TI = TI+(W(J)*W(J))/A(JJ)
   10 CONTINUE
      GO TO 40
C                                  SOLVE L*W = Z
   15 DO 20 J=1,N
         W(J) = Z(J)
   20 CONTINUE
      DO 35 J=1,N
         JJ = JJ+J
         V = W(J)
         IF (A(JJ).GT.ZERO) GO TO 25
         W(J) = ZERO
         GO TO 35
   25    TI = TI+(V*V)/A(JJ)
         IF (J.EQ.N) GO TO 35
         IJ = JJ
         JP1 = J+1
         DO 30 I=JP1,N
            IJ = IJ+I-1
            W(I) = W(I)-V*A(IJ)
   30    CONTINUE
   35 CONTINUE
C                                  SET TI, TIM AND W
   40 IF (IR.LE.0) GO TO 45
      IF (TI.GT.ZERO) GO TO 50
      IF (MK-1) 65,65,55
   45 TI = ZERO
      IR = -IR-1
      GO TO 55
   50 TI = EPS/SIG
      IF (EPS.EQ.ZERO) IR = IR-1
   55 TIM = TI
      II = JJ
      I = N
      DO 60 J=1,N
         IF (A(II).NE.ZERO) TIM = TI-(W(I)*W(I))/A(II)
         W(I) = TI
         TI = TIM
         II = II-I
         I = I-1
   60 CONTINUE
      MM = 1
      GO TO 70
   65 MM = 0
      TIM = ONE/SIG
   70 JJ = 0
C                                  UPDATE A
      DO 110 J=1,N
         JJ = JJ+J
         IJ = JJ
         JP1 = J+1
C                                  UPDATE A(J,J)
         V = Z(J)
         IF (A(JJ).GT.ZERO) GO TO 85
C                                  A(J,J) .EQ. ZERO
         IF (IR.GT.0.OR.SIG.LT.ZERO.OR.V.EQ.ZERO) GO TO 80
         IR = 1-IR
         A(JJ) = (V*V)/TIM
         IF (J.EQ.N) GO TO 9005
         DO 75 I=JP1,N
            IJ = IJ+I-1
            A(IJ) = Z(I)/V
   75    CONTINUE
         GO TO 9005
   80    TI = TIM
         GO TO 110
C                                  A(J,J) .GT. ZERO
   85    AL = V/A(JJ)
         TI = W(J)
         IF (MM.EQ.0) TI = TIM+V*AL
         R = TI/TIM
         A(JJ) = R*A(JJ)
         IF (R.EQ.ZERO) GO TO 115
         IF (J.EQ.N) GO TO 115
C                                  UPDATE REMAINDER OF COLUMN J
         B = AL/TI
         IF (R.GT.FOUR) GO TO 95
         DO 90 I=JP1,N
            IJ = IJ+I-1
            Z(I) = Z(I)-V*A(IJ)
            A(IJ) = A(IJ)+B*Z(I)
   90    CONTINUE
         GO TO 105
   95    GM = TIM/TI
         DO 100 I=JP1,N
            IJ = IJ+I-1
            Y = A(IJ)
            A(IJ) = B*Z(I)+Y*GM
            Z(I) = Z(I)-V*Y
  100    CONTINUE
  105    TIM = TI
  110 CONTINUE
  115 IF (IR.LT.0) IR = -IR
 9005 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C
      LOGICAL*1          UNPAKD(1),PACKED(1),LBYTE,LBLANK
      INTEGER*2          IBYTE,IBLANK
      EQUIVALENCE (LBYTE,IBYTE)
      DATA               LBLANK /1H /
      DATA               IBYTE /1H /
      DATA               IBLANK /1H /
C                                  INITIALIZE NCHMTB
      NCHMTB = 0
C                                  RETURN IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) RETURN
C                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
      NWORDS = NC*4
      J = 1
      DO 110 I = 1,NWORDS,4
      UNPAKD(I) = PACKED(J)
      UNPAKD(I+1) = LBLANK
      UNPAKD(I+2) = LBLANK
      UNPAKD(I+3) = LBLANK
  110 J = J+1
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB
C                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NWORDS,4
         NN = NWORDS - N - 2
         LBYTE = UNPAKD(NN)
         IF(IBYTE .NE. IBLANK) GO TO 210
  200 CONTINUE
  210 NCHMTB = (NN + 3) / 4
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
