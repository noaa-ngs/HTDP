      PROGRAM HTDP

**************************************************************
*  NAME:       HTDP (Horizontal Time-Dependent Positioning)
*
*  WRITTEN BY: Richard A. Snay & Chris Pearson
*
*  PURPOSE:    Transform coordinates across time
*              and between reference frames
*
****************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      character    HTDP_version*10
      CHARACTER    OPTION*1
      character    cont*5

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /VERSION/ HTDP_version

C  You must change HTDP version here if necessary

      HTDP_version = 'v3.2.8'

*** Introduce variables for file id's

      LUIN = 5
*       interactive input
      LUOUT = 6
*       interactive output
      I1 = 11
*       input of velocity grid in GETVEL      
*       input of earthquake parameters in GETEQ
*       input of blue-book BFILE in DLACE, VELOC, UPDATE and TRFPOS
*       input of blue-book GFILE in UPDATE
      I2 = 12
*       output of predicted displacements in DLACE
*       output of predicted velocities in VELOC
*       output of updated blue-book BFILE in UPDATE
*       output of updated blue-book GFILE in UPDATE
*       output of updated coordinates in UPDATE
      I3 = 13
*       output of point-velocity records in VELOC
      I4 = 14
*       storage of reference latitude & longitude & region
      I5 = 15
*       storage of earthquake parameters
C     OPEN(I5, FILE='TEMPQK' , ACCESS = 'DIRECT', RECL = 72,
C    1      FORM='UNFORMATTED')
      I6 = 16
*       output of transformed blue-book BFILE in TRFPOS

*** Obtain parameters defining crustal motion model
      CALL MODEL
       
*** Initialize transformation parameters between reference frames
      CALL SETTP 

*** Initialize conversion table between reference frame identifiers
      CALL SETRF

      WRITE(LUOUT,5) HTDP_version
    5 FORMAT(
     1 '********************************************************'/
     1 '*   HTDP (Horizontal Time-Dependent Positioning)       *'/
     1 '*   SOFTWARE VERSION ',a10                               / 
     1 '*                                                      *'/)
      WRITE(LUOUT,501) 
  501 FORMAT(
     1 '*   AUTHORS: Richard Snay, Chris Pearson & Jarir Saleh *'/
     1 '*            Email: ngs.cors.htdp@noaa.gov             *'/
     1 '*                                                      *'/
     1 '********************************************************'/)
      WRITE(LUOUT,10)
   10 FORMAT( 
     1 ' This software incorporates numerical models that',/
     3 ' characterize continuous crustal motion as well as ',/
     3 ' the episodic motion associated with earthquakes.'/)
      WRITE(LUOUT,11)
   11 FORMAT(
     5 ' The User Guide contains additional information and a set'/
     5 ' of exercises to familiarize users with the software.'//
     5 ' Hit ENTER or RETURN to continue.  ')
      read(luin, '(a5)',err=51,iostat=ios) cont
      if (ios /= 0) goto 51

   25 WRITE(LUOUT,26)
   26 FORMAT(' ***************************************'/
     1 ' MAIN MENU:',/
     6 '    0... Exit software.',/
     7 '    1... Estimate displacements between two dates.'/
     8 '    2... Estimate velocities.'/
     9 '    3... Update positions and/or observations '
     &           ,'to a specified date.'/
     & '    4... Transform positions between reference frames.  '/
     & '    5... Transform velocities between reference frames. ')
   30 READ(LUIN,35,err=52,iostat=ios) OPTION
      if (ios /= 0) goto 52
   35 FORMAT(A1)
      IF(OPTION .EQ. '0') THEN
        GO TO 50                   
      ELSEIF(OPTION .EQ. '1') THEN
        CALL DPLACE
      ELSEIF(OPTION .EQ. '2') THEN
        CALL VELOC
      ELSEIF(OPTION .EQ. '3') THEN
        CALL UPDATE(HTDP_version)
      elseif(option .eq. '4') then
        call TRFPOS
      elseif(option .eq. '5') then
        call TRFVEL
      ELSE
        WRITE(LUOUT,40)
   40   FORMAT(' Improper entry--select again  ')
        GO TO 30
      ENDIF
      GO TO 25
   50 CONTINUE
      stop
C     CLOSE(I5, STATUS = 'DELETE')

   51 write (*,'(/)')
      write (*,*) 'You did not hit enter      : ios = ',ios
      write (*,*) "ABNORMAL TERMINATION"
      STOP

   52 write (*,'(/)')
      write (*,*) 'Problem with reading OPTION: ios = ',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      STOP

      END
*****************************************************************************
      SUBROUTINE MODEL

*** Obtain parameters defining crustal motion model

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF

      A = 6.378137D06
      F = 1.D0 / 298.257222101D0
      E2 = 0.6694380022903146D-2
      AF = A / (1.D0 -F)
      EPS = F*(2.D0 - F) / ((1.D0 -F)**2)
      PI = 4.D0 * DATAN(1.D0)
      RHOSEC = (180.D0 * 3600.D0) / PI
      TWOPI = PI + PI

C*** Set default reference epoch to Jan. 1, 2010
      IYRREF = 2010
      IMOREF = 1
      IDYREF = 1
      CALL IYMDMJ (IYRREF, IMOREF, IDYREF, MJD)
      ITREF = MJD * 24 * 60

      CALL GETBDY

      RETURN
      END

*******************************************************************
      SUBROUTINE GETBDY

*** Obtain coordinates for vertices that form the polygons
*** that correspond to the boundaries for the regions.       
*** Region 1 is the San Andreas fault in central California         
*** Region 2 is southern California
*** Region 3 is Northern California
*** Region 4 is the Pacific Noerthwest
*** Region 5 is western CONUS
*** Region 6 is CONUS
*** Region 7 is St. Elias, Alaska
*** Region 8 is south-central Alaska
*** Region 9 is southeast Alaska
*** Region 10 is All Mainland Alaska
*** Region 11 is the North American plate 
*** Region 12 is the Caribbean plate
*** Region 13 is the Pacific plate
*** Region 14 is the Juan de Fuca plate
*** Region 15 is the Cocos plate
*** Region 16 is the Mariana plate
*** REGION 17 is the Philippine Sea plate

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 17)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /BNDRY/ X(4000), Y(4000), NPOINT(30)
 
      IEND = NPOINT(NMREGN + 1) - 1  
      DO 10 J = 1, IEND              
        X(J) = (X(J) * 3600.D0)/RHOSEC
        Y(J) = (Y(J) * 3600.D0)/RHOSEC
   10 CONTINUE
      RETURN
      END
*******************************************************************
      SUBROUTINE GETREG(X0,YKEEP,JREGN)

*** Determine the region containing a given point.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 17)
      COMMON /BNDRY/ X(4000), Y(4000), NPOINT(30)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

      Y0 = TWOPI - YKEEP
      IF (Y0 .lt. 0.d0) Y0 = Y0 + TWOPI
      IR = 0     
    1 IR = IR + 1
      IF(IR .GT. NMREGN) THEN
        JREGN = 0
        RETURN
      ENDIF
      IBEGIN = NPOINT(IR)
      NUMVER = NPOINT(IR + 1) - IBEGIN
      CALL POLYIN(X0,Y0,X(IBEGIN),Y(IBEGIN), NUMVER, NTEST)
      IF(NTEST .EQ. 0) GO TO 1
      JREGN = IR       

      RETURN
      END
********************************************************
      SUBROUTINE  POLYIN (X0,Y0,X,Y,N,NPC)
C     SUBROUTINE TO DETERMINE IF A POINT AT (X0,Y0) IS INSIDE OR
C     OUTSIDE OF A CLOSED FIGURE DESCRIBED BY A SEQUENCE OF CONNECTED
C     STRAIGHT LINE SEGMENTS WITH VERTICES AT X, Y.
C
C     INPUT -
C         X0, Y0    COORDINATES OF A POINT TO BE TESTED
C                    Y0 corresponds to longitude and must be a number
C                    between 0.0 and 2*PI
C         X, Y      ARRAYS CONTAINING THE VERTICES, IN ORDER, OF A
C                   CLOSED FIGURE DESCRIBED BY STRAIGHT LINE SEGMNENTS.
C                   FOR EACH 'I', THE STRAIGHT LINE FROM (XI),Y(I)) TO
C                   TO (X(I+1),Y(I+1)), IS AN EDGE OF THE FIGURE.
C         N         DIMENSION OF X AND Y, NUMBER OF VERTICES, AND NUMBER
C                   OF STRAIGHT LINE SEGMENTS IN FIGURE.
C     OUTPUT -
C         NPC       NPC=0 WHEN X0,Y0 IS OUTSIDE OF FIGURE DESCRIBED
C                   BY X,Y
C                   NPC=1 WHEN X0,Y0 IS INSIDE FIGURE
C                   NPC=2 WHEN X0,Y0 IS ON BORDER OF FIGURE
C     METHOD -
C     A COUNT IS MADE OF THE NUMBER OF TIMES THE LINE FROM (X0,Y0) TO
C     (X0,+ INFINITY) CROSSES THE BORDER OF THE FIGURE. IF THE COUNT
C     IS ODD, THE POINT IS INSIDE; IF THE COUNT IS EVEN THE POINT
C     IS OUTSIDE.
C     LIMITATIONS -
C     NONE. THE PROGRAM LOGIC IS VALID FOR ALL CLOSED FIGURES,
C     NO MATTER HOW COMPLEX.
C     ACCURACY -
C     MAINTAINS FULL ACCURACY OF INPUT COORDINATES.
C
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N)
      DATA I6/6/
      IS=0
      NPC=0
C
C     FIND STARTING POINT WHERE X(I).NE.X0
      IP=0
   10 IP=IP+1
      IF ((X(IP)-X0) == 0) goto 12 
      IF ((X(IP)-X0)  < 0) goto 15 
      IF ((X(IP)-X0)  > 0) goto 16 
c     IF (X(IP)-X0) 15,12,16
   12 IF(IP.LE.N) GO TO 10
      WRITE(I6,6001)
 6001 FORMAT('0  POLYGON INPUT ERROR - ALL POINTS ON LINE X = X0')
      STOP
   15 IL=-1
      GO TO 20
   16 IL=1
   20 XL=X(IP)
      YL=Y(IP)
C
C     SET UP SEARCH LOOP
C
      IP1=IP+1
      IPN=IP+N
      DO 100 II=IP1,IPN
      I=II
      IF(I.GT.N) I=I-N
      IF(IL == 0) goto 50
      IF(IL  < 0) goto 30
      IF(IL  > 0) goto 40
c     IF(IL) 30,50,40

   30 IF(X(I)-X0 == 0) goto 32
      IF(X(I)-X0  < 0) goto 90
      IF(X(I)-X0  > 0) goto 34
c  30 IF(X(I)-X0) 90,32,34
   32 IS=-1
      GO TO 60
   34 IL=1
      GO TO 80
   40 IF(X(I)-X0 == 0) goto 44
      IF(X(I)-X0  < 0) goto 42
      IF(X(I)-X0  > 0) goto 90
c  40 IF(X(I)-X0) 42,44,90
   42 IL=-1
      GO TO 80
   44 IS=1
      GO TO 60
   50 IF(X(I)-X0 == 0) goto 55
      IF(X(I)-X0  < 0) goto 52
      IF(X(I)-X0  > 0) goto 54
c  50 IF(X(I)-X0) 52,55,54
   52 IL=-1
      IF(IS == 0) goto 140
      IF(IS  < 0) goto 90
      IF(IS  > 0) goto 80
c     IF(IS) 90,140,80
   54 IL=1
      IF(IS == 0) goto 140
      IF(IS  < 0) goto 80
      IF(IS  > 0) goto 90
c     IF(IS) 80,140,90

   55 IF(Y(I)-Y0 == 0) goto 120
      IF(Y(I)-Y0  < 0) goto 57
      IF(Y(I)-Y0  > 0) goto 58
c  55 IF(Y(I)-Y0) 57,120,58

   57 IF(YL-Y0 == 0) goto 120
      IF(YL-Y0  < 0) goto 90
      IF(YL-Y0  > 0) goto 120
c  57 IF(YL-Y0) 90,120,120

   58 IF(YL-Y0 == 0) goto 120
      IF(YL-Y0  < 0) goto 120
      IF(YL-Y0  > 0) goto  90
c  58 IF(YL-Y0) 120,120,90
C
   60 IL=0
      IF(Y(I)-Y0 == 0) goto 120
      IF(Y(I)-Y0  < 0) goto 90
      IF(Y(I)-Y0  > 0) goto 90
c     IF(Y(I)-Y0) 90,120,90

   80 IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL) == 0) goto 120
      IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)  < 0) goto 90
      IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)  > 0) goto 85
c  80 IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)) 90,120,85
   85 NPC=NPC+1
   90 XL=X(I)
      YL=Y(I)
  100 CONTINUE
      NPC=MOD(NPC,2)
      RETURN
  120 NPC=2
      RETURN
 140  WRITE(I6,6002)
 6002 FORMAT('0  POLYGON LOGIC ERROR - PROGRAM SHOULD NOT REACH THIS',
     .              ' POINT')
      RETURN
      END
*****************************************************************
      SUBROUTINE RADR8T (YLAT,VN,VE,VNR,VER)

C Convert horizontal velocities from mm/yr to rad/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CALL RADII (YLAT,RADMER,RADPAR)
c     write (*,*) "From RADR8T",RADMER,RADPAR
      VNR = VN / (1000.D0 * RADMER)
      VER = VE / (1000.D0 * RADPAR)
c     write (*,*) "From RADR8T",VNR,VER

      RETURN
      END
*****************************************************************
      SUBROUTINE COMVEL(YLAT,YLON,JREGN,VN,VE,VU)
C
C Compute the NAD_83(CORS96) velocity at a point in mm/yr    !Not anymore since 09/12/2014
                                                             !Now the velocity refer to ITRF2008
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 10)
      parameter (NMREGN = 17)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1,I2,I3,I4,I5,I6
      COMMON /VGRID/ B(210000)
      DIMENSION WEI(2,2), VEL(2,2,3)

c     WRITE (6, 1001) JREGN
c1001 FORMAT( 'JREGN = ', I6)

      IF(JREGN .GT. NUMGRD .AND. JREGN .LE. NMREGN) THEN
*** Use tectonic plate model to compute velocity relative
***    to ITRF2008
        IPLATE = JREGN - NUMGRD
        ELON = - YLON
        HT = 0.D0
        CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
        CALL PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
        VX = VX * 1000.D0
        VY = VY * 1000.D0
        VZ = VZ * 1000.D0
*** Convert ITRF2008 velocity to NAD_83(CORS96) velocity
c       CALL VTRANF(X, Y, Z, VX, VY, VZ, 15, 1)          !No Do not. Leave the velocity in ITRF2008
        CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)

      ELSEIF(JREGN .GE. 1 .AND. JREGN .LE. NUMGRD) THEN

C*** Get indices for the lower left hand corner of the grid
C*** and get the weights for the four corners
        CALL GRDWEI (YLON, YLAT, JREGN, I, J, WEI)
c       write (*,*) 'HHHHHHHHHHHHHHHHHHHHHHHHH',JREGN,WEI
c       write (*,*) 'HHHHHHHHHHHHHHHHHHHHHHHHH',i,j

C*** Get the velocity vectors at the four corners
        CALL GRDVEC (JREGN, I, J, VEL, B)

        VN = WEI(1,1) * VEL(1,1,1) + WEI(1,2) * VEL(1,2,1)
     *     + WEI(2,1) * VEL(2,1,1) + WEI(2,2) * VEL(2,2,1)

        VE = WEI(1,1) * VEL(1,1,2) + WEI(1,2) * VEL(1,2,2)
     *     + WEI(2,1) * VEL(2,1,2) + WEI(2,2) * VEL(2,2,2)
  
        VU = WEI(1,1) * VEL(1,1,3) + WEI(1,2) * VEL(1,2,3)
     *     + WEI(2,1) * VEL(2,1,3) + WEI(2,2) * VEL(2,2,3)

c       write (*,*) 'From COMVEL   ',VN,VE,VU
c       write (*,*) 'From COMVEL   ',I,J          

C*** If the point in one of the four Alaskan regions,
C*** then set its vertical velocity to 0.0
        IF(JREGN .GE. 7 .AND. JREGN .LE. 10) THEN
           VU = 0.D0
        ENDIF

C*** If the point is in one of the first ten regions, then
c*** the velocity grids contain the ITRF2008 velocity.
c*** Hence, the following code transforms this ITRF2008 velocity
c*** to the corresponding NAD 83 (CORS96) velocity. 
C
c       IF(JREGN .LE. NUMGRD) THEN                            !Starting09/12/2014, the velocities
c          ELON = - YLON                                      !coming out of this routine are in ITRF2008
c          HT = 0.D0
c          CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
c          CALL TOVXYZ(YLAT, ELON, VN, VE, VU, VX, VY, VZ)
c          CALL VTRANF( X, Y, Z, VX, VY, VZ, 15, 1)
c          CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)
c       ENDIF

      ELSE
        WRITE(LUOUT,100) JREGN
  100   FORMAT(' Improper region identifier ',I4,'in COMVEL.')
        STOP
      ENDIF
      RETURN
      END
****************************************
      SUBROUTINE PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
*** Compute the ITRF2008 velocity at point on plate = IPLATE
***    with coordinates X, Y, Z (in meters)
***    The resulting velocities--VX, VY, and VZ--will be in meters/yr
***    References 
***     Altamimi et al. 2012 = JGR (Paper on ITRF2008 plate motion)
***     DeMets et al. 2010 = Geophysical Journal Int'l, vol 181, 
***     Snay 2003 = SALIS, Vol 63, No 1 (Paper on Frames for Pacific)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION WX(7), WY(7), WZ(7)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

*** IPLATE = 1 --> North America (from Altamimi et al. 2012)
***          2 --> Caribbean (from Altamimi et al. 2012)
***          3 --> Pacific (from Altamimi et al. 2012)
***          4 --> Juan de Fuca (from DeMets et al. 2010)
***          5 --> Cocos (from DeMets et al. 2010)
***          6 --> Mariana (from Snay, 2003)
***          7 --> Philippine Sea (from DeMets et al. 2010)
      DATA WX /0.170D-9, 0.238D-9, -1.993D-9, 
     1         6.626D-9, -10.390D-9, -.097D-9, -0.841D-9/
      DATA WY /-3.209D-9, -5.275D-9, 5.023D-9,
     1         11.708D-9, -14.954D-9,  .509D-9, 3.989D-9/
      DATA WZ /-0.485D-9, 3.219D-9,-10.501D-9, 
     1        -10.615D-9,  9.148D-9,-1.682D-9, -10.626D-9/

      IF (IPLATE .LE. 0 .OR. IPLATE .GT. 7) THEN
          WRITE (LUOUT, 1) IPLATE
    1     FORMAT(' Improper plate ID in PLATVL = ', I6)
          STOP
      ENDIF

      VX = -WZ(IPLATE) * Y + WY(IPLATE) * Z
      VY =  WZ(IPLATE) * X - WX(IPLATE) * Z
      VZ = -WY(IPLATE) * X + WX(IPLATE) * Y

*** The parameters--WX, WY, and WZ--refer to ITRF2000
*** for the Mariana Plate (Snay, 2003). Hence,
*** for this plate, VX, VY, and VZ, correspond to ITRF2000.
*** The following code converts these to ITRF2008 velocities for
*** this plate.
      IF (IPLATE .EQ. 6) THEN
         VX = VX*1000.d0
         VY = VY*1000.d0
         VZ = VZ*1000.d0
         CALL VTRANF(X, Y, Z, VX, VY, VZ, 11, 15)
         VX = VX/1000.d0
         VY = VY/1000.d0
         VZ = VZ/1000.d0
*** The following translations rates are added per Altamimi et al. (2012)
*** for the other six plates
      ELSE
         VX = 0.00041d0 + VX
         VY = 0.00022d0 + VY
         VZ = 0.00041d0 + VZ
      ENDIF

      RETURN
      END
**********************************************************************
      SUBROUTINE PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)

***  Print out a point-velocity (PV) record

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT,I1,I2,I3,I4,I5, I6

      LATS = IDINT(SLAT*100.D0 + 0.5D0)
      LONS = IDINT(SLON*100.D0 + 0.5D0)
      IVN  = IDINT(  VN*100.D0 + 0.5D0)
      IVE  = IDINT(  VE*100.D0 + 0.5D0)
      IVU  = IDINT(  VU*100.D0 + 0.5D0)
      JVN  = 300 
      JVE  = 300 
      JVU  = 500

      WRITE(I3,10) LATD,LATM,LATS,LOND,LONM,LONS,
     1             IVN,JVN,IVE,JVE,IVU,JVU
   10 FORMAT('PV',I3,I2.2,I4.4,'N',I3,I2.2,I4.4,'W',6I6)
      RETURN
      END
****************************************************************************
      SUBROUTINE DSDA(LOCISN, LOCJSN, MIN1, MIN2, DS, DA)

** Compute change in distance and change in azimuth.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      HTF = 0.0D0
      HTT = 0.0D0

      READ(I4,REC=LOCISN,err=50,iostat=ios) YLATF,YLONF,VNF,VEF,VUF
      if (ios /= 0) goto 50
      READ(I4,REC=LOCJSN,err=51,iostat=ios) YLATT,YLONT,VNT,VET,VUT
      if (ios /= 0) goto 51

      CALL COMPSN(YLATF1,YLONF1,HTF1,YLATF,YLONF,HTF,
     1            MIN1, VNF, VEF, VUF)
      CALL COMPSN(YLATF2,YLONF2,HTF2,YLATF,YLONF,HTF,
     1            MIN2, VNF, VEF, VUF)
      CALL COMPSN(YLATT1,YLONT1,HTT1,YLATT,YLONT,HTT,
     1            MIN1, VNT, VET, VUT)
      CALL COMPSN(YLATT2,YLONT2,HTT2,YLATT,YLONT,HTT,
     1            MIN2, VNT, VET, VUT)

      CALL HELINV(YLATF1,YLONF1,YLATT1,YLONT1,FAZ1,BAZ1,S1)
      CALL HELINV(YLATF2,YLONF2,YLATT2,YLONT2,FAZ2,BAZ2,S2)
      DS = S2 - S1
      DA = FAZ2 - FAZ1
      RETURN

 50   write (*,'(/)') 
      write (*,*) "wrong input for the 1st reading in DSDA:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop
    
 51   write (*,'(/)') 
      write (*,*) "wrong input for the 2nd reading in DSDA:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop
    
      END
*********************************************************************
      SUBROUTINE HELINV(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
C                                                                       
C *** SOLUTION OF THE GEODETIC INVERSE PROBLEM AFTER T.VINCENTY.       
C *** MODIFIED RAINSFORD WITH HELMERT ELLIPTICAL TERMS.                 
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL. 
C *** STANDPOINT/FOREPOINT MUST NOT BE THE GEOGRAPHIC POLE .           
C
C   INPUT
C      GLAT1 = Latitude of from point (radians, positive north)
C      GLON1 = Longitude of from point (radians, positive west)
C      GLAT2 = Latitude of to point (radians, positive north)
C      GLON2 = Longitude of to point (radians, positive west)
C                        
C   OUTPUT
C      FAZ = Foward azimuth (radians, clockwise from north)
C      BAZ = Back azimuth   (radians, clockwise from north)
C      S   = Distance (meters)
C                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      IMPLICIT INTEGER*4 (I-N)
      COMMON/CONST/A,F,EPS2,EPS,AF,PI,TWOPI,RHOSEC                      
      DATA TOL/0.5D-14/                                                 
      R = 1.0D0-F                                                       
      TU1 = R*DSIN(GLAT1)/DCOS(GLAT1)                                   
      TU2 = R*DSIN(GLAT2)/DCOS(GLAT2)                                   
      CU1 = 1.0D0/DSQRT(TU1*TU1+1.0D0)                                  
      SU1 = CU1*TU1                                                     
      CU2 = 1.0D0/DSQRT(TU2*TU2+1.0D0)                                  
      S = CU1*CU2                                                       
      BAZ = S*TU2                                                       
      FAZ = BAZ*TU1                                                     
      X = GLON1-GLON2                                                   
  100 SX = DSIN(X)                                                      
      CX = DCOS(X)                                                      
      TU1 = CU2*SX                                                      
      TU2 = SU1*CU2*CX-BAZ                                              
      SY = DSQRT(TU1*TU1+TU2*TU2)                                       
      CY = S*CX+FAZ                                                     
      Y = DATAN2(SY,CY)                                                 
      SA = S*SX/SY                                                      
      C2A = -SA*SA+1.0D0                                                
      CZ = FAZ+FAZ                                                      
      IF(C2A .GT. 0.0D0)CZ = -CZ/C2A+CY                                 
      E = CZ*CZ*2.0D0-1.0D0                                             
      C = ((-3.0D0*C2A+4.0D0)*F+4.0D0)*C2A*F/16.0D0                     
      D = X                                                             
      X = ((E*CY*C+CZ)*SY*C+Y)*SA                                       
      X = (1.0D0-C)*X*F+GLON1-GLON2                                     
      IF(DABS(D-X) .GT. TOL)GO TO 100                                   
      FAZ = DATAN2(-TU1,TU2)                                            
      BAZ = DATAN2(CU1*SX,BAZ*CX-SU1*CU2)     
      FAZ = FAZ + PI
      BAZ = BAZ + PI                      
      IF(FAZ .LT. 0.0D0)FAZ = FAZ+TWOPI                             
      IF(BAZ .LT. 0.0D0)BAZ = BAZ+TWOPI   
      IF(FAZ .GT. TWOPI)FAZ = FAZ - TWOPI
      IF(BAZ .GT. TWOPI)BAZ = BAZ - TWOPI                            
      X = DSQRT((1.0D0/R/R-1.0D0)*C2A+1.0D0)+1.0D0                      
      X = (X-2.0D0)/X                                                   
      C = 1.0D0-X                                                       
      C = (X*X/4.0D0+1.0D0)/C                                           
      D = (0.375D0*X*X-1.0D0)*X                                         
      X = E*CY                                                          
      S = 1.0D0-E-E                                                     
      S = ((((SY*SY*4.0D0-3.0D0)*S*CZ*D/6.0D0-X)*D/4.0D0+CZ)*SY*D+Y)
     1    *C*A*R                                                        
      RETURN                                                            
      END
*************************************************************************

c     SUBROUTINE TODMSA(val,id,im,s)
 
*** convert position radians to deg,min,sec
*** range is [-twopi to +twopi]
 
c     implicit double precision(a-h,o-z)
c     IMPLICIT INTEGER*4 (I-N)
c     common/CONST/A,F,E2,EP2,AF,PI,TWOPI,RHOSEC
 
c   1 if(val.gt.twopi) then
c       val=val-twopi
c       go to 1
c     endif
 
c   2 if(val.lt.-twopi) then
c       val=val+twopi
c       go to 2
c     endif
 
c     if(val.lt.0.d0) then
c       isign=-1
c     else
c       isign=+1
c     endif
 
c     s=dabs(val*RHOSEC/3600.D0)
c     id=idint(s)
c     s=(s-id)*60.d0
c     im=idint(s)
c     s=(s-im)*60.d0
 
*** account for rounding error
 
c     is=idnint(s*1.d5)
c     if(is.ge.6000000) then
c       s=0.d0
c       im=im+1
c     endif
c     if(im.ge.60) then
c       im=0
c       id=id+1
c     endif
 
c     id = isign*id
c     im = isign*im
c     s  = isign*s

c     return
c     end
*********************************************************
      SUBROUTINE TRFDAT(CARD,DATE,IREC12,IYEAR1,IYEAR2,MINS)

C Convert blue-book date to time in minutes

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER    CARD*80
      CHARACTER   DATE*6
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

      IF (IREC12 .EQ. 0) THEN
         WRITE (LUOUT, 5)
    5    FORMAT(' ABORT: The blue-book needs a valid *12* record.')
         STOP
      ENDIF

      READ (DATE, 10,err=50,iostat=ios) IYEAR, MONTH, IDAY
      if (ios /= 0) goto 50
   10 FORMAT (3I2)

      IF ( IYEAR1 .LE. (1900 + IYEAR) .AND.
     1     (1900 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 1900 + IYEAR
      ELSEIF ( IYEAR1 .LE. (2000 + IYEAR) .AND.
     1     (2000 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 2000 + IYEAR
      ELSEIF ( IYEAR1 .LE. (1800 + IYEAR) .AND.
     1     (1800 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 1800 + IYEAR
      ELSE
         WRITE (LUOUT, 20) CARD
   20    FORMAT(' ABORT: The following record has a date'/
     1          ' which is inconsistent with the *12* record'/
     1          3x, A80)
      ENDIF

      IF (IYEAR .LE. 1906) THEN
         WRITE (LUOUT, 30)
   30    FORMAT(' ***WARNING***'/
     1   ' The blue-book file contains an observation that'/
     1   ' predates 1906.  The TDP model may not be valid'/
     1   ' and the computed corrections may be erroneous.')
      ENDIF
   
      IF (IDAY .EQ. 0) IDAY = 15

      IF (MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF

C     CALL TOTIME(IYEAR,MONTH,IDAY, MINS)
      CALL IYMDMJ(IYEAR,MONTH,IDAY, MJD)
      MINS = MJD * 24 * 60
      RETURN

  50  write (*,'(/)') 
      write (*,*) 'wrong IYEAR, IMONTH, IDAY:ios =',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop


      END
*********************************************************
      SUBROUTINE TNFDAT(DATE,MINS)

C Convert blue-book date to time in years

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      character DATE*8
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      READ(DATE,10) IYEAR,MONTH,IDAY
   10 FORMAT(I4,I2,I2)
      IF(IYEAR .LE. 1906) THEN
        WRITE(LUOUT,20)
   20     FORMAT(' ***WARNING***'/
     1    ' The blue-book file contains an observation that'/
     2    ' predates 1906.  The TDP model is not valid and the'/
     3    ' computed correction may be erroneous.')
      ENDIF
      IF(IDAY .EQ. 0) IDAY = 15
      IF(MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF
C     CALL TOTIME(IYEAR,MONTH,IDAY,MINS)
      CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
      MINS = MJD * 24 * 60
      RETURN
      END
***************************************************************
C     SUBROUTINE GETTIM(MONTH, IDAY, IYEAR, DATE, MINS, TEST)
C
*** Read month-day-year and convert to decimal years
*** and Julian time in minutes      
***    MONTH      input - number from 1 to 12
***    IDAY       input - number from 1 to 31
***    IYEAR      input - must be after 1906
***    DATE       output - corresponding time in decimal years
***    MINS       output - corresponding julian time in minutes
***    TEST       output - if (true) then there is an error
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT INTEGER*4 (I-N)
C     LOGICAL TEST
C     COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

C     READ(LUIN,*) MONTH,IDAY,IYEAR
    
C     IF(IYEAR .le. 1906) THEN
C         WRITE(LUOUT,10)
C  10     FORMAT(' The model is not valid for dates prior ',
C    1           'to 1906.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     IF(MONTH .le. 0 .or. MONTH .gt. 12) THEN
C         WRITE(LUOUT,20)
C  20     FORMAT(' Improper month specified.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     IF(IDAY .le. 0 .or. IDAY .gt. 31) THEN
C         WRITE(LUOUT,30)
C  30     FORMAT(' Improper day specified.'/)
C         TEST = .TRUE.
C         RETURN
C     ENDIF

C     CALL TOTIME(IYEAR, MONTH, IDAY, MINS)
C     CALL TOTIME(IYEAR, 1, 1, MIN00)
C     DATE = DBLE(IYEAR) + DBLE(MINS - MIN00)/525600.D0
C     TEST = .FALSE.
C     RETURN
C     END

C************************************************************************
      SUBROUTINE TOXYZ(glat,glon,eht,x,y,z)
 
*** compute x,y,z
*** ref p.17 geometric geodesy notes vol 1, osu, rapp
 
      implicit double precision(a-h,o-z)
      common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
 
      slat=dsin(glat)
      clat=dcos(glat)
      w=dsqrt(1.d0-e2*slat*slat)
      en=a/w
 
      x=(en+eht)*clat*dcos(glon)
      y=(en+eht)*clat*dsin(glon)
      z=(en*(1.d0-e2)+eht)*slat
 
      return
      end
C************************************************************************
      logical function FRMXYZ(x,y,z,glat,glon,eht)
 
*** convert x,y,z into geodetic lat, lon, and ellip. ht
*** ref: eq a.4b, p. 132, appendix a, osu #370
*** ref: geom geod notes gs 658, rapp
 
      implicit double precision(a-h,o-z)
      parameter(maxint=10,tol=1.d-13)
      common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
 
      ae2=a*e2
 
*** compute initial estimate of reduced latitude  (eht=0)
 
      p=dsqrt(x*x+y*y)
      icount=0
      tgla=z/p/(1.d0-e2)
 
*** iterate to convergence, or to max # iterations
 
    1 if(icount.le.maxint) then
        tglax=tgla
        tgla=z/(p-(ae2/dsqrt(1.d0+(1.d0-e2)*tgla*tgla)))
        icount=icount+1
        if(dabs(tgla-tglax).gt.tol) go to 1
 
*** convergence achieved
 
        frmxyz=.true.
        glat=datan(tgla)
        slat=dsin(glat)
        clat=dcos(glat)
        glon=datan2(y,x)
        w=dsqrt(1.d0-e2*slat*slat)
        en=a/w
        if(dabs(glat).le.0.7854d0) then
          eht=p/clat-en
        else
          eht=z/slat-en+e2*en
        endif
        glon=datan2(y,x)
 
*** too many iterations
 
      else
        frmxyz=.false.
        glat=0.d0
        glon=0.d0
        eht=0.d0
      endif
 
      return
      end
*********************************************************************
      SUBROUTINE RADII(YLAT,RADMER,RADPAR)
C
C  Computes the radius of curvature in the meridian
C  and the radius of curvature in a parallel of latitude
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COSLAT = DCOS(YLAT)
      DENOM = DSQRT(1.D0 + EPS*COSLAT*COSLAT)
      RADMER = AF/(DENOM**3)
      RADPAR = AF*COSLAT/DENOM
      RETURN
      END
*********************************************************************
      SUBROUTINE GETGRD(NAMEG,MINLAT,MAXLAT,ILAT,MINLON,MAXLON,ILON)

*** Interactively obtain specifications for a grid.

      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER   NAMEG*10
      WRITE(LUOUT,200)
  200     FORMAT(' Enter name for grid (10 character max).  ')
      READ(LUIN,210,err=50,iostat=ios) NAMEG
          if (ios /= 0) goto 50
  210     FORMAT(A10)
      WRITE(LUOUT,220)
  220     FORMAT(' Enter minimum latitude for grid in deg-min-sec'/
     1           ' in free format (integer values only).  ')
      READ(LUIN,*,err=51,iostat=ios) ID1, IM1, IS1
          if (ios /= 0) goto 51
      WRITE(LUOUT,230)
  230     FORMAT(' Enter maximum latitude in same format.  ')
      READ(LUIN,*,err=52,iostat=ios) ID2, IM2, IS2
          if (ios /= 0) goto 52
      WRITE(LUOUT,240)
  240     FORMAT(' Enter latitude increment in seconds ',
     1           ' (integer value).  ')
      READ(LUIN,*,err=53,iostat=ios) ILAT
          if (ios /= 0) goto 53
      WRITE(LUOUT,250)
  250     FORMAT(' Enter minimum longitude with positive being west.  ')
      READ(LUIN,*,err=54,iostat=ios) JD1, JM1, JS1
          if (ios /= 0) goto 54
          WRITE(LUOUT,260)
  260     FORMAT(' Enter maximum longitude.  ')
      READ(LUIN,*,err=55,iostat=ios) JD2, JM2, JS2
          if (ios /= 0) goto 55
          WRITE(LUOUT,270)
  270     FORMAT(' Enter longitude increment in seconds ',
     1           ' (integer value).  ')
      READ(LUIN,*,err=56,iostat=ios) ILON
          if (ios /= 0) goto 56
      MINLAT = 3600*ID1 + 60*IM1 + IS1
      MAXLAT = 3600*ID2 + 60*IM2 + IS2
      MINLON = 3600*JD1 + 60*JM1 + JS1
      MAXLON = 3600*JD2 + 60*JM2 + JS2
      RETURN

  50      write (*,'(/)') 
          write (*,*) "Wrong grid name in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  51      write (*,'(/)') 
          write (*,*) "Wrong min lat in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  52      write (*,'(/)') 
          write (*,*) "Wrong max lat in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  53      write (*,'(/)') 
          write (*,*) "Wrong ILAT in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  54      write (*,'(/)') 
          write (*,*) "Wrong min lon in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  55      write (*,'(/)') 
          write (*,*) "Wrong max lon in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  56      write (*,'(/)') 
          write (*,*) "Wrong ILON in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

      END
*********************************************************************
      SUBROUTINE GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)

*** Interactively obtain the specifications for a line.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT,I1,I2,I3,I4,I5,I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      CHARACTER   NAMEG*10

      WRITE(LUOUT,100)
  100 FORMAT(' Enter name for line (10 character max.)  ')
      READ(LUIN,110,err=50,iostat=ios) NAMEG
      if (ios /= 0) goto 50
  110 FORMAT(A10)
      WRITE(LUOUT,120)
  120 Format(' Specify the latitude for the origin of the line'/
     1       ' in deg-min-sec in free format. For example'/
     2       '         35 17 28.3'/
     3       ' Positive is north.  ')
      READ(LUIN,*,err=51,iostat=ios) LATD,LATM,SLAT
      if (ios /= 0) goto 51
      XLAT = (DBLE(3600*LATD + 60*LATM)+SLAT)/RHOSEC
      WRITE(LUOUT,130)
  130 FORMAT(' Specify the longitude for the origin of the line'/
     1       ' in free format with west being positive.  ')
      READ(LUIN,*,err=52,iostat=ios) LOND,LONM,SLON
      if (ios /= 0) goto 52
      XLON = (DBLE(3600*LOND+60*LONM)+SLON)/RHOSEC
      WRITE(LUOUT,140)
  140 FORMAT(' Specify the orientation of the line clockwise from'/
     1       ' north in decimal degrees, eg., 43.7.  ')
      READ(LUIN,*,err=53,iostat=ios) FAZ
      if (ios /= 0) goto 53
      FAZ = 3600.D0*FAZ/RHOSEC
      BAZ = FAZ + PI
      IF(BAZ .GT. TWOPI) BAZ = BAZ - TWOPI
      WRITE(LUOUT,150)
  150 FORMAT(' Specify minimum and maximum distance from origin in'/
     1       ' meters, eg., -40000. 30000. '/
     2       ' NOTE: negative distance corresponds to distance in'/
     3       '       the opposite direction.  ')
      READ(LUIN,*,err=54,iostat=ios) XMIN,XMAX
      if (ios /= 0) goto 54
      WRITE(LUOUT,160)
  160 FORMAT(' Specify distance increment in meters.  ')
      READ(LUIN,*,err=55,iostat=ios) XINC
      if (ios /= 0) goto 55
      RETURN

  50  write (*,'(/)') 
      write (*,*) 'Wrong name of line in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  51  write (*,'(/)') 
      write (*,*) 'Wrong LATD,LATM,SLAT in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  52  write (*,'(/)') 
      write (*,*) 'Wrong LOND,LONM,SLON in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  53  write (*,'(/)') 
      write (*,*) 'Wrong FAZ in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  54  write (*,'(/)') 
      write (*,*) 'Wrong min and max D in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  55  write (*,'(/)') 
      write (*,*) 'Wrong increment in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*******************************************************************
      SUBROUTINE DIRCT1(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
C
C *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
C *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
C
C *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
C *** FINV IS THE FLATTENING OF THE REFERENCE ELLIPSOID
C *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
C *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
C *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
C
C *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
C *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,FINV,E2,EPI,AF,PI,TWOPI,RHOSEC
      DATA EPS/0.5D-13/
      R=1.D0-FINV
      TU=R*DSIN(GLAT1)/DCOS(GLAT1)
      SF=DSIN(FAZ)
      CF=DCOS(FAZ)
      BAZ=0.D0
      IF(CF.NE.0.D0) BAZ=DATAN2(TU,CF)*2.D0
      CU=1.D0/DSQRT(TU*TU+1.D0)
      SU=TU*CU
      SA=CU*SF
      C2A=-SA*SA+1.D0
      X=DSQRT((1.D0/R/R-1.D0)*C2A+1.D0)+1.D0
      X=(X-2.D0)/X
      C=1.D0-X
      C=(X*X/4.D0+1.D0)/C
      D=(0.375D0*X*X-1.D0)*X
      TU=S/R/A/C
      Y=TU
  100 SY=DSIN(Y)
      CY=DCOS(Y)
      CZ=DCOS(BAZ+Y)
      E=CZ*CZ*2.D0-1.D0
      C=Y
      X=E*CY
      Y=E+E-1.D0
      Y=(((SY*SY*4.D0-3.D0)*Y*CZ*D/6.D0+X)*D/4.D0-CZ)*SY*D+TU
      IF(DABS(Y-C).GT.EPS)GO TO 100
      BAZ=CU*CY*CF-SU*SY
      C=R*DSQRT(SA*SA+BAZ*BAZ)
      D=SU*CY+CU*SY*CF
      GLAT2=DATAN2(D,C)
      C=CU*CY-SU*SY*CF
      X=DATAN2(SY*SF,C)
      C=((-3.D0*C2A+4.D0)*FINV+4.D0)*C2A*FINV/16.D0
      D=((E*CY*C+CZ)*SY*C+Y)*SA
      GLON2=GLON1+X-(1.D0-C)*D*FINV
      IF (GLON2.GE.TWOPI) GLON2=GLON2-TWOPI
      IF(GLON2.LT.0.D0) GLON2=GLON2+TWOPI
      BAZ=DATAN2(SA,BAZ)+PI
      IF (BAZ.GE.TWOPI) BAZ=BAZ-TWOPI
      IF (BAZ.LT.0.D0) BAZ=BAZ+TWOPI
      RETURN
      END
***********************************************************
      SUBROUTINE TOCHAR(ORIG,CHAR14)

*** Convert double precision real number to character*14

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER   CHAR14*14

      ORIG = ORIG*10000.D0
      WRITE(CHAR14,10) ORIG
   10 FORMAT(F14.0)
      RETURN
      END
*********************************************************

      SUBROUTINE DDXYZ(ISN, JSN, MIN1, MIN2, 
     1                 DDX, DDY, DDZ)

** Compute change in DX,DY,DZ-vector from time MIN1 to MIN2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      character      PIDs*6
      parameter (nbbdim = 10000)
      COMMON /ARRAYS/ HT(nbbdim), LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      READ(I4,REC=LOC(ISN),err=50,iostat=ios) GLATI,GLONI,VNI,VEI,VUI
      if (ios /= 0) goto 50
      READ(I4,REC=LOC(JSN),err=51,iostat=ios) GLATJ,GLONJ,VNJ,VEJ,VUJ
      if (ios /= 0) goto 51

      CALL COMPSN(YLATI1,YLONI1,HTI1,GLATI,GLONI,HT(ISN),
     1            MIN1,VNI, VEI, VUI)
      CALL COMPSN(YLATI2,YLONI2,HTI2,GLATI,GLONI,HT(ISN),
     1            MIN2, VNI, VEI, VUI)
      CALL COMPSN(YLATJ1,YLONJ1,HTJ1,GLATJ,GLONJ,HT(JSN),
     1            MIN1, VNJ, VEJ, VUJ)
      CALL COMPSN(YLATJ2,YLONJ2,HTJ2,GLATJ,GLONJ,HT(JSN),
     1            MIN2, VNJ, VEJ, VUJ)

      XLONI1 = -YLONI1
      XLONI2 = -YLONI2
      XLONJ1 = -YLONJ1
      XLONJ2 = -YLONJ2

      CALL TOXYZ(YLATI1,XLONI1,HTI1,XI1,YI1,ZI1)
      CALL TOXYZ(YLATI2,XLONI2,HTI2,XI2,YI2,ZI2)
      CALL TOXYZ(YLATJ1,XLONJ1,HTJ1,XJ1,YJ1,ZJ1)
      CALL TOXYZ(YLATJ2,XLONJ2,HTJ2,XJ2,YJ2,ZJ2)

      DDX = (XJ2 - XI2) - (XJ1 - XI1)
      DDY = (YJ2 - YI2) - (YJ1 - YI1)
      DDZ = (ZJ2 - ZI2) - (ZJ1 - ZI1)

      RETURN

 50   write (*,'(/)') 
      write (*,*) "Failed in 1st read statement in DDXYZ:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 51   write (*,'(/)') 
      write (*,*) "Failed in 2nd read statement in DDXYZ:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*************************************************************
      SUBROUTINE DISLOC (YLAT,YLON,STRIKE,HL,EQLAT,EQLON,
     &          SS,DS,DIP,DEPTH,WIDTH,DNORTH,DWEST,DUP)

*** Compute 3-dimensional earthquake displacement at point
*** using dislocation theory
*
*   INPUT:
*         YLAT = Latitude in radians (positive north)
*         YLON = Longitude in radians (positive west)
*         STRIKE = strike in radians clockwise from north such
*                  that the direction of dip is pi/2 radians
*                  counterclockwise from the direction of strike
*         HL   = Half-length in meters
*         EQLAT = Latitude in radians of midpoint of the
*                 rectangle's upper edge (positive north)
*         EQLON = Longitude in radians of midpoint of the
*                 rectangle's upper edge (positive west)
*         SS = strike slip in meters (positive = right lateral)
*         DS = dip slip  in meters (positive = normal faulting)
*         DIP = dip in radians
*         DEPTH = Vertical depth of rectangle's upper edge
*                 in meters
*         WIDTH = width of rectangle in meters
*
*   OUTPUT:
*         DNORTH = northward displacement in radians
*         DWEST = westward displacement in radians
*         DUP = upward displacement in meters
************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

*** Compute radii of curvature at fault center
      CALL RADII (EQLAT, RMER, RPAR)

*** Compute planar coordinates in meters
      DLAT = (YLAT - EQLAT) * RMER
      DLON = (YLON - EQLON) * RPAR
      COSSTR = DCOS(STRIKE)
      SINSTR = DSIN(STRIKE)
      X1 = COSSTR*DLAT - SINSTR*DLON
      X2 = SINSTR*DLAT + COSSTR*DLON

*** Compute displacements in fault-oriented coordinates
      CALL OKADA(X1,X2,HL,DEPTH,WIDTH,DIP,U1SS,U2SS,
     &      U3SS,U1DS,U2DS,U3DS)
      U1 = U1SS*SS + U1DS*DS
      U2 = U2SS*SS + U2DS*DS
      DUP = U3SS*SS + U3DS*DS

*** Convert horizontal displacements to radians
*** in north-west coordinate system
      DNORTH = ( COSSTR*U1 + SINSTR*U2) / RMER
      DWEST  = (-SINSTR*U1 + COSSTR*U2) / RPAR

      RETURN
      END
****************************************************
      SUBROUTINE OKADA(X1,X2,XL,DU,W,DIP,
     1     U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

************************************************************
*  This subroutine computes displacements at the point X1,X2
*  on the Earth's surface due to 1.0 meter of right-lateral 
*  strike slip (SS) and 1.0 meter of normal dip slip (DS) 
*  along a rectangular fault.
*
*  The rectangular fault dips in the direction of the positive
*  X2-axis.  The rectangle's strike parallels the X1-axis.
*  With the X3-axis directed upward out of the Earth, the X1-,
*  X2-, and X3-axes form a right-handed system.
*
*  The equations of dislocation theory are employed whereby
*  Earth is represented an a homogeneous, isotropic half-space
*  with a Poisson ratio of PNU.
*
*  REFERENCE: Okada, Y., Surface deformation due to shear and
*    tensile faults in a half-space, Bulletin of the 
*    Seismological Society of America, vol. 75, pp. 1135-1154 (1985)
*
*  The X3 = 0 plane corresponds to the Earth's surface. The plane's
*  origin is located directly above the midpoint of the rectangle's
*  upper edge.
*
*  INPUT:
*    X1,X2 - Location in meters
*    XL    - Rectangle's half-length in meters
*    DU    - Vertical depth to rectangle's upper edge in meters
*            (always positive or zero)
*    W     - Rectangle's width in meters
*    DIP   - Rectangle's dip in radians (always between 0 and PI/2
*
*  OUTPUT
*    U1SS  - Displacement in X1-direction due to 1.0 meters
*            of right-lateral strike slip
*    U2SS  - Displacement in X2-direction due to 1.0 meters
*            of right-lateral strike slip
*    U3SS  - Displacement in X3-direction due to 1.0 meters
*            of right-lateral strike slip
*    U1DS  - Displacement in X1-direction due to 1.0 meters
*            of normal dip slip
*    U2DS  - Displacement in X2-direction due to 1.0 meters
*            of normal dip slip
*    U3DS  - Displacement in X3-direction due to 1.0 meters
*            of normal dip slip
*******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL VERT   
      PI = 3.141593D0
      TWOPI = PI + PI
      PNU = 0.25D0
      RATIO = 1.D0 - 2.D0*PNU

      IF(DABS(PI/2.D0 - DIP) .LT. .01D0)THEN
               DIPK = -PI/2.D0
               VERT = .TRUE.
      ELSE
               DIPK = -DIP
               VERT = .FALSE.
      ENDIF

      SDIP = DSIN(DIPK)
      CDIP = DCOS(DIPK)
      P = X2*CDIP + DU*SDIP
      Q = X2*SDIP - DU*CDIP

      PSI = X1 + XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1    VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

      PSI = X1 + XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS 
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS + C1SS
      U2SS = U2SS + C2SS
      U3SS = U3SS + C3SS
      U1DS = U1DS + C1DS
      U2DS = U2DS + C2DS
      U3DS = U3DS + C3DS
      RETURN
      END
*************************************************************
      SUBROUTINE OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1           VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL VERT   

      YBAR = ETA*CDIP + Q*SDIP
      DBAR = ETA*SDIP - Q*CDIP
      R = DSQRT(PSI*PSI + ETA*ETA + Q*Q)
      X = DSQRT(PSI*PSI + Q*Q)
      IF(DABS(Q) .LE. 0.1d0) THEN
         TERM = 0.D0
      ELSE
         TERM = DATAN(PSI*ETA/(Q*R))
      ENDIF

      IF(VERT) THEN
         F5 = -RATIO*PSI*SDIP/(R + DBAR)
         F4 = -RATIO*Q/(R + DBAR)
         F3 = 0.5D0*RATIO*(ETA/(R + DBAR)
     1          + YBAR*Q/((R + DBAR)*(R + DBAR))
     2          - DLOG(R + ETA))
         F1 = -0.5D0*RATIO*PSI*Q/
     1        ((R + DBAR)*(R + DBAR))
      ELSE
         IF(DABS(PSI) .LE. 0.1D0) then
            F5 = 0.d0
         ELSE
            F5 = 2.D0*RATIO*
     1      DATAN((ETA*(X+Q*CDIP)+X*(R+X)*SDIP)/(PSI*(R+X)*CDIP))
     2          /CDIP
         ENDIF
         F4 = RATIO*(DLOG(R+DBAR)-SDIP*DLOG(R+ETA))/CDIP
         F3 = RATIO*(YBAR/(CDIP*(R+DBAR)) - DLOG(R+ETA))
     1         + SDIP*F4/CDIP
         F1 = -RATIO*(PSI/(CDIP*(R+DBAR))) - SDIP*F5/CDIP
      ENDIF
         F2 = -RATIO*DLOG(R+ETA) - F3

      U1SS = -(PSI*Q/(R*(R+ETA))
     1         + TERM + F1*SDIP)/TWOPI
      U2SS = -(YBAR*Q/(R*(R+ETA))
     1         + Q*CDIP/(R+ETA)
     2         + F2*SDIP)/TWOPI
      U3SS = -(DBAR*Q/(R*(R+ETA))
     1         + Q*SDIP/(R+ETA)
     2         + F4*SDIP)/TWOPI
      U1DS = -(Q/R - F3*SDIP*CDIP)/TWOPI
      U2DS = -(YBAR*Q/(R*(R+PSI))
     1         + CDIP*TERM - F1*SDIP*CDIP)/TWOPI
      U3DS = -(DBAR*Q/(R*(R+PSI))
     1         + SDIP*TERM - F5*SDIP*CDIP)/TWOPI
      RETURN
      END
*******************************************************************

C     SUBROUTINE GRDCHK (POSX, POSY, INSIDE)

C
C ROUTINE CHECKS IF THE POINT HAVING COORDINATES (POSX, POSY)
C IS WITHIN THE REGION SPANNED BY THE GRID
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT INTEGER*4 (I-N)
C     LOGICAL INSIDE
C     COMMON /CDGRID/ GRDLX, GRDUX, GRDLY, GRDUY, ICNTX, ICNTY

C     INSIDE = .TRUE.

C     IF (POSX .LT. GRDLX .OR. POSX .GT. GRDUX) THEN
C        INSIDE = .FALSE.
C     ENDIF
C     IF (POSY .LT. GRDLY .OR. POSY .GT. GRDUY) THEN
C        INSIDE = .FALSE.
C     ENDIF

C     RETURN
C     END

*******************************************************************

      SUBROUTINE GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDWEI
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C YLON         LONGITUDE OF POINT IN RADIANS, POSITIVE WEST
C YLAT         LATITUDE OF POINT IN RADIANS, POSITIVE NORTH
C JREGN        ID OF GEOGRAPHIC REGION CONTAINING POINT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
C              FOR THE CORNER VECTORS
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   COMVEL
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /CDGRID/, /CONST/
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9511.09, RAS, MODIFIED FOR HTDP
C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
C********1*********2*********3*********4*********5*********6*********7**
    
C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 10)
      DIMENSION WEI(2,2)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD), 
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Convert input coordinates to degrees
      POSX = (TWOPI - YLON) * 180.D0 / PI
      POSY = YLAT * 180.D0 / PI

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (GRDUX(JREGN) - GRDLX(JREGN)) / ICNTX(JREGN)
      STEPY = (GRDUY(JREGN) - GRDLY(JREGN)) / ICNTY(JREGN)
      I = IDINT((POSX - GRDLX(JREGN))/STEPX) + 1
      J = IDINT((POSY - GRDLY(JREGN))/STEPY) + 1
c     write(6,1001) JREGN, I, J
c1001 format(1x, 'jregn = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = GRDLX(JREGN) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = GRDLY(JREGN) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END

C*********************************************************************
C
      SUBROUTINE GRDVEC (JREGN, I, J, VEL, B)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDVEC
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE APPROXIMATE VALUES OF THE
C              GRID NODE VELOCITIES FOR GRID (I,J) 
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C JREGN        ID OF GEOGRAPHIC REGION CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C B            THE ARRAY CONTAINING ALL THE APPROXIMATE VALUES
C              FOR THE ADJUSTMENT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C VEL          A TWO BY TWO ARRAY CONTAINING THE VELOCITY VECTORS
C              FOR THE CORNERS OF THE GRID
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   COMVEL
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9712.05, RAS, MODIFIED FOR HTDP (version 2.2)
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION VEL(2,2,3), B(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IUNGRD(JREGN, I + II, J + IJ, IVEC)
               VEL(II + 1, IJ + 1, IVEC) = B(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

      RETURN
      END

C***************************************************

      SUBROUTINE RDEG (INPUT,VAL,CNEG)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER  INPUT*10
      CHARACTER  CNEG*1
      INTEGER      DEG, MIN, ISEC

      DO 10 I = 1, 9
         IF (INPUT(I:I).EQ.' ') THEN
             INPUT(I:I) = '0'
         ENDIF
   10 CONTINUE

      READ (INPUT,20) DEG, MIN, ISEC

   20 FORMAT(I3,I2,I4)

      SEC = ISEC/100.D0

      VAL = (DEG + (MIN/60.D0) + (SEC/3600.D0) ) 

      IF (INPUT(10:10).EQ.CNEG) THEN
         VAL = -VAL
      ENDIF

      RETURN
      END

C***************************************************

      INTEGER FUNCTION IUNGRD(IREGN, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 10)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)

      IUNGRD = NBASE(IREGN) +
     1      3 * ((J - 1) * (ICNTX(IREGN) + 1) +  (I - 1)) + IVEC

      RETURN
      END
     
C*********************************************************
C
      SUBROUTINE TOMNT( IYR, IMON, IDAY, IHR, IMN, MINS )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       TOMNT (ORIGINALLY IYMDMJ)
C VERSION:    9004.17
C WRITTEN BY: M. SCHENEWERK
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE PLUS UT
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IDAY              DAY
C IMON              MONTH
C IYR               YEAR
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C MINS              MODIFIED JULIAN DATE IN MINUTES
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C A                 TEMPORARY STORAGE
C B                 TEMPORARY STORAGE
C C                 TEMPORARY STORAGE
C D                 TEMPORARY STORAGE
C IMOP              TEMPORARY STORAGE
C IYRP              TEMPORARY STORAGE
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:     DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:       
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
C                              I.E. 1992 RATHER THAN 92.  
C
C********1*********2*********3*********4*********5*********6*********7**
C::LAST MODIFICATION
C::8909.06, MSS, DOC STANDARD IMPLIMENTED
C::9004.17, MSS, CHANGE ORDER YY MM DD
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      INTEGER*4     A, B, C, D

      IYRP = IYR
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IF( IMON .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMON + 12
      ELSE
        IMOP= IMON
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP*0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0*IYRP
      D=  30.6001D0*(IMOP + 1)
      MINS =  (B + C + D + IDAY - 679006) * (24 * 60)
     &       + (60 * IHR) + IMN
C      
      RETURN
      END
*****************************************************
      SUBROUTINE TOVNEU(GLAT,GLON,VX,VY,VZ,VN,VE,VU)

*** Convert velocities from vx,vy,vz to vn,ve,vu

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VN = -SLAT*CLON*VX - SLAT*SLON*VY + CLAT*VZ
      VE = -SLON*VX + CLON*VY
      VU = CLAT*CLON*VX + CLAT*SLON*VY + SLAT*VZ

      RETURN
      END
***************************************************
      SUBROUTINE TOVXYZ(GLAT,GLON,VN,VE,VU,VX,VY,VZ)
*** Convert velocities from vn,ve,vu to vx,vy,vz
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VX = -SLAT*CLON*VN - SLON*VE + CLAT*CLON*VU
      VY = -SLAT*SLON*VN + CLON*VE + CLAT*SLON*VU
      VZ =  CLAT*VN + SLAT*VU

      RETURN
      END
*****************************************************
      SUBROUTINE DPLACE

*** Predict displacements between two dates.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 16)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER   CARD*80
      CHARACTER   record*120
      CHARACTER   NAMEF*30,NAME*30,NAMEBB*30, NAMEIF*30
      CHARACTER   NAME24
      CHARACTER   BLAB*17
      CHARACTER   NAMEG*10
      CHARACTER   TYPE*4
      CHARACTER   OPTION*1,JN*1,JW*1, VOPT*1
      CHARACTER   LATDIR*1, LONDIR*1
      CHARACTER   ANSWER*1
      character   frame1*24
      LOGICAL TEST

      BLAB = 'OUTSIDE OF REGION'

      WRITE(LUOUT,10)
   10 FORMAT(
     1   ' Displacements will be predicted from time T1 to time T2.')
      WRITE(LUOUT,* ) ' Please enter T1 '
   20 CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
      IF(TEST) then
         write(luout,*) ' Do you wish to re-enter T1? (y/n)'
         read (luin,21,err=600,iostat=ios) ANSWER
         if (ios /= 0) goto 600
   21    format( A1 )
         IF (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 20
         RETURN
      ENDIF
      WRITE(LUOUT,*) ' Please enter T2 '
   35 CALL GETMDY(MONTH2,IDAY2,IYEAR2,DATE2,MIN2,TEST)
      IF(TEST) then
         write(luout,*) ' Do you wish to re-enter T2? (y/n) '
         read(luin,21,err=600,iostat=ios) ANSWER
         if (ios /= 0) goto 600
         if (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 35
         RETURN
      ENDIF

      WRITE(LUOUT,45)
   45 FORMAT(
     1    ' Please enter the name for the file to contain'/
     2    ' the predicted displacements.  ')
      READ(LUIN,50,err=601,iostat=ios) NAMEF
      if (ios /= 0) goto 601
   50 FORMAT(A30)
      OPEN(I2,FILE=NAMEF,STATUS='UNKNOWN')
      CALL HEADER

*** Choosing reference frame for displacements
   56 WRITE(LUOUT,55)
   55 FORMAT(' ************************************************'/
     1   ' Select the reference frame to be used for specifying'/
     2   ' positions and displacements.  '/)
      call MENU1(iopt, frame1)

      IF(IOPT .GE. 1 .AND . IOPT .LE. numref) THEN
          WRITE(I2,57) frame1
   57     FORMAT(' DISPLACEMENTS IN METERS RELATIVE TO ', a24)
      ELSE
          WRITE(LUOUT,70)
   70     FORMAT(' Improper selection--try again.  ')
          GO TO 56
      ENDIF
            
      WRITE(I2,71) MONTH1,IDAY1,IYEAR1,MONTH2,IDAY2,IYEAR2,
     1             DATE1, DATE2
   71 FORMAT (
     1  ' FROM ',I2.2,'-',I2.2,'-',I4,' TO ',
     2           I2.2,'-',I2.2,'-',I4,' (month-day-year)'/
     2  ' FROM ',F8.3, ' TO ',F8.3, ' (decimal years)'//
     3  'NAME OF SITE             LATITUDE          LONGITUDE    ',
     4        '        NORTH    EAST    UP ')
   75 WRITE(LUOUT,80)
   80 FORMAT(' ********************************'/
     1   ' Displacements will be predicted at each point whose',/
     2   ' horizontal position is specified.',/
     4   ' Please indicate how you wish to supply positions.'/ )
   85 WRITE(LUOUT,86)
   86 FORMAT(
     5   '    0. No more points. Return to main menu.'/
     6   '    1. Individual points entered interactively.'/
     7   '    2. Points on a specified grid.'/
     8   '    3. The *80* records in a specified blue-book file.'/
     9   '    4. Points on a specified line.  ' /
     1   '    5. Batch file of delimited records of form: ' /
     2   '       LAT,LON,TEXT ' /
     4   '       LAT = latitude in degrees (positive north/DBL PREC)' /
     5   '       LON = longitude in degrees (positive west/DBL PREC)' /
     7   '       TEXT = Descriptive text (CHARACTER*24) ' /
     8   '       Example:  ' /
     9   '       40.731671553,112.212671753,SALT AIR '/)
      READ(LUIN,'(A1)',err=602,iostat=ios) OPTION
      if (ios /= 0) goto 602

      IF(OPTION .EQ. '0') THEN
        GO TO 510               
      ELSEIF(OPTION .EQ. '1') THEN
        CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1        LONDIR,NAME24, X,Y,Z,YLAT,YLON,EHT)
        ELON = - YLON
          call GETVLY(YLAT,ELON,VX,VY,VZ,VN,VE,VU,VOPT,210)
        if (vopt .eq. '0') then
          call PREDV( ylat, ylon, eht, date1, iopt,
     1                jregn, vn, ve, vu)
          if (jregn .eq. 0) then
            write(luout, 140)
            go to 85 
          endif
        endif
        call NEWCOR(ylat, ylon, eht, min1, min2, 
     1     ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
  140   FORMAT(' ****************************************'/
     1         ' A displacement can not be predicted because'/
     1         ' the point is outside of the modeled region.'/
     2         ' For additional displacements, please indicate how'/
     3         ' you wish to supply the horizontal coordinates.'/)
        WRITE(LUOUT,150) DN, DE,DU
  150   FORMAT(' *************************************'/
     1         ' Northward displacement = ',F7.3,' meters.'/
     1         ' Eastward displacement  = ',F7.3,' meters.'/
     1         ' Upward displacement    = ',F7.3,' meters.'/
     1         ' ****************************************'//
     2         ' For additional displacements, please indicate how'/
     3         ' you wish to supply the horizontal coordinates.'/)
        WRITE(I2,160) NAME24,LATD,LATM,SLAT,LATDIR,
     1        LOND,LONM,SLON,LONDIR,DN,DE,DU
  160   FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,
     1         I3,1X,I2,1X,F8.5,1X,A1,1X,3F8.3)
      ELSEIF(OPTION .EQ. '2') THEN
        CALL GETGRD(NAMEG,MINLAT,MAXLAT,IDS,MINLON,MAXLON,JDS)
        I = -1
  280   I = I + 1
        LAT = MINLAT + I*IDS
        IF(LAT .GT. MAXLAT) GO TO 296
        XLAT = DBLE(LAT)/RHOSEC
        CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
        LATDIR = 'N'
        IF (ISIGN .eq. -1) LATDIR = 'S'
        J = -1
  290   J = J + 1
        YLAT = XLAT
        LON = MINLON + J*JDS
        IF(LON .GT. MAXLON) GO TO 280
        YLON = DBLE(LON)/RHOSEC
        CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
        LONDIR = 'W'
        IF (ISIGN .eq. -1) LONDIR = 'E'
        EHT = 0.D0
        call PREDV(ylat,ylon, eht, date1, iopt,
     1             jregn, vn, ve, vu)
        IF(JREGN .eq. 0) THEN     
          WRITE(I2,291)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,LOND,
     1       LOND,LONM,SLON,LONDIR,BLAB
  291        FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,
     1              1X,I2,1X,F8.5,1X,A1,1X,A17)
        ELSE
          call NEWCOR( ylat, ylon, eht, min1, min2,
     1      ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
          WRITE(I2,295)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,
     1          LOND,LONM,SLON,LONDIR,DN,DE,DU
  295     FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1           1X,F8.5,1X,A1,1X,3F8.3)
        ENDIF
        GO TO 290
  296   WRITE(LUOUT,297)
  297   FORMAT(' ***********************************'/
     1  ' Displacements have been calculated for the specified'/
     2  ' grid.  If you wish to calculate additional displacements,'/
     3  ' please indicate how you will supply the coordinates.'/)
      ELSEIF(OPTION .EQ. '3') THEN
        VOPT = '0'
        WRITE(LUOUT,300)
  300   FORMAT(' Enter name of blue-book file  ')
        READ(LUIN,310,err=603,iostat=ios) NAMEBB
        if (ios /= 0) goto 603
  310   FORMAT(A30)
        OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  320   READ(I1,330,END=350,err=604,iostat=ios) CARD
        if (ios /= 0) goto 604
  330   FORMAT(A80)
        TYPE = CARD(7:10)
        IF(TYPE .EQ. '*80*') THEN
          READ(CARD,340,err=605,iostat=ios) NAME,LATD,LATM,SLAT,JN,
     1              LOND,LONM,SLON,JW
          if (ios /= 0) goto 605
  340     FORMAT(BZ,14X,A30,I2,I2,F7.5,A1,I3,I2,F7.5,A1)
          NAME24 = NAME(1:24)
C         IF(JN.EQ.'S' .OR. JW.EQ.'E')GO TO 320
          YLAT =(DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
          IF (JN .eq. 'S') YLAT = -YLAT
          YLON =(DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
          IF (JW .eq. 'E') YLON = -YLON
          EHT = 0.D0
          call PREDV( ylat, ylon, eht, date1, iopt,
     1        jregn, vn, ve, vu)
          IF(JREGN .EQ. 0) THEN      
            WRITE(I2,345)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,SLON,
     1                   JW,BLAB
  345       FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1             F8.5,1X,A1,1X,A17)
          ELSE
            call NEWCOR( ylat, ylon, eht, min1, min2,
     1        ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
            WRITE(I2,160)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,SLON,
     1            JW,DN,DE,DU
          ENDIF
        ENDIF
        GO TO 320
  350   CLOSE(I1,STATUS='KEEP')
        WRITE(LUOUT,360)
  360   FORMAT(' ************************************'/
     1  ' Displacements have been calculated for the specified'/
     2  ' blue-book file.  If you wish to calculate additional'/
     3  ' displacements, please indicate how you will supply'/
     4  ' the horizontal coordinates.'/)
      ELSEIF(OPTION .EQ. '4') THEN
        VOPT = '0'
        CALL GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)
        XLON = TWOPI - XLON
        I = -1
  400   I = I + 1
        S = XMIN + I*XINC
        IF(S .GT. XMAX) GO TO 430
        IF(S .LT. 0.0D0) THEN
          S1 = -S
          AZ = BAZ
        ELSE
          S1 = S
          AZ = FAZ
        ENDIF
        CALL DIRCT1(XLAT,XLON,YLAT,YLON,AZ,AZ1,S1)
        YLON = TWOPI - YLON
        CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
        LATDIR = 'N'
        IF (ISIGN .eq. -1) LATDIR = 'S'
        CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
        LONDIR = 'W'
        IF (ISIGN .eq. -1) LONDIR = 'E'
        EHT = 0.D0
        call PREDV( ylat, ylon, eht, date1, iopt,
     1       jregn, vn, ve, vu)
        IF(JREGN .eq. 0) THEN
          WRITE(I2,405)NAMEG,I,LATD,LATM,SLAT,LATDIR,
     1          LOND,LONM,SLON,LONDIR, BLAB
  405     FORMAT(A10,I4,11X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1           F8.5,1X,A1,1X,A17)  
        ELSE
          call NEWCOR( ylat, ylon, eht, min1, min2,
     1        ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
          WRITE(I2,410)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                 SLON,LONDIR,DN,DE,DU
  410     FORMAT(A10,I4,11X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1           1X,F8.5,1X,A1,1X,3F8.3)
        ENDIF
        GO TO 400
  430   WRITE(LUOUT,440)
  440   FORMAT(' ***************************************'/
     1  ' Displacements have been calculated for the specified'/
     2  ' line.  If you wish to calculate additional displacements,'/
     3  ' please indicate how you will supply the coordinates.'/)
      ELSEIF(OPTION .EQ. '5') THEN
        VOPT = '0'
        EHT = 0.0D0
        write(luout,450)
  450   format(' Enter name of batch file ')
        read(luin, 451,err=606,iostat=ios) NAMEIF
        if (ios /= 0) goto 606
  451   format(a30)
        open(I1,FILE=NAMEIF,STATUS='OLD')
  455   read(I1,'(a)',END=460,err=607,iostat=ios) record
        if (ios /= 0) goto 607
c       write(i2, 457) record
        call interprate_latlon_record (record,XLAT,XLON,name24)
c       write (i2,456) xlat, xlon
c       write(i2,457) record
c 457   format (a50)
c 456   format(1x,' xlat = ', f15.3, 'xlon = ', f15.3)
        YLAT = (XLAT*3600.D0)/RHOSEC
        YLON = (XLON*3600.D0)/RHOSEC
        CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
        IF (ISIGN .EQ. 1) THEN
          JN = 'N'
        ELSE
          JN = 'S'
        ENDIF
        CALL TODMSS(YLON, LOND,LONM,SLON,ISIGN)
        IF (ISIGN .EQ. 1) THEN
          JW = 'W'
        ELSE
          JW = 'E'
        ENDIF
        CALL PREDV(YLAT,YLON,EHT,DATE1,IOPT,
     1             JREGN,VN,VE,VU)
        IF (JREGN .EQ. 0) THEN
          Write(I2,345) NAME24,LATD,LATM,SLAT, JN,
     1                  LOND,LONM,SLON,JW, BLAB
        ELSE
          CALL NEWCOR (YLAT, YLON, EHT, MIN1, MIN2,
     1         YLATT, YLONT, EHTNEW, DN,DE,DU,VN,VE,VU)
          WRITE(I2, 160) NAME24, LATD, LATM, SLAT, JN,
     1          LOND, LONM, SLON, JW, DN, DE , DU
        ENDIF
        GO TO 455
  460   CLOSE(I1, STATUS = 'KEEP')
        write(LUOUT, 470)
  470   format(' ******************************************'/
     1  ' Displacements have been calculated for the specified'/
     1  ' file.  If you wish to calculate additional displacements,'/
     1  ' please indicate how you will supply the horizontal '/
     1  ' coordinates.'/)
      ELSE
        WRITE(LUOUT,500)
  500   FORMAT(' Improper entry--select again.  ')
      ENDIF
      GO TO 85
  510 CONTINUE
      CLOSE(I2, STATUS = 'KEEP')
      RETURN

  600 write (*,'(/)') 
      write (*,*) 'Wrong answer in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)') 
      write (*,*) 'Wrong file name in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)') 
      write (*,*) 'Wrong OPTION in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)') 
      write (*,*) 'Wrong bbname in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  604 write (*,'(/)') 
      write (*,*) 'Wrong CARD from bbfile in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  605 write (*,'(/)') 
      write (*,*) 'Wrong CARD80 from bbfile in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  606 write (*,'(/)') 
      write (*,*) 'Wrong batch file name in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  607 write (*,'(/)') 
      write (*,*) 'Wrong record from batch in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
********************************************************************
      SUBROUTINE VELOC

*** Compute velocities at specified locations

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 16)
      parameter (nrsrch = 0)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER   CARD*80
      CHARACTER   NAMEF*30,NAME*30,NAMEBB*30,PVFILE*30
      CHARACTER   NAME24*24
      CHARACTER   BLAB*17
      CHARACTER   NAMEG*10
      CHARACTER   TYPE*4
      CHARACTER   OPTION*1,JN*1,JW*1,PVOUT*1
      CHARACTER   LATDIR*1, LONDIR*1
      character   frame1*24
      character   record*120

*** temporary code to plot results ****************
      character TEMPNA

      TEMPNA = '        '
      DUMMY = 0.0D0
*** end of temporary code *************************

      BLAB = 'OUTSIDE OF REGION'

      WRITE(LUOUT,10)
   10 FORMAT(' Please enter name for the file to contain the'/
     1       ' predicted velocities.  ')
      READ(LUIN,20,err=700,iostat=ios) NAMEF
      if (ios /= 0) goto 700
   20 FORMAT(A30)
      OPEN(I2,FILE=NAMEF,STATUS='UNKNOWN')
      CALL HEADER

      if (nrsrch .eq. 1) then
         WRITE(LUOUT,21)
   21    FORMAT( 
     1      ' Do you want to create a file of point-velocity (PV)'/
     1      ' records for use with the DYNAPG software (y/n)?   ')
         READ(LUIN,'(A1)',err=701,iostat=ios) PVOUT
         if (ios /= 0) goto 701
         IF(PVOUT .EQ. 'Y' .OR. PVOUT .EQ. 'y') THEN
            PVOUT = 'Y'
            WRITE(LUOUT,22)
   22       FORMAT(' Enter name for the file that is to contain'/
     1          ' the PV records.  ')
            READ(LUIN,'(A30)',err=702,iostat=ios) PVFILE
            if (ios /= 0) goto 702
            OPEN(I3,FILE = PVFILE, STATUS = 'UNKNOWN')
         ENDIF
      else
         PVOUT = 'N'
      endif

*** Choosing reference system for velocities
 1000 WRITE(LUOUT,1001)
 1001 FORMAT(' **************************************************'/
     1   ' Select the reference frame to be used for specifying'/
     2   ' positions and velocities. '/
c    6   '  0...Positions in NAD 83 and velocities relative to'/
c    6   '         a specified point having a specified velocity.')
     6      )
      call MENU1(iopt, frame1)
      if (iopt .ge. 1 .and. iopt .le. numref) then
        WRITE(I2,1002) frame1
 1002   FORMAT('VELOCITIES IN MM/YR RELATIVE TO ', a24 /)
        RVN = 0.0D0
        RVE = 0.0D0
        RVU = 0.0D0
      ELSEIF(IOPT .EQ. 0) THEN
        WRITE(LUOUT,1010)
        CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1       LONDIR, NAME24, X,Y,Z,YLAT,YLON,EHT)
        CALL GETREG(YLAT,YLON,JREGN)
        IF(JREGN .EQ. 0 ) THEN
          WRITE(LUOUT,1008)
 1008     FORMAT(' **********************************'/
     1    ' Point is not in modeled region and it can not'/
     1    ' be used as the reference location.'/)
          GO TO 1000
        ENDIF
        WRITE(I2,1009)
 1009   FORMAT('VELOCITIES IN MM/YR IN A REFERENCE FRAME FOR ',
     1  'WHICH THE FIRST POINT LISTED'/
     2  'BELOW HAS THE VELOCITY AS LISTED'//
     4  '**************REFERENCE SITE')
 1010   FORMAT(' For the reference point...'/)
        WRITE(LUOUT,1030)
 1030   FORMAT(' Enter northward velocity of reference location'/
     1         ' in mm/yr.  ')
        READ(LUIN,*,err=703,iostat=ios) SVN
        if (ios /= 0) goto 703
        WRITE(LUOUT,1040)
 1040   FORMAT(' Enter eastward velocity of reference location'/
     1         ' in mm/yr.  ')
        READ(LUIN,*,err=703,iostat=ios) SVE
        if (ios /= 0) goto 703
        WRITE(LUOUT,1050)
 1050   FORMAT(' Enter upward velocity of the reference location'/
     1         ' in mm/yr.  ')
        READ(LUIN,*,err=703,iostat=ios) SVU 
        if (ios /= 0) goto 703
        CALL COMVEL(YLAT,YLON,JREGN,VN,VE,VU)
        RVN = SVN - VN
        RVE = SVE - VE
        RVU = SVU - VU
        GLON = -YLON
        CALL TOVXYZ(YLAT,GLON,SVN,SVE,SVU,SVX,SVY,SVZ)
        WRITE(I2,1060) NAME24,LATD,LATM,SLAT,LATDIR,SVN,LOND,LONM,
     1      SLON,LONDIR,SVE,EHT,SVU,X,SVX,Y,SVY,Z,SVZ
 1060   FORMAT(/10X,A24,/
     1'LATITUDE   = ',2I3,F9.5,1X,A1,'  NORTH VELOCITY =',F7.2,' mm/yr'/
     2'LONGITUDE  = ',2I3,F9.5,1X,A1,'  EAST VELOCITY  =',F7.2,' mm/yr'/
     3'ELLIPS. HT. = ',F10.3,' m', 6X,'UP VELOCITY    =',F7.2,' mm/yr'/
     4'X =',F13.3,' m',14X,            'X VELOCITY     =',F7.2,' mm/yr'/
     5'Y =',F13.3,' m',14X,            'Y VELOCITY     =',F7.2,' mm/yr'/
     6'Z =',F13.3,' m',14X,            'Z VELOCITY     =',F7.2,' mm/yr')
        WRITE(I2,1070)
 1070   FORMAT('******************************')
      ELSE
        write(luout, 1080)
 1080   format(' Improper selection -- try again.  ')
        GO TO 1000
      ENDIF

*** Choosing input format for locations where velocities are to be predicted
      WRITE(LUOUT,30)
   30 FORMAT(' ************************************************'/
     1   ' Velocities will be predicted at each point whose'/
     2   ' horizontal position is specified.  Please indicate'/
     4   ' how you wish to supply positions.'/)
   40 WRITE(LUOUT,50)
   50 FORMAT(
     1 '    0... No more points. Return to main menu.'/
     2 '    1... Individual points entered interactively.'/
     3 '    2... Points on a specified grid.'/
     4 '    3... The *80* records in a specified blue-book file.'/
     5 '    4... Points on a specified line.  '/
     6 '    5... Batch file of delimited records of form: '/
     7 '         LAT,LON,TEXT '/
     8 '         LAT = latitude in degrees (positive north/DBL PREC) '/
     9 '         LON = longitude in degrees (positive west/DBL PREC) '/
     1 '         TEXT = Descriptive text (CHARACTER*24) '/
     1 '         Example:  '/
     2 '         40.731671553,112.212671753,SALT AIR '/)

      READ(LUIN,60,err=704,iostat=ios) OPTION
      if (ios /= 0) goto 704
   60 FORMAT(A1)

      IF(OPTION .EQ. '0') THEN
        GO TO 610               
      ELSEIF(OPTION .EQ. '1') THEN
        CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1         LONDIR, NAME24, X,Y,Z,YLAT,YLON,EHT)
        CALL GTOVEL(YLAT,YLON,EHT,RVN,RVE,RVU,
     1              VN,VE,VU,VX,VY,VZ,JREGN,IOPT)
        IF(JREGN .eq. 0 ) THEN
          WRITE(LUOUT,80)
   80     FORMAT(' *************************************'/
     1    ' A velocity can not be predicted because'/
     1    ' the point is outside of the modeled region.'/
     2    ' For additional velocities, please indicate how'/
     3    ' you wish to supply the horizontal coordinates.'/)
        ELSE
          WRITE(LUOUT,90) VN, VE, VU, VX, VY, VZ
   90     FORMAT(' **************************************'/
     1           ' Northward velocity = ',F6.2,' mm/yr'/
     1           ' Eastward velocity  = ',F6.2,' mm/yr'/
     1           ' Upward velocity    = ',F6.2,' mm/yr'//
     1           ' X-dim. velocity    = ',F6.2,' mm/yr'/
     1           ' Y-dim. velocity    = ',F6.2,' mm/yr'/
     1           ' Z-dim. velocity    = ',F6.2,' mm/yr'/
     1           ' **************************************'//
     2         ' For additional velocities, please indicate how'/
     3         ' you wish to specify the horizontal coordinates.'/)
          WRITE(I2,1060)NAME24,LATD,LATM,SLAT,LATDIR,VN,LOND,LONM, 
     1              SLON, LONDIR, VE, EHT,VU,X,VX,Y,VY,Z,VZ
  100     FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,F8.5,
     1        1X,A1,1X,3F8.2)
          IF(PVOUT .EQ. 'Y') THEN
            CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
          ENDIF
        ENDIF
      ELSEIF(OPTION .EQ. '2') THEN
           EHT = 0.D0
           WRITE(I2,105)
  105      FORMAT(
     1     'NAME OF SITE',14X,'LATITUDE          LONGITUDE    ',
     2     '       NORTH    EAST    UP'/)
  106      FORMAT(
     1     'NAME OF LINE', 6X,'LATITUDE          LONGITUDE    ',
     2     '       NORTH    EAST    UP'/)
       CALL GETGRD(NAMEG,MINLAT,MAXLAT,IDS,MINLON,MAXLON,JDS)
       I = -1
  110  I = I + 1
       LAT = MINLAT + I*IDS
       IF(LAT .GT. MAXLAT) GO TO 150
       XLAT = DBLE(LAT)/RHOSEC
       CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
       LATDIR = 'N'
       IF (ISIGN .eq. -1) LATDIR = 'S'
       J = -1
  120  J = J + 1
       YLAT = XLAT
       ON = MINLON + J*JDS
       IF(LON .GT. MAXLON) GO TO 110
       YLON = DBLE(LON)/RHOSEC
       CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
       LONDIR = 'W'
       IF (ISIGN .eq. -1) LONDIR = 'E'
       CALL GTOVEL(YLAT,YLON,EHT,RVN,RVE,RVU,
     1             VN,VE,VU,VX,VY,VZ,JREGN,IOPT)
       IF(JREGN .EQ. 0 ) THEN
         WRITE(I2,125)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                SLON,LONDIR,BLAB
  125    FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1                1X,F8.5,1X,A1, 1X,A17)
       ELSE

*** temporary code to plot results
         ZLAT = DBLE(LATD) + DBLE(LATM)/60.D0 + SLAT/3600.D0
         ZLON = DBLE(LOND) + DBLE(LONM)/60.D0 + SLON/3600.D0
         ZLON = 360.D0 - ZLON
         TOTVEL = DSQRT(VN*VN + VE*VE)

*** code to create vectors to be plotted
c	     IF (TOTVEL .GE. 5.0D0 .AND. TOTVEL .LE. 50.D0) THEN
c	       TVE = VE/10.D0
c	       TVN = VN/10.D0
c	       WRITE(I2,129) ZLON,ZLAT,TVE,TVN,DUMMY,DUMMY,DUMMY,
c    1                   TEMPNA
c 129      FORMAT(F10.6, 2x, F9.6, 2X, 5(F7.2, 2x), A8)
c 	     ENDIF

*** code to create contour plots
C	     WRITE(I2,129) ZLON, ZLAT, TOTVEL
C 129    FORMAT(F6.2, 1X, F6.2, 1X, F5.1)
 
*** end of temporary code

*** beginning of original code
         WRITE(I2,130)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,
     1              LOND,LONM,SLON,LONDIR,VN,VE,VU
  130    FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,
     1          I3,1X,I2,1X,F8.5,1X,A1,1X,3F8.2)
*** end of original code

         IF(PVOUT .EQ. 'Y') THEN
           CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
         ENDIF
       ENDIF
       GO TO 120
  150      WRITE(LUOUT,160)
  160      FORMAT(' **********************************************'/
     1     ' Velocities have been calculated for the specified grid.'/
     2     ' If you wish to calculate additional velocities, please'/
     3     ' indicate how you will specify positional coordinates.'/)
      ELSEIF(OPTION .EQ. '3') THEN
        EHT = 0.D0
        WRITE(I2,105)
        WRITE(LUOUT,200)
  200   FORMAT(' Enter name of blue-book file.  ')
        READ(LUIN,210,err=705,iostat=ios) NAMEBB
        if (ios /= 0) goto 705
  210   FORMAT(A30)
        OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  220   READ(I1,230,END=250,err=706,iostat=ios) CARD
        if (ios /= 0) goto 706
  230   FORMAT(A80)
        TYPE = CARD(7:10)
        IF(TYPE .EQ. '*80*') THEN
          READ(CARD,240,err=707,iostat=ios)NAME,LATD,LATM,SLAT,JN,
     1                                     LOND,LONM,SLON,JW
          if (ios /= 0) goto 707
  240     FORMAT(BZ,14X,A30,I2,I2,F7.5,A1,I3,I2,F7.5,A1)
          NAME24 = NAME(1:24)
          YLAT =(DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
          YLON =(DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
          if (jn .eq. 'S') ylat = -ylat
          if (jw .eq. 'E') ylon = -ylon
          CALL GTOVEL(YLAT,YLON,EHT,RVN,RVE,RVU,
     1                VN,VE,VU,VX,VY,VZ,JREGN,IOPT)
          IF(JREGN .EQ. 0 ) THEN
            WRITE(I2,245)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,
     1                   SLON,JW,BLAB
  245       FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1             1X,F8.5,1X,A1,1X,A17)
          ELSE
            WRITE(I2,100)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,
     1                   SLON,JW,VN,VE,VU  
            IF(PVOUT .EQ. 'Y') THEN
              CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
            ENDIF
          ENDIF
        ENDIF
        GO TO 220
  250   CLOSE(I1,STATUS='KEEP')
        WRITE(LUOUT,260)
  260      FORMAT(' ****************************************'/
     1     ' Velocities have been calculated for the specified '/
     2     ' blue-book file.  If you wish to calculate additional'/
     3     ' velocities, please indicate how you will supply the'/
     4     ' horizontal coordinates.'/)
      ELSEIF(OPTION .EQ. '4') THEN
        EHT = 0.D0
        WRITE(I2,106)
        CALL GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)
        XLON = TWOPI - XLON
        I = -1
  300   I= I+1
        S = XMIN + I*XINC
        IF(S .GT. XMAX) GO TO 330
        IF(S .LT. 0.0D0) THEN
          S1 = -S
          AZ = BAZ
        ELSE
          S1 = S
          AZ = FAZ
        ENDIF
        CALL DIRCT1(XLAT,XLON,YLAT,YLON,AZ,AZ1,S1)
        YLON = TWOPI - YLON
        CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
        LATDIR = 'N'
        IF (ISIGN .eq. -1) LATDIR = 'S'
        CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
        LONDIR = 'W'
        IF (ISIGN .eq. -1) LONDIR = 'E'
        CALL GTOVEL(YLAT,YLON,EHT,RVN,RVE,RVU,
     1              VN,VE,VU,VX,VY,VZ,JREGN,IOPT)
        IF(JREGN .EQ. 0 ) THEN
          WRITE(I2,305)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                 SLON,LONDIR,BLAB
  305     FORMAT(A10,I4,3X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,
     1           I2,1X,F8.5,1X,A1,1X,A17)
        ELSE
          WRITE(I2,310)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1                 LONDIR,VN,VE,VU
  310     FORMAT(A10,I4,3X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1           F8.5,1X,A1,1X,3F8.2)
          IF(PVOUT .EQ. 'Y') THEN
            CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
          ENDIF
        ENDIF
        GO TO 300
  330   WRITE(LUOUT,340)
  340   FORMAT(' ***********************************************'/
     1  ' Velocities have been calculated for the specified line.'/
     2  ' If you wish to calculate additional velocities, please'/
     3  ' indicate how you will specify positional coordinates.'/)
      ELSEIF(OPTION .EQ. '5') THEN
           EHT = 0.0D0
           WRITE(I2, 105)
           WRITE(LUOUT,400)
  400      FORMAT(' Enter name of input file.  ')
           READ(LUIN,410,err=708,iostat=ios) NAMEBB
           if (ios /= 0) goto 708
  410      FORMAT(A30)
           OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  420      READ(I1,'(a)',END=450,err=709,iostat=ios) record
           call interprate_latlon_record (record,XLAT,XLON,NAME24)
           if (ios /= 0) goto 709
           YLAT = (XLAT*3600.D0)/RHOSEC
           YLON = XLON*3600.d0/RHOSEC
           CALL TODMSS(YLAT, LATD, LATM, SLAT, ISIGN)
           IF (ISIGN .eq. 1) then
              JN = 'N'
           Else
              JN = 'S'
           endif
           call TODMSS(YLON, LOND, LONM, SLON, ISIGN)
           IF (ISIGN .eq. 1) then
              JW = 'W'
           else
              JW = 'E'
           endif
           CALL GTOVEL(YLAT,YLON,EHT,RVN,RVE, RVU,
     1                 VN,VE,VU,VX,VY,VZ,JREGN,IOPT)
           IF (JREGN .eq. 0) then
              write(I2,245)NAME24,LATD,LATM,SLAT,JN,
     1                     LOND,LONM,SLON,JW,BLAB
           ELSE
              write(I2,100)NAME24,LATD,LATM,SLAT,JN,
     1                      LOND,LONM,SLON,JW,VN,VE,VU
              IF(PVOUT .EQ. 'Y') then
                call PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
              ENDIF
           ENDIF
           GO TO 420
  450      CLOSE(I1, STATUS = 'KEEP')
           write(LUOUT, 460)
  460      format(' ********************************'/
     1     ' Velocities have been calculated for the specified'/
     2     ' file.  if you wish to calculate additional velocities, '/
     3     ' please indicate how you will supply the coordinates.'/)
      ELSE
           WRITE(LUOUT,600)
  600      FORMAT(' Improper entry--select again.  ')
      ENDIF
      GO TO 40
  610 CONTINUE
      CLOSE(I2, STATUS = 'KEEP')
      IF( PVOUT .EQ. 'Y') CLOSE(I3, STATUS = 'KEEP')
      RETURN

  700 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  701 write (*,'(/)') 
      write (*,*) 'Failed to read answer: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  702 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  703 write (*,'(/)') 
      write (*,*) 'Failed to read velocity: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  704 write (*,'(/)') 
      write (*,*) 'Failed to read OPTION: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  705 write (*,'(/)') 
      write (*,*) 'Failed to read name of bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  706 write (*,'(/)') 
      write (*,*) 'Failed to read card from bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  707 write (*,'(/)') 
      write (*,*) 'Failed to read card80 from bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  708 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  709 write (*,'(/)') 
      write (*,*) 'Failed to read bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop


      END
**************************************************************************
      SUBROUTINE GTOVEL(YLAT,YLON,EHT,RVN,RVE,RVU,
     1               VN,VE,VU,VX,VY,VZ,JREGN,IOPT)

*** Compute velocity in appropriate reference frame for point with
*** latitude YLAT (radians), longitude YLON (radians, positive west)
*** and ellipsoid height EHT (meters) in this reference frame.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      logical  Is_iopt_NAD83

      Is_iopt_NAD83 = (IOPT == 1)

*** Get reference latitude RLAT and reference longitude RLON in NAD 83
c     IF(IOPT .EQ. 0 .OR. IOPT .EQ. 1) THEN   !Not NAD83 antmore
      IF(IOPT .EQ. 15) THEN
          RLAT = YLAT
          RLON = YLON
      ELSE
          ELON = -YLON
          CALL TOXYZ(YLAT,ELON,EHT,X,Y,Z)
          DATE = 2010.0d0
c         CALL XTONAD(X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !Removed on 09/12/2014. The velocity grids are not converted to NAD83 anymore.
          CALL XTO08 (X,Y,Z,RLAT,RLON,EHT08,DATE,IOPT)   !They are in ITRF2008
      ENDIF

*** Get velocity in NAD 83  !Not anymore since 09/12/2014. Now it is in ITRF2008
      CALL GETREG(RLAT,RLON,JREGN)
      IF (JREGN .EQ. 0) RETURN
      CALL COMVEL(RLAT,RLON,JREGN,VN,VE,VU)

      VN = VN + RVN
      VE = VE + RVE
      VU = VU + RVU
      ELON = -YLON
      CALL TOVXYZ(YLAT,ELON,VN,VE,VU,VX,VY,VZ)

*** Convert velocity into another reference frame if needed
c     IF(IOPT .NE. 0 .AND. IOPT .NE. 1) THEN    !Commented out on 09/12/2014
      IF(IOPT .NE. 15) THEN
        IF(Is_iopt_NAD83) THEN
c         CALL VTRANF(X,Y,Z,VX,VY,VZ, 1, IOPT)  !No longer NAD83
          CALL VTRANF(X,Y,Z,VX,VY,VZ, 15, IOPT)  !No longer NAD83
        else
          CALL VTRANF_IERS(X,Y,Z,VX,VY,VZ, 15, IOPT) !Now ITRF2008
        endif
          CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)
      ENDIF
     
      RETURN
      END
****************************************************************************
      SUBROUTINE XTO08 (X,Y,Z,RLAT,WLON,EHT08,DATE,IOPT)

*** Converts X,Y,Z in specified datum to latitude and
*** longitude (in radians) and height (meters) in ITRF2008
*** datum with longitude positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL FRMXYZ
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

*** Convert to cartesian coordinates in ITRF2008   
c     if (iopt .eq. 0 .or. iopt .eq. 1) then
      if (iopt .eq. 15) then
         x2 = x
         y2 = y
         z2 = z
      elseif (iopt .eq. 1) then
         call toit94(x,y,z,x1,y1,z1,date,iopt)
         call frit94(x1,y1,z1,x2,y2,z2,date,15)
      else 
         call toit94_iers(x,y,z,x1,y1,z1,date,iopt)
         call frit94_iers(x1,y1,z1,x2,y2,z2,date,15)
      endif

***Convert to geodetic coordinates
      IF(.NOT.FRMXYZ(X2,Y2,Z2,RLAT,ELON,EHT08))STOP 666

      WLON = -ELON
 100  IF(WLON .LT. 0.D0) THEN
          WLON = WLON + TWOPI
          GO TO 100
      ENDIF

c     write (*,*) "FROM XT08 ",RLAT*180.d0/pi,WLON*180.d0/pi,EHT08
      RETURN
      END
 
****************************************************************************
      SUBROUTINE XTONAD (X,Y,Z,RLAT,WLON,EHTNAD,DATE,IOPT)

*** Converts X,Y,Z in specified datum to latitude and
*** longitude (in radians) and height (meters) in NAD 83 
*** datum with longitude positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL FRMXYZ
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

*** Convert to cartesian coordinates in NAD 83
      if (iopt .eq. 0 .or. iopt .eq. 1) then
         x2 = x
         y2 = y
         z2 = z
      else
         call toit94(x,y,z,x1,y1,z1,date,iopt)
         jopt = 1
         call frit94(x1,y1,z1,x2,y2,z2,date,jopt)
      endif

***Convert to geodetic coordinates
      IF(.NOT.FRMXYZ(X2,Y2,Z2,RLAT,ELON,EHTNAD))STOP 666

      WLON = -ELON
 100  IF(WLON .LT. 0.D0) THEN
          WLON = WLON + TWOPI
          GO TO 100
      ENDIF
      RETURN
      END
 
******************************************************
      subroutine TRFPOS

*** Transform position between reference frames

      implicit double precision (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)
      parameter (nbbdim = 10000)
      parameter (rad2deg = 180.d0/3.14159265358979d0)

      double precision  lat,lon
      character    card*80,namebb*80,nameif*80,name24*80
      character    record*120                    
      character    namef*30
      character    frame1*24, frame2*24
      character    jn*1,jw*1,LATDIR*1,LONDIR*1
      character    option*1, answer*1, vopt*1
      character    PID*6,PIDs*6
      character    HTDP_version*10
      LOGICAL      FRMXYZ
      LOGICAL      TEST
      LOGICAL      Is_inp_NAD83,Is_out_NAD83
      LOGICAL      Is_inp_NAD83PAC,Is_out_NAD83PAC
      LOGICAL      Is_inp_NAD83MAR,Is_out_NAD83MAR
      logical      am_I_in_or_near_AK
      logical      am_I_in_or_near_AS
      logical      am_I_in_or_near_CONUS
      logical      am_I_in_or_near_CQ  
      logical      am_I_in_or_near_Guam
      logical      am_I_in_or_near_HI  
      logical      am_I_in_or_near_PR  
      logical      am_I_in_or_near_VQ  
      logical      am_I_in_or_near_KW  

      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /ARRAYS/ HT(nbbdim), LOC(nbbdim),PIDs(nbbdim)
      COMMON /VERSION/ HTDP_version

      write(luout,80)
   80 format(
     1  ' Please enter the name for the file to contain'/
     2  ' the transformed positions.  ')
      read(luin,'(a30)',err=600,iostat=ios) namef
      if (ios /= 0) goto 600
      open(i2, file = namef, status = 'unknown')
      CALL HEADER

   95 write(luout,100)
  100 format (' *******************************************'/
     1  ' Enter the reference frame of the input positions')
      call MENU1(iopt1, frame1)
      if (iopt1 .lt. 1 .or. iopt1 .gt. numref) then
         write(luout,*) ' Improper selection -- try again.  '
         go to 95
      endif
      Is_inp_NAD83    = (iopt1 == 1 .and. frame1(1:3) /= "WGS") 
      Is_inp_NAD83PAC = (iopt1 == 12) 
      Is_inp_NAD83MAR = (iopt1 == 13) 

  105 write(luout,110)
  110 format (/' Enter the reference frame for the output positions')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
         write(luout,*) ' Improper selection -- try again.  '
         go to 105
      endif
      Is_out_NAD83    = (iopt2 == 1 .and. frame1(1:3) /= "WGS")
      Is_out_NAD83PAC = (iopt2 == 12)
      Is_out_NAD83MAR = (iopt2 == 13)

      write(luout,120)
  120 format (/
     1   ' Enter the reference date of the input positions.')
  121 CALL GETMDY(MONTH1, IDAY1, IYEAR1, DATE1, MIN1, TEST)
      IF(TEST) then
         write(luout,*) ' Do you wish to re-enter the date? (y/n)'
         read(luin, '(A1)',err=601,iostat=ios) ANSWER
         if (ios /= 0) goto 601
         if(ANSWER .eq. 'y' .or. ANSWER .eq. 'Y')GO TO 121
         RETURN
      ENDIF

      write(luout,130)
  130 format(/
     1  '  Enter the reference date for the output positions.  ')
  131 CALL GETMDY(MONTH2, IDAY2, IYEAR2, DATE2, MIN2, TEST)
      IF(TEST) then
        write(luout,*) ' Do you wish to re-enter the date? (y/n)'
        read(luin, '(A1)',err=601,iostat=ios) ANSWER
        if (ios /= 0) goto 601
        if(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y') GO TO 131
        RETURN
      ENDIF

      write(i2, 150) frame1, month1, iday1, iyear1, date1,
     1               frame2, month2, iday2, iyear2, date2
  150 format(' TRANSFORMING POSITIONS FROM ',A24,' (EPOCH = ',
     1  I2.2,'-',I2.2,'-',I4,' (',F9.4,'))'/26X,
     1  'TO ', A24, ' (EPOCH = ',I2.2,'-',I2.2,'-',I4,' (',F9.4,'))'/)
c    1  13X,'INPUT COORDINATES   OUTPUT COORDINATES',
c    1     4x, 'INPUT VELOCITY' /)
 
  160 write(luout, 170)
  170 format (/' ***************************************'/
     1  ' Coordinates will be transformed at each specified point.'/
     1  ' Please indicate how you wish to supply positions.'/
     1  '    0... No more points.   Return to main menu.'/
     1  '    1... Individual points entered interactively.'/
     1  '    2... The *80* records in a specified blue-book file.  '/
     1  '    3... Transform positions contained in batch file '/
     1  '         of delimited records of the form: '/
     1  '         LAT,LON,EHT,TEXT' /
     1  '         LAT = latitude in degrees (positive north/DBL PREC)'/
     1  '         LON = longitude in degrees (positive west/DBL PREC)'/
     1  '         EHT = ellipsoid height in meters (DBL PREC)'/
     1  '         TEXT = Descriptive text (CHARACTER*24) '/
     1  '         Example:   '/
     1  '         40.731671553,112.212671753,34.241,SALT AIR  '/ 
     1  '    4... Transform positions contained in batch file '/
     1  '         of delimited records of the form: '/
     1  '         X, Y, Z,TEXT' /
     1  '         X, Y, Z are Cartesian coordinates of the station'/
     1  '         TEXT = Descriptive text (CHARACTER*24) '/
     1  '         Example:   '/
     1  '           -86682.104,-5394026.861,3391189.647,SALT AIR  '/
     1  '         Necessary velocities predicted by HTDP          '/
     1  '    5... Transform positions using velocities contained in a '/
     1  '         batch file of delimited records of the form: '/
     1  '         X, Y, Z, Vx, Vy, Vz, TEXT' /
     1  '         X, Y, Z are Cartesian coordinates of the station'/
     1  '         Vx,Vy,Vz are Cartesian velocities in m/year     '/
     1  '         TEXT = Descriptive text (CHARACTER*24) '/
     1  '         Example:   '/
     1'-86682.104,-5394026.861,3391189.647,-0.012,-0.001,-0.001,SALT'/)

      read(luin,'(a1)',err=602,iostat=ios) option
      if (ios /= 0) goto 602

c     write (*,*) "Starting to look at options"

      if (option .eq. '0') then
         go to 500
      elseif (option .eq. '1') then
        vxsave = 0.d0
        vysave = 0.d0
        vzsave = 0.d0
        vnsave = 0.d0
        vesave = 0.d0
        vusave = 0.d0
  180   call GETPNT(latd, latm, slat, LATDIR, lond, lonm, slon, 
     1        LONDIR, name24, x, y, z, ylat, ylon, eht)
        
C  Added by JS on 07/22/2015 to limit NAD83 to the US*****************************************************
C  This was reversed on 05/11/2017                   *****************************************************

c        if (Is_inp_NAD83 .or. Is_out_NAD83) then
c          lat = ylat*rad2deg ; lon = 360.d0 - ylon*rad2deg    
c          if (lon < 0.d0) lon = lon + 360.d0
c          if(.not. am_I_in_or_near_AK  (lat,lon) .and.
c    &       .not. am_I_in_or_near_CONUS(lat,lon) .and.
c    &       .not. am_I_in_or_near_PR   (lat,lon) .and.
c    &       .not. am_I_in_or_near_VQ   (lat,lon)) then 
c            write (luout,'(10x,a,8x,a)') 
c    &       "NAD83 (2011) is defined only in US territories",name24
c            stop
c          endif
c        endif

c       if(.not. am_I_in_or_near_AK  (lat,lon) .and.
c    &     .not. am_I_in_or_near_AS   (lat,lon) .and.
c    &     .not. am_I_in_or_near_CONUS(lat,lon) .and.
c    &     .not. am_I_in_or_near_CQ   (lat,lon) .and.
c    &     .not. am_I_in_or_near_Guam (lat,lon) .and.
c    &     .not. am_I_in_or_near_HI   (lat,lon) .and.
c    &     .not. am_I_in_or_near_PR   (lat,lon) .and.
c    &     .not. am_I_in_or_near_VQ   (lat,lon) .and.
c    &     .not. am_I_in_or_near_KW   (lat,lon)) then
c          write (luout,'(10x,a,8x,a)') 
c    &     "FYI, NAD83 (2011) is defined only in US territories",name24
c          stop
c       endif

c        if (Is_inp_NAD83MAR .or. Is_out_NAD83MAR) then
c          lat = ylat*rad2deg ; lon = 360.d0 - ylon*rad2deg    
c          if (lon < 0.d0) lon = lon + 360.d0
c          if(.not. am_I_in_or_near_GUAM (lat,lon) .and.
c    &        .not. am_I_in_or_near_CQ   (lat,lon)) then 
c            write (luout,'(10x,a,8x,a)') 
c    &       "NAD83 (2011) is defined only in US territories",name24
c            stop
c          endif
c        endif

C  Until here on 07/22/2015 to limit NAD83 to the US******************************************************

C  Added by JS on 09/10/2014

        if (Is_inp_NAD83 .or. Is_out_NAD83) then
          call TOIT94(x, y, z, x1, y1, z1, date1, iopt1)
          call FRIT94(x1, y1, z1, x2, y2, z2, date1, iopt2)
        else
          call TOIT94_IERS(x, y, z, x1, y1, z1, date1, iopt1)
          call FRIT94_IERS(x1, y1, z1, x2, y2, z2, date1, iopt2)
        endif
        if (min1 .ne. min2) then                                  !i.e., if the input and output dates are different
          if(.not.FRMXYZ(x2,y2,z2,ylat1,elon1,eht1)) STOP 666
          ylon1 = -elon1
          if(ylon1 .lt. 0.d0) ylon1 = ylon1 + twopi
c         write (*,*) ylon1*180.d0/pi,ylat1*180.d0/pi,eht1
          call GETVLY(ylat1,elon1,vx,vy,vz,vn,ve,vu,vopt,210)
c         write (*,*) vx,vy,vz,vn,ve,vu
          if (vopt .eq. '0') then
            call PREDV(ylat1,ylon1,eht1,date1,iopt1,jregn,vn,ve,vu)
c           write (*,*) jregn,vn,ve,vu
            if (jregn .eq. 0) then
              write(luout, 200)
  200         format(/' ****************************************'/
     1  ' These coordinates can not be transformed to a different'/
     1  ' date because the point is outside of the modeled region.'/)
              go to 220
            else
              call TOVXYZ( ylat1, elon1, vn, ve, vu, vx, vy, vz)
c             write (*,*) vx,vy,vz,vn,ve,vu
            endif
          endif
          vxsave = vx
          vysave = vy
          vzsave = vz
          vnsave = vn
          vesave = ve
          vusave = vu
          if (Is_inp_NAD83 .or. Is_out_NAD83) then
            call VTRANF( x, y, z, vx, vy, vz, iopt1, iopt2)
          else
            call VTRANF_IERS( x, y, z, vx, vy, vz, iopt1, iopt2)
          endif
          call TOVNEU( ylat1, elon1, vx, vy, vz, vn, ve, vu)
c         write (*,*) "Passed TOVNEU",vn, ve, vu
          call NEWCOR( ylat1, ylon1, eht1, min1, min2, 
     1            ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
c         write (*,*) "Passed NEWCOR",vn, ve, vu
          call TOXYZ(ylatt,-ylont, ehtnew, xt, yt, zt)
c         write (*,*) "Passed TOXYZ"
        else
          xt = x2
          yt = y2
          zt = z2
          if(.not.(FRMXYZ(xt,yt,zt,ylatt,elont,ehtnew))) STOP 666
          ylont = -elont
          if(ylont .lt. 0.0d0) ylont = ylont + twopi
c         write (*,*) ylont*180.d0/pi,ylatt*180.d0/pi,ehtnew
        endif
         
        call PRNTTP(x, y, z, xt, yt, zt, ylat, ylatt,
     1      ylon, ylont, eht, ehtnew, name24, 1,
     2      vxsave, vysave, vzsave, vnsave, vesave, vusave)

  220   write(luout,*) ' Transform more positions? (y/n)  '
        read(luin,'(a1)',err=601,iostat=ios) answer
        if (ios /= 0) goto 601
        if(answer .eq. 'Y' .or. answer .eq. 'y') go to 180
        close (i2, status = 'keep')
      elseif (option .eq. '2') then     !Blue book input and output
         vxsave = 0.d0
         vysave = 0.d0
         vzsave = 0.d0
         vnsave = 0.d0
         vesave = 0.d0
         vusave = 0.d0
         write(luout, 300)
  300    format(/ 
     1   ' Enter name of the blue-book file:  ')
         read(luin, '(a)',err=603,iostat=ios) namebb
         if (ios /= 0) goto 603
         open (i1, file = namebb, status = 'old')

*** Obtaining the ellipsoid heights from the blue book file

         do i = 1, nbbdim
           ht(i) = 0.d0
         enddo    

  302    read (i1, '(a80)', end = 309,err=604,iostat=ios) card
         if (ios /= 0) goto 604

         if (card(7:10) .eq. '*80*') then
c           read (card, 303) isn, eht
c 303       format (bz, 10x, i4, t70, f6.2)
            read (card, 303,err=605,iostat=ios) PID,isn,eht
            if (ios /= 0) goto 605
  303       format (a6,4x,i4, t70, f6.2)
            if (ht(isn) .eq. 0.d0) ht(isn) = eht
            PIDs(isn) = PID

         elseif (card(7:10) .eq. '*86*') then
           if (card(46:52) .ne. '       ') then
             read (card, 304,err=606,iostat=ios) isn, eht
             if (ios /= 0) goto 606
  304        format (bz, 10x, i4, t46, f7.3)
             ht(isn) = eht
           else
               read (card, 305,err=606,iostat=ios) isn, oht, ght
               if (ios /= 0) goto 606
  305          format (bz, 10x, i4, t17, f7.3, t36, f7.3)
               ht(isn) = oht + ght
           endif
         endif
         go to 302
  309    rewind i1
c        open (i6,file='transformed_'//namebb)

C  write some comments in the transformed files

         call extract_name (namebb,iii)
c        write (I2,1309) namebb(1:iii)
c        write (I2,1309) trim(namebb)
 1309    format (/'The file, transformed_',a,
     &   ', contains the transformed input file.'/)

         write (i2,1310) HTDP_version
 1310    format (' ***CAUTION: This file was processed using HTDP',
     &           ' version ',a10, '***')
         write (i2,1311) frame2       
 1311    format (' ***CAUTION: Coordinates in this file are in ',    
     &           a24, '***')
         WRITE (i2, 1312) MONTH2, IDAY2, IYEAR2, DATE2
 1312    FORMAT(' ***CAUTION: Coordinates in this file have been ',
     *   'updated to ',I2,'-',I2.2,'-',I4, '=(',F8.3,') ***'/)

*** Reread blue book file to get latitudes and longitudes
*** and to perform transformations
*** and to write transformed lat and lon to the *80* records of the output bfile
*** and to write transformed ellip. h to the *86* records of the output bfile

  310    read(i1, '(a80)', end = 390,err=604,iostat=ios) card
         if (ios /= 0) goto 604
         if (card(7:10) .eq. '*80*') then
c          write (*,'(a80)'                              ) card
           read(card,320,err=605,iostat=ios) isn,name24(1:30),latd,
     &                           latm,slat,jn,lond, lonm, slon, jw
c          write (*,321) latd, latm, slat, jn, lond, lonm, slon, jw

           if (ios /= 0) goto 605
  320      format(10x, i4, a30, i2, i2, f7.5, a1,
     &            i3, i2, f7.5, a1)
  321      format(i2, i2, f7.5, a1,3x, i3, i2, f7.5, a1)
           ylat = (dble((latd*60 + latm)*60) + slat) / rhosec
           ylon = (dble((lond*60 + lonm)*60) + slon) / rhosec
           if (jn .eq. 'S') ylat = -ylat
           if (jw .eq. 'E') ylon = -ylon
           eht = ht(isn)
           call TOXYZ(ylat, -ylon, eht, x, y, z)

C  Modified in 09/15/2014 by JS to accommodate the IERS 96-97 trans. par.
           if (Is_inp_NAD83 .or. Is_out_NAD83) then
             call TOIT94(x, y, z, x1, y1, z1, date1, iopt1)
             call FRIT94(x1, y1, z1, x2, y2, z2, date1, iopt2)
           else
             call TOIT94_IERS(x, y, z, x1, y1, z1, date1, iopt1)
             call FRIT94_IERS(x1, y1, z1, x2, y2, z2, date1, iopt2)
           endif
C  End changes

           if (min1 .ne. min2) then
              if(.not.FRMXYZ(x2,y2,z2,ylat1,elon1,eht1)) STOP 666
              ylon1 = -elon1
              if (ylon1 .lt. 0.0d0) ylon1 = ylon1 + twopi
              call PREDV(ylat1, ylon1, eht1, date1, iopt1,
     1                   jregn, vn, ve, vu)
              if (jregn .eq. 0) then
                write(i2, 330) name24
  330           format (/ 1x,a,
     1              ' is outside of the modeled region.'/)
                go to 310
              else
                call TOVXYZ(ylat1,elon1,vn,ve,vu,vx,vy,vz)
                vxsave = vx
                vysave = vy
                vzsave = vz
                vnsave = vn
                vesave = ve
                vusave = vu

C  Modified on 09/15/2014 by JS to accomodate the IERS trans. par.
                 if (Is_inp_NAD83 .or. Is_out_NAD83) then
                   call VTRANF(x,y,z,vx,vy,vz,iopt1,iopt2)
                 else
                   call VTRANF_IERS( x, y, z, vx, vy, vz, iopt1, iopt2)
                 endif
C  End changes
                 call TOVNEU(ylat1,elon1,vx,vy,vz,vn,ve,vu)
                 call NEWCOR(ylat1,ylon1,eht1,min1,min2,
     1                       ylatt,ylont,ehtnew,dn,de,dui,vn,ve,vu)
                 call TOXYZ(ylatt,-ylont,ehtnew,xt,yt,zt)
               endif
            else
               xt = x2
               yt = y2
               zt = z2
               if(.not.FRMXYZ(xt,yt,zt,ylatt,elont,ehtnew))STOP 666
               ylont = -elont
               if (ylont .lt. 0.0d0) ylont = ylont + twopi
            endif
c           call PRNTTP(x,y,z, xt,yt, zt, ylat, ylatt,
c    1         ylon, ylont, eht, ehtnew, name24, 0,
c    1         vxsave, vysave, vzsave, vnsave, vesave, vusave)

C  Now write output *80* records to the output bfile

            call TODMSS (ylatt,latdeg,latmin,seclat,isign)
            latsec = nint(100000*seclat)
            latdir = 'N'
            if (isign .eq. -1) latdir = 'S'
            call TODMSS (ylont,londeg,lonmin,seclon,isign)
            lonsec = nint(100000*seclon)
            londir = 'W'
            if (isign .eq. -1) londir = 'E'
            write (card(45:46),'(i2)') latdeg
            write (card(47:48),'(i2)') latmin
            write (card(49:55),'(i7)') latsec
            do i=45,55
              if (card(i:i) == ' ') card(i:i) = '0'
            enddo
            write (card(56:56),'(a1)') latdir
 
            write (card(57:59),'(i3)') londeg
            write (card(60:61),'(i2)') lonmin
            write (card(62:68),'(i7)') lonsec
            do i=57,68
              if (card(i:i) == ' ') card(i:i) = '0'
            enddo
            write (card(69:69),'(a1)') londir
            write (i2,'(a)') card 
c           write (*,'(a)') card 

C  Now write output *86* records

         elseif (card(8:9) == '86') then
c           write (*,'(a)') card 
            write (card(46:52),'(i7)') nint(ehtnew*1000)
            write (i2,'(a)') card 
c           write(*, '(a80)') card
         else
            write (i2,'(a)') card 
         endif
         go to 310
  390    close (i1, status = 'keep')
c        close (i6, status = 'keep')
         close (i2, status = 'keep')

      elseif (option .eq. '3') then
         vxsave = 0.d0
         vysave = 0.d0
         vzsave = 0.d0
         vnsave = 0.d0
         vesave = 0.d0
         vusave = 0.d0
         write (luout, 400)
  400    format (' Enter name of input file: ')     
         read (luin,'(a)',err=608,iostat=ios) nameif
         if (ios /= 0) goto 608
         open (i1, file = nameif, status = 'old')
C        open (i6,file='transformed_'//nameif)

C  write some comments in the transformed files

         call extract_name (nameif,iii)
c        write (I2,1309) nameif(1:iii)
c        write (I2,1309) trim(nameif)
         write (i2,1310) HTDP_version
         write (i2,1311) frame2       
         WRITE (I2, 1312) MONTH2, IDAY2, IYEAR2, DATE2
c        const = 180.d0/PI

  410    read (i1,'(a)',end=450,err=607,iostat=ios) record
         if (ios /= 0) goto 607
         call interprate_XYZ_record (record,xlat,xlon,eht,name24)
         ylat = (xlat*3600.d0) / rhosec
         ylon = (xlon*3600.d0) / rhosec
         elon = -ylon
         call TOXYZ(ylat, elon, eht, x, y, z)

C  Modified in 09/15/2014 by JS to accommodate the IERS 96-97 trans. par.

         if (Is_inp_NAD83 .or. Is_out_NAD83) then
           call TOIT94 (x,y,z,x1,y1,z1,date1,iopt1)
           call FRIT94 (x1,y1,z1,x2,y2,z2,date1,iopt2)
         else
           call TOIT94_IERS(x, y, z, x1, y1, z1, date1, iopt1)
           call FRIT94_IERS(x1, y1, z1, x2, y2, z2, date1, iopt2)
         endif

C  End changes

         if (min1 .ne. min2) then
            if(.not.FRMXYZ(x2,y2,z2,ylat1,elon1,eht1)) STOP 666
            ylon1 = -elon1
            if (ylon1 .lt. 0.d0) ylon1 = ylon1 + twopi
            call PREDV (ylat1, ylon1, eht1, date1, iopt1, jregn, 
     1         vn, ve, vu)
            if (jregn .eq. 0) then
               write(i2, 330) name24
               go to 410
            else
               call TOVXYZ(ylat1,elon1,vn,ve,vu,vx,vy,vz)
               vxsave = vx
               vysave = vy
               vzsave = vz
               vnsave = vn
               vesave = ve
               vusave = vu

C  Modified on 09/15/2014 by JS to accomodate the IERS trans. par.

               if (Is_inp_NAD83 .or. Is_out_NAD83) then
                 call VTRANF(x,y,z,vx,vy,vz,iopt1,iopt2)
               else
                 call VTRANF_IERS( x, y, z, vx, vy, vz, iopt1, iopt2)
               endif
C  End changes

               call TOVNEU(ylat1, elon1,vx,vy,vz,vn,ve,vu)
               call NEWCOR (ylat1, ylon1, eht1, min1, min2,
     1            ylatt,ylont,ehtnew,dn,de,du,vn,ve,vu)
               call TOXYZ(ylatt,-ylont,ehtnew,xt,yt,zt)
            endif
         else
            xt = x2
            yt = y2
            zt = z2
            if(.not.FRMXYZ(xt,yt,zt,ylatt,elont,ehtnew)) STOP 666
            ylont = -elont
            if (ylont .lt. 0.d0) ylont = ylont + twopi
         endif
c        call PRNTTP(x,y,z,xt,yt,zt,ylat,ylatt,ylon,ylont,eht,ehtnew,
c    1      name24,0,vxsave,vysave,vzsave,vnsave,vesave,vusave)
         outlat = ylatt*rad2deg 
         outlon = ylont*rad2deg 
         call extract_name (name24,iii)

C  Added by JS on 07/22/2015 to limit the use in NAD83 to the US
c  Reversed on 05/11/2017

c        if (Is_inp_NAD83 .or. Is_out_NAD83) then
c          lat = ylat*rad2deg ; lon = 360.d0 - ylon*rad2deg    
c          if(.not. am_I_in_or_near_AK  (lat,lon) .and.
c    &       .not. am_I_in_or_near_AS   (lat,lon) .and.
c    &       .not. am_I_in_or_near_CONUS(lat,lon) .and.
c    &       .not. am_I_in_or_near_CQ   (lat,lon) .and.
c    &       .not. am_I_in_or_near_Guam (lat,lon) .and.
c    &       .not. am_I_in_or_near_HI   (lat,lon) .and.
c    &       .not. am_I_in_or_near_PR   (lat,lon) .and.
c    &       .not. am_I_in_or_near_VQ   (lat,lon) .and.
c    &       .not. am_I_in_or_near_KW   (lat,lon)) then
c            write (i2,'(10x,a,8x,a)') 
c    &     "FYI, NAD83(2011) is defined only in US territories",name24
c          endif
c        endif
         write (i2,449) outlat,outlon,ehtnew,name24(1:iii)

C  Until here on 07/22/2015 to limit the use in NAD83 to the US

c        write (i6,449) outlat,outlon,ehtnew,trim(name24)
 449     format (2f16.10,f10.3,4x,a)
         go to 410
      
  450    close(i1, status = 'keep')
c        close(i6, status = 'keep')
         close(i2, status = 'keep')

      elseif (option .eq. '4') then
         vxsave = 0.d0
         vysave = 0.d0
         vzsave = 0.d0
         vnsave = 0.d0
         vesave = 0.d0
         vusave = 0.d0
         write (luout, 4001)
4001     format (' Enter name of input file: ')
         read (luin, '(a)',err=600,iostat=ios) nameif
         if (ios /= 0) goto 600
         open (i1, file = nameif, status = 'old')
c        open (i6,file='transformed_'//nameif)

C  write some comments in the transformed files

         call extract_name (nameif,iii)
c        write (i2,1309) nameif(1:iii)
c        write (i2,1309) trim(nameif)
c        write (i6,1310) HTDP_version
c        write (i6,1311) frame2       
c        WRITE (i6, 1312) MONTH2, IDAY2, IYEAR2, DATE2
c        const = 180.d0/PI
  411    read (i1,'(a)',end=451,err=607,iostat=ios) record         
         if (ios /= 0) goto 607
         call interprate_XYZ_record (record,x,y,z,name24)
         if(.not.FRMXYZ(x,y,z,ylat,elon,eht)) STOP 666
         ylon = -elon
         if (ylon .lt. 0.d0) ylon = ylon + twopi

C  Modified in 09/15/2014 by JS to accommodate the IERS 96-97 trans. par.

         if (Is_inp_NAD83 .or. Is_out_NAD83) then
           call TOIT94 (x,y,z,x1,y1,z1,date1,iopt1)
           call FRIT94 (x1,y1,z1,x2,y2,z2,date1,iopt2)
         else
           call TOIT94_IERS(x, y, z, x1, y1, z1, date1, iopt1)
           call FRIT94_IERS(x1, y1, z1, x2, y2, z2, date1, iopt2)
         endif

C  End changes

         if (min1 .ne. min2) then
           if(.not.FRMXYZ(x2,y2,z2,ylat1,elon1,eht1)) STOP 666
           ylon1 = -elon1
           if (ylon1 .lt. 0.d0) ylon1 = ylon1 + twopi
           call PREDV (ylat1,ylon1,eht1,date1,iopt1,jregn,vn,ve,vu)
           if (jregn .eq. 0) then
             write(i2, 330) name24
             go to 411
           else
             call TOVXYZ(ylat1,elon1,vn,ve,vu,vx,vy,vz)
             vxsave = vx
             vysave = vy
             vzsave = vz
             vnsave = vn
             vesave = ve
             vusave = vu
             write (*,*) vn,ve,vu
             write (*,*) vx,vy,vz

C  Modified on 09/15/2014 by JS to accomodate the IERS trans. par.

             if (Is_inp_NAD83 .or. Is_out_NAD83) then
               call VTRANF(x,y,z,vx,vy,vz,iopt1,iopt2)
             else
               call VTRANF_IERS( x, y, z, vx, vy, vz, iopt1, iopt2)
             endif
C  End changes

             call TOVNEU(ylat1, elon1,vx,vy,vz,vn,ve,vu)
             write (*,*) vn,ve,vu
             call NEWCOR (ylat1, ylon1, eht1, min1, min2,
     1          ylatt,ylont,ehtnew,dn,de,du,vn,ve,vu)
             call TOXYZ(ylatt,-ylont,ehtnew,xt,yt,zt)
           endif
         else
           xt = x2
           yt = y2
           zt = z2
           if(.not.FRMXYZ(xt,yt,zt,ylatt,elont,ehtnew)) STOP 666
           ylont = -elont
           if (ylont .lt. 0.d0) ylont = ylont + twopi
         endif
c        call PRNTTP(x,y,z,xt,yt,zt,ylat,ylatt,ylon,ylont,eht,ehtnew,
c    1      name24,0,vxsave,vysave,vzsave,vnsave,vesave,vusave)
         call extract_name (name24,iii)

C  Added by JS on 07/22/2015 to limit the use in NAD83 to the US

c        if (Is_inp_NAD83 .or. Is_out_NAD83) then
c          lat = ylatt*rad2deg ; lon = 360.d0 - ylont*rad2deg    
c          if(.not. am_I_in_or_near_AK  (lat,lon) .and.
c    &       .not. am_I_in_or_near_AS   (lat,lon) .and.
c    &       .not. am_I_in_or_near_CONUS(lat,lon) .and.
c    &       .not. am_I_in_or_near_CQ   (lat,lon) .and.
c    &       .not. am_I_in_or_near_Guam (lat,lon) .and.
c    &       .not. am_I_in_or_near_HI   (lat,lon) .and.
c    &       .not. am_I_in_or_near_PR   (lat,lon) .and.
c    &       .not. am_I_in_or_near_VQ   (lat,lon) .and. 
c    &       .not. am_I_in_or_near_KW   (lat,lon)) then
c            write (i2,'(10x,a,8x,a)') 
c    &       "FYI, NAD83(2011) is defined only in US territories",name24
c            write (i2,1449) xt,yt,zt,name24(1:iii)
c          endif
c        else
           write (i2,1449) xt,yt,zt,name24(1:iii)
c        endif

C  Until here on 07/22/2015 to limit the use in NAD83 to the US

c        write (i6,1449) xt,yt,zt,trim(name24)
 1449    format (3f20.3,4x,a)
         go to 411
      
  451    close(i1, status = 'keep')
         close(i2, status = 'keep')

      elseif (option .eq. '5') then
         vxsave = 0.d0
         vysave = 0.d0
         vzsave = 0.d0
         vnsave = 0.d0
         vesave = 0.d0
         vusave = 0.d0
         write (luout, 5001)
5001     format (' Enter name of input file: ')
         read (luin, '(a)',err=600,iostat=ios) nameif
         if (ios /= 0) goto 600
         open (i1, file = nameif, status = 'old')
c        open (i6,file='transformed_'//nameif)

C  write some comments in the transformed files

         call extract_name (nameif,iii)
c        write (i2,1309) nameif(1:iii)
c        write (i2,1309) trim(nameif)
c        write (i6,1310) HTDP_version
c        write (i6,1311) frame2       
c        WRITE (i6, 1312) MONTH2, IDAY2, IYEAR2, DATE2
c        const = 180.d0/PI
  511    read (i1,'(a)',end=551,err=607,iostat=ios) record         
         if (ios /= 0) goto 607
         call interprate_XYZVxVyVz_record (record,x,y,z,Vx,Vy,Vz,name24)
c        write (*,*) X,Y,Z
c        write (*,*) Vx,Vy,Vz
         if(.not.FRMXYZ(x,y,z,ylat,elon,eht)) STOP 666
         ylon = -elon
         if (ylon .lt. 0.d0) ylon = ylon + twopi

C  Modified in 09/15/2014 by JS to accommodate the IERS 96-97 trans. par.

         if (Is_inp_NAD83 .or. Is_out_NAD83) then
           call TOIT94 (x,y,z,x1,y1,z1,date1,iopt1)
           call FRIT94 (x1,y1,z1,x2,y2,z2,date1,iopt2)
         else
           call TOIT94_IERS(x, y, z, x1, y1, z1, date1, iopt1)
           call FRIT94_IERS(x1, y1, z1, x2, y2, z2, date1, iopt2)
         endif

C  End changes on 09/15/2014

         if (min1 .ne. min2) then
           if(.not.FRMXYZ(x2,y2,z2,ylat1,elon1,eht1)) STOP 666
           ylon1 = -elon1
           if (ylon1 .lt. 0.d0) ylon1 = ylon1 + twopi
c          call PREDV (ylat1,ylon1,eht1,date1,iopt1,jregn,vn,ve,vu)
c          if (jregn .eq. 0) then
c            write(i2, 330) name24
c            go to 411
c          else
c            call TOVXYZ(ylat1,elon1,vn,ve,vu,vx,vy,vz)

             call TOVNEU(ylat1,elon1,vx,vy,vz,vn,ve,vu)
             vxsave = vx
             vysave = vy
             vzsave = vz
             vnsave = vn
             vesave = ve
             vusave = vu
c            write (*,*) vn,ve,vu

C  Modified on 09/15/2014 by JS to accomodate the IERS trans. par.

             if (Is_inp_NAD83 .or. Is_out_NAD83) then
               call VTRANF(x,y,z,vx,vy,vz,iopt1,iopt2)
             else
               call VTRANF_IERS( x, y, z, vx, vy, vz, iopt1, iopt2)
             endif
C  End changes

             call TOVNEU(ylat1, elon1,vx,vy,vz,vn,ve,vu)
             call NEWCOR (ylat1, ylon1, eht1, min1, min2,
     1          ylatt,ylont,ehtnew,dn,de,du,vn,ve,vu)
             call TOXYZ(ylatt,-ylont,ehtnew,xt,yt,zt)
c          endif
         else
           xt = x2
           yt = y2
           zt = z2
           if(.not.FRMXYZ(xt,yt,zt,ylatt,elont,ehtnew)) STOP 666
           ylont = -elont
           if (ylont .lt. 0.d0) ylont = ylont + twopi
         endif
c        call PRNTTP(x,y,z,xt,yt,zt,ylat,ylatt,ylon,ylont,eht,ehtnew,
c    1      name24,0,vxsave,vysave,vzsave,vnsave,vesave,vusave)
         call extract_name (name24,iii)

C  Added by JS on 07/22/2015 to limit the use in NAD83 to the US

c        if (Is_inp_NAD83 .or. Is_out_NAD83) then
c          lat = ylatt*rad2deg ; lon = 360.d0 - ylont*rad2deg    
c          if(.not. am_I_in_or_near_AK  (lat,lon) .and.
c    &       .not. am_I_in_or_near_AS   (lat,lon) .and.
c    &       .not. am_I_in_or_near_CONUS(lat,lon) .and.
c    &       .not. am_I_in_or_near_CQ   (lat,lon) .and.
c    &       .not. am_I_in_or_near_Guam (lat,lon) .and.
c    &       .not. am_I_in_or_near_HI   (lat,lon) .and.
c    &       .not. am_I_in_or_near_PR   (lat,lon) .and.
c    &       .not. am_I_in_or_near_VQ   (lat,lon) .and. 
c    &       .not. am_I_in_or_near_KW   (lat,lon)) then
c            write (i2,'(10x,a,8x,a)') 
c    &       "FYI, NAD83(2011) is defined only in US territories",name24
c            write (i2,1449) xt,yt,zt,name24(1:iii)
c          endif
c        else
           write (i2,1449) xt,yt,zt,name24(1:iii)
c        endif

C  Until here on 07/22/2015 to limit the use in NAD83 to the US

c        write (i6,1449) xt,yt,zt,trim(name24)
         go to 511
      
  551    close(i1, status = 'keep')
         close(i2, status = 'keep')
      else
         write(luout,*) ' Improper selection -- try again.  '
         go to 160
      endif

  500 continue
c     close(i2, status = 'keep')
      return

  600 write (*,'(/)') 
      write (*,*) "Wrong output file name in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)') 
      write (*,*) "Wrong answer in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)') 
      write (*,*) "Wrong option in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)') 
      write (*,*) "Wrong bbname in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  604 write (*,'(/)') 
      write (*,*) "Wrong card in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  605 write (*,'(/)') 
      write (*,*) "Failed to read 80 record in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  606 write (*,'(/)') 
      write (*,*) "Failed to read 86 record in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  607 write (*,*) "Failed reading input file in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  608 write (*,'(/)') 
      write (*,*) "Wrong input file name in TRFPOS: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      end
C*******************************************8*****************************************************
      subroutine PRNTTP(x, y, z, x1, y1, z1, ylat,            
     1  ylatt, ylon, ylont, eht, ehtnew, name24, iprint,
     2  vx, vy, vz, vn, ve, vu )

** Print updated and/or transformed parameters               

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER    NAME24*24
      CHARACTER    LATDIR*1, LONDIR*1, LATDR*1, LONDR*1

            CALL TODMSS(YLATT,LATDN,LATMN,SLATN,ISIGN)
            LATDIR = 'N'     
            if (isign .eq. -1) LATDIR = 'S'
            CALL TODMSS(YLONT,LONDN,LONMN,SLONN,ISIGN)
            LONDIR = 'W'
            if (isign .eq. -1) LONDIR = 'E'

      if (iprint .eq. 1) then
            WRITE(LUOUT,1065)LATDN,LATMN,SLATN,LATDIR,LONDN,LONMN,
     1                    SLONN,LONDIR, EHTNEW, X1,Y1,Z1
 1065       FORMAT(' ****************************************'/
     1         ' New latitude   = ',I3,1X,I2,1X,F8.5,1X,A1  /
     2         ' New longitude  = ',I3,1x,I2,1X,F8.5,1X,A1  /
     2         ' New Ellip. Ht. = ', F12.3 ,' meters' /
     3         ' New X          = ',F12.3  ,' meters' /
     4         ' New Y          = ',F12.3  ,' meters' /
     5         ' New Z          = ',F12.3  ,' meters' /
     3         ' ****************************************'/)
      endif

            call TODMSS(ylat,latd,latm,slat,isign)
            latdr = 'N'
            if(isign .eq. -1) latdr = 'S'
            call TODMSS(ylon, lond,lonm,slon,isign)
            londr = 'W'
            if(isign .eq. -1) londr = 'E'
            WRITE(I2,1070)NAME24,
     1          LATD,LATM,SLAT,latdr,LATDN,LATMN,SLATN,latdir,vn,
     1          LOND,LONM,SLON,londr,LONDN,LONMN,SLONN,londir,ve,
     1          EHT, EHTNEW, vu, X, X1, vx,
     1          Y, Y1, vy,  Z, Z1, vz
 1070       FORMAT(1X,A24,/
     1   2X, 'LATITUDE  ',2(2X,I3,1X,I2.2,1X,F8.5,1X,A1,2X),
     1         2X, f8.2, ' mm/yr  north' /
     1   2X, 'LONGITUDE ',2(2X,I3,1X,I2.2,1X,F8.5,1X,A1,2X),
     1         2X, f8.2, ' mm/yr  east' /
     1   2X, 'ELLIP. HT.',2(6X,F14.3),' m  ',
     1         f8.2, ' mm/yr  up' /
     1   2X, 'X         ',2(6X,F14.3),' m  ', f8.2, ' mm/yr' /
     1   2X, 'Y         ',2(6X,F14.3),' m  ', f8.2, ' mm/yr' /
     1   2X, 'Z         ',2(6X,F14.3),' m  ', f8.2, ' mm/yr'/ )
         RETURN
         END

**********************************************************
      subroutine SETTP

*** Specify transformation parameters from ITRF94
*** to other reference frames
**************************************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97
C  The latter parameters were added to HTDP in 09/2014. They will be used to transform between
C  ITRF systems. They will not be used if the transformation involves NAD83 or WGS84 (transit).
C  They will be used for the Pacific branches of NAD83.


      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /const/ a, f, e2, eps, af, pi, twopi, rhosec
      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)
      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

C  Parameters computed with the IGS values of ITRF96==>ITRF97

*** From ITRF94 to NAD 83
      tx(1) = 0.9910d0        
      ty(1) = -1.9072d0      
      tz(1) = -.5129d0      
      dtx(1) = 0.d0        
      dty(1) = 0.d0       
      dtz(1) = 0.d0      
      rx(1) = 1.25033d-7
      ry(1) = 0.46785d-7
      rz(1) = 0.56529d-7
      drx(1) = 0.00258d-7
      dry(1) = -.03599d-7
      drz(1) = -.00153d-7
      scale(1) = 0.d0
      dscale(1) = 0.0d0
      refepc(1) = 1997.0d0

*** From ITRF94 to ITRF88
      tx(2) = 0.018d0
      ty(2) = 0.000d0
      tz(2) = -.092d0
      dtx(2) = 0.0d0    
      dty(2) = 0.0d0
      dtz(2) = 0.0d0
      rx(2) = -.0001d0 / rhosec
      ry(2) = 0.0d0                 
      rz(2) = 0.0d0 
      drx(2) = 0.0d0                  
      dry(2) = 0.0d0              
      drz(2) = 0.0d0                   
      scale(2) = 0.74d-8
      dscale(2) = 0.0d0
      refepc(2) = 1988.0d0

*** From ITRF94 to ITRF89
      tx(3) = 0.023d0
      ty(3) = 0.036d0
      tz(3) = -.068d0
      dtx(3) = 0.0d0    
      dty(3) = 0.0d0
      dtz(3) = 0.0d0
      rx(3) = 0.0d0
      ry(3) = 0.0d0                 
      rz(3) = 0.0d0 
      drx(3) = 0.0d0                  
      dry(3) = 0.0d0              
      drz(3) = 0.0d0                   
      scale(3) = 0.43d-8
      dscale(3) = 0.0d0
      refepc(3) = 1988.0d0

*** From ITRF94 to ITRF90
      tx(4) = 0.018d0
      ty(4) = 0.012d0
      tz(4) = -.030d0
      dtx(4) = 0.0d0    
      dty(4) = 0.0d0
      dtz(4) = 0.0d0
      rx(4) = 0.0d0
      ry(4) = 0.0d0                 
      rz(4) = 0.0d0 
      drx(4) = 0.0d0                  
      dry(4) = 0.0d0              
      drz(4) = 0.0d0                   
      scale(4) = 0.09d-8
      dscale(4) = 0.0d0
      refepc(4) = 1988.0d0

*** From ITRF94 to ITRF91
      tx(5) = 0.020d0
      ty(5) = 0.016d0
      tz(5) = -.014d0
      dtx(5) = 0.0d0    
      dty(5) = 0.0d0
      dtz(5) = 0.0d0
      rx(5) = 0.0d0
      ry(5) = 0.0d0                 
      rz(5) = 0.0d0 
      drx(5) = 0.0d0                  
      dry(5) = 0.0d0              
      drz(5) = 0.0d0                   
      scale(5) = 0.06d-8
      dscale(5) = 0.0d0
      refepc(5) = 1988.0d0

*** From ITRF94 to ITRF92
      tx(6) = 0.008d0
      ty(6) = 0.002d0
      tz(6) = -.008d0
      dtx(6) = 0.0d0    
      dty(6) = 0.0d0
      dtz(6) = 0.0d0
      rx(6) = 0.0d0
      ry(6) = 0.0d0                 
      rz(6) = 0.0d0 
      drx(6) = 0.0d0                  
      dry(6) = 0.0d0              
      drz(6) = 0.0d0                   
      scale(6) = -.08d-8
      dscale(6) = 0.0d0
      refepc(6) = 1988.0d0

*** From ITRF94 to ITRF93
      tx(7) = 0.006d0
      ty(7) = -.005d0
      tz(7) = -.015d0
      dtx(7) = -.0029d0
      dty(7) = 0.0004d0
      dtz(7) = 0.0008d0
      rx(7) = 0.00039d0 / rhosec
      ry(7) = -.00080d0 / rhosec
      rz(7) = 0.00096d0 / rhosec
      drx(7) = .00011d0 / rhosec
      dry(7) = .00019d0 / rhosec
      drz(7) =-.00005d0 / rhosec
      scale(7) = 0.04d-8
      dscale(7) = 0.0d0
      refepc(7) = 1988.0d0

*** From ITRF94 to ITRF96
      tx(8) = 0.d0
      ty(8) = 0.d0
      tz(8) = 0.d0
      dtx(8) = 0.d0
      dty(8) = 0.d0
      dtz(8) = 0.d0
      rx(8) = 0.d0
      ry(8) = 0.d0
      rz(8) = 0.d0
      drx(8) = 0.d0
      dry(8) = 0.d0
      drz(8) = 0.d0
      scale(8) = 0.d0
      dscale(8) = 0.0d0
      refepc(8) = 1996.0d0

*** From ITRF94 to ITRF97 (based on IGS adopted values)
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx(9) = 0.00207d0
      ty(9) = 0.00021d0
      tz(9) = -0.00995d0
      dtx(9) = -0.00069d0
      dty(9) = 0.00010d0
      dtz(9) = -0.00186d0
      rx(9) = -0.00012467d0 / rhosec
      ry(9) = 0.00022355d0 / rhosec
      rz(9) = 0.00006065d0 / rhosec
      drx(9) = -0.00001347d0 / rhosec
      dry(9) = 0.00001514d0 / rhosec
      drz(9) = -0.00000027d0 / rhosec
      scale(9) = 0.93496d-9
      dscale(9) = 0.19201d-9
      refepc(9) = 1997.0d0

*** From ITRF94 to WGS 72 (composition of ITRF94 -> NAD_83 -> WGS_72)
      tx(10) = 0.9910d0
      ty(10) = -1.9072d0
      tz(10) = -.5129d0 - 4.5d0
      dtx(10) = 0.d0
      dty(10) = 0.d0
      dtz(10) = 0.d0
      rx(10) = 1.25033d-7
      ry(10) = 0.46785d-7
      rz(10) = 0.56529d-7 + 26.85868d-7
      drx(10) = 0.00258d-7
      dry(10) = -.03599d-7
      drz(10) = -.00153d-7
      scale(10) = 0.d0 - 0.2263d-6
      dscale(10) = 0.0d0
      refepc(10) = 1997.0d0

*** From ITRF94 to ITRF00
*** assumes that ITRF94 = ITRF96 and
*** uses IGS values for ITRF96 -> ITRF97
*** and IERS values for ITRF97 -> ITRF00
      tx(11) = -.00463d0
      ty(11) = -.00589d0
      tz(11) = +.00855d0
      dtx(11) = -0.00069d0
      dty(11) = 0.00070d0
      dtz(11) = -0.00046d0
      rx(11) = -.00012467d0 / rhosec
      ry(11) = 0.00022355d0 / rhosec
      rz(11) = 0.00006065d0 / rhosec
      drx(11) = -0.00001347d0 / rhosec
      dry(11) = 0.00001514d0 / rhosec
      drz(11) = 0.00001973d0 / rhosec
      scale(11) = -0.61504d-9
      dscale(11) = 0.18201d-9
      refepc(11) = 1997.0d0

*** From ITRF94 to PACP00
*** use PA/ITRF00 rotation rates from Beavan et al., (2002)
      tx(12) = 0.9056d0
      ty(12) = -2.0200d0
      tz(12) = -0.5516d0
      dtx(12) = -.00069d0
      dty(12) = 0.00070d0
      dtz(12) = -0.00046d0
      rx(12) = 0.027616d0 / rhosec
      ry(12) = 0.013692d0 / rhosec
      rz(12) = 0.002773d0 / rhosec
      drx(12) = -.000397d0 / rhosec
      dry(12) = 0.001022d0 / rhosec
      drz(12) = -.002166d0 / rhosec
      scale(12) = -0.61504d-9
      dscale(12) = 0.18201d-9
      refepc(12) = 1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
      tx(13) = 0.9056d0
      ty(13) = -2.0200d0
      tz(13) = -0.5516d0
      dtx(13) = -0.00069d0
      dty(13) = 0.00070d0
      dtz(13) = -0.00046d0
      rx(13) = .028847d0 / rhosec
      ry(13) = .010644d0 / rhosec
      rz(13) = 0.008989d0 / rhosec
      drx(13) = -0.000033d0 / rhosec
      dry(13) =  0.000120d0 / rhosec
      drz(13) = -0.000327d0 / rhosec
      scale(13) = -0.61504d-9
      dscale(13) = 0.18201d-9
      refepc(13) = 1997.00d0

*** From ITRF94 to ITRF2005
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
      tx(14) = -0.00533d0
      ty(14) = -0.00479d0
      tz(14) =  0.00895d0
      dtx(14) = -0.00049d0
      dty(14) =  0.00060d0
      dtz(14) =  0.00134d0
      rx(14) = -.00012467d0 / rhosec
      ry(14) = 0.00022355d0 / rhosec
      rz(14) = 0.00006065d0 / rhosec
      drx(14) = -0.00001347d0 / rhosec
      dry(14) = 0.00001514d0 / rhosec
      drz(14) = 0.00001973d0 / rhosec
      scale(14) = -0.77504d-9
      dscale(14) = 0.10201d-9
      refepc(14) = 1997.0d0

*** From ITRF94 to ITRF2008 (also IGS08 and IGB08)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)

      tx(15) = -0.00243d0
      ty(15) = -0.00389d0
      tz(15) =  0.01365d0
      dtx(15) = -0.00079d0
      dty(15) =  0.00060d0
      dtz(15) =  0.00134d0
      rx(15) = -0.00012467d0 / rhosec
      ry(15) =  0.00022355d0 / rhosec
      rz(15) =  0.00006065d0 / rhosec
      drx(15) = -0.00001347d0 / rhosec
      dry(15) =  0.00001514d0 / rhosec
      drz(15) =  0.00001973d0 / rhosec
      scale(15) = -1.71504d-9
      dscale(15) = 0.10201d-9
      refepc(15) = 1997.0d0

*** From ITRF94 to ITRF2014
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)
*** uses IERS values for ITRF2008 -> ITRF2014

c     tx(16)     = -0.00403d0                 !Differs from ITRF2008
      tx(16)     = -0.01430d0
c     ty(16)     = -0.00579d0                 !Differs from ITRF2008
      ty(16)     =  0.00201d0
c     tz(16)     =  0.01125d0                 !Differs from ITRF2008
      tz(16)     =  0.02867d0

      dtx(16)    = -0.00079d0                 !Like ITRF2008
      dty(16)    =  0.00060d0                 !Like ITRF2008
      dtz(16)    =  0.00144d0                 !Differs from ITRF2008

c     rx(16)     = -0.00012467d0 / rhosec     !Like ITRF2008
      rx(16)     = -0.00029978d0 / rhosec
c     ry(16)     =  0.00022355d0 / rhosec     !Like ITRF2008
      ry(16)     =  0.00042037d0 / rhosec
c     rz(16)     =  0.00006065d0 / rhosec     !Like ITRF2008
      rz(16)     =  0.00031714d0 / rhosec
      drx(16)    = -0.00001347d0 / rhosec     !Like ITRF2008
      dry(16)    =  0.00001514d0 / rhosec     !Like ITRF2008
      drz(16)    =  0.00001973d0 / rhosec     !Like ITRF2008

C     scale(16)  = -1.69504d-9                !Differs from ITRF2008, but only insignificantly
      scale(16)  = -0.36891d-9
      dscale(16) =  0.07201d-9                !Differs from ITRF2008, but only insignificantly
      refepc(16) =  2010.0d0                  !Differs from ITRF2008

C*************************************************************************************************************************
C  Parameters computed with the IERS values of ITRF96==>ITRF97

*** From ITRF94 to NAD 83
      tx1(1) = 0.9910d0        
      ty1(1) = -1.9072d0      
      tz1(1) = -.5129d0      
      dtx1(1) = 0.d0        
      dty1(1) = 0.d0       
      dtz1(1) = 0.d0      
      rx1(1) = 1.25033d-7
      ry1(1) = 0.46785d-7
      rz1(1) = 0.56529d-7
      drx1(1) = 0.00258d-7
      dry1(1) = -.03599d-7
      drz1(1) = -.00153d-7
      scale(1) = 0.d0
      dscale1(1) = 0.0d0
      refepc1(1) = 1997.0d0

*** From ITRF94 to ITRF88
      tx1(2) = 0.018d0
      ty1(2) = 0.000d0
      tz1(2) = -.092d0
      dtx1(2) = 0.0d0    
      dty1(2) = 0.0d0
      dtz1(2) = 0.0d0
      rx1(2) = -.0001d0 / rhosec
      ry1(2) = 0.0d0                 
      rz1(2) = 0.0d0 
      drx1(2) = 0.0d0                  
      dry1(2) = 0.0d0              
      drz1(2) = 0.0d0                   
      scale1(2) = 0.74d-8
      dscale1(2) = 0.0d0
      refepc1(2) = 1988.0d0

*** From ITRF94 to ITRF89
      tx1(3) = 0.023d0
      ty1(3) = 0.036d0
      tz1(3) = -.068d0
      dtx1(3) = 0.0d0    
      dty1(3) = 0.0d0
      dtz1(3) = 0.0d0
      rx1(3) = 0.0d0
      ry1(3) = 0.0d0                 
      rz1(3) = 0.0d0 
      drx1(3) = 0.0d0                  
      dry1(3) = 0.0d0              
      drz1(3) = 0.0d0                   
      scale1(3) = 0.43d-8
      dscale1(3) = 0.0d0
      refepc1(3) = 1988.0d0

*** From ITRF94 to ITRF90
      tx1(4) = 0.018d0
      ty1(4) = 0.012d0
      tz1(4) = -.030d0
      dtx1(4) = 0.0d0    
      dty1(4) = 0.0d0
      dtz1(4) = 0.0d0
      rx1(4) = 0.0d0
      ry1(4) = 0.0d0                 
      rz1(4) = 0.0d0 
      drx1(4) = 0.0d0                  
      dry1(4) = 0.0d0              
      drz1(4) = 0.0d0                   
      scale1(4) = 0.09d-8
      dscale1(4) = 0.0d0
      refepc1(4) = 1988.0d0

*** From ITRF94 to ITRF91
      tx1(5) = 0.020d0
      ty1(5) = 0.016d0
      tz1(5) = -.014d0
      dtx1(5) = 0.0d0    
      dty1(5) = 0.0d0
      dtz1(5) = 0.0d0
      rx1(5) = 0.0d0
      ry1(5) = 0.0d0                 
      rz1(5) = 0.0d0 
      drx1(5) = 0.0d0                  
      dry1(5) = 0.0d0              
      drz1(5) = 0.0d0                   
      scale1(5) = 0.06d-8
      dscale1(5) = 0.0d0
      refepc1(5) = 1988.0d0

*** From ITRF94 to ITRF92
      tx1(6) = 0.008d0
      ty1(6) = 0.002d0
      tz1(6) = -.008d0
      dtx1(6) = 0.0d0    
      dty1(6) = 0.0d0
      dtz1(6) = 0.0d0
      rx1(6) = 0.0d0
      ry1(6) = 0.0d0                 
      rz1(6) = 0.0d0 
      drx1(6) = 0.0d0                  
      dry1(6) = 0.0d0              
      drz1(6) = 0.0d0                   
      scale1(6) = -.08d-8
      dscale1(6) = 0.0d0
      refepc1(6) = 1988.0d0

*** From ITRF94 to ITRF93
      tx1(7) = 0.006d0
      ty1(7) = -.005d0
      tz1(7) = -.015d0
      dtx1(7) = -.0029d0
      dty1(7) = 0.0004d0
      dtz1(7) = 0.0008d0
      rx1(7) = 0.00039d0 / rhosec
      ry1(7) = -.00080d0 / rhosec
      rz1(7) = 0.00096d0 / rhosec
      drx1(7) = .00011d0 / rhosec
      dry1(7) = .00019d0 / rhosec
      drz1(7) =-.00005d0 / rhosec
      scale1(7) = 0.04d-8
      dscale1(7) = 0.0d0
      refepc1(7) = 1988.0d0

*** From ITRF94 to ITRF96
      tx1(8) = 0.d0
      ty1(8) = 0.d0
      tz1(8) = 0.d0
      dtx1(8) = 0.d0
      dty1(8) = 0.d0
      dtz1(8) = 0.d0
      rx1(8) = 0.d0
      ry1(8) = 0.d0
      rz1(8) = 0.d0
      drx1(8) = 0.d0
      dry1(8) = 0.d0
      drz1(8) = 0.d0
      scale1(8) = 0.d0
      dscale1(8) = 0.0d0
      refepc1(8) = 1996.0d0

*** From ITRF94 to ITRF97 (based on IERS adopted values)
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx1(9)     = 0.00000d0
      ty1(9)     = 0.00000d0
      tz1(9)     = 0.00000d0
      dtx1(9)    = 0.00000d0
      dty1(9)    = 0.00000d0
      dtz1(9)    = 0.00000d0
      rx1(9)     = 0.00000000d0 
      ry1(9)     = 0.00000000d0
      rz1(9)     = 0.00000000d0 
      drx1(9)    = 0.00000000d0 
      dry1(9)    = 0.00000000d0 
      drz1(9)    = 0.00000000d0 
      scale1(9)  = 0.d0          
      dscale1(9) = 0.d0          
      refepc1(9) = 2000.0d0

*** From ITRF94 to WGS 72 (composition of ITRF94 -> NAD_83 -> WGS_72)
      tx1(10) = 0.9910d0
      ty1(10) = -1.9072d0
      tz1(10) = -.5129d0 - 4.5d0
      dtx1(10) = 0.d0
      dty1(10) = 0.d0
      dtz1(10) = 0.d0
      rx1(10) = 1.25033d-7
      ry1(10) = 0.46785d-7
      rz1(10) = 0.56529d-7 + 26.85868d-7
      drx1(10) = 0.00258d-7
      dry1(10) = -.03599d-7
      drz1(10) = -.00153d-7
      scale1(10) = 0.d0 - 0.2263d-6
      dscale1(10) = 0.0d0
      refepc1(10) = 1997.0d0

*** From ITRF94 to ITRF00
*** assumes that         ITRF94 = ITRF96 and
*** uses IERS values for ITRF96 -> ITRF97
*** and  IERS values for ITRF97 -> ITRF00
      tx1(11)     = -.00670d0
      ty1(11)     = -.00430d0
      tz1(11)     = +.02270d0
      dtx1(11)    = 0.00000d0
      dty1(11)    = 0.00060d0
      dtz1(11)    = 0.00140d0
      rx1(11)     = 0.00000000d0 / rhosec
      ry1(11)     = 0.00000000d0 / rhosec
      rz1(11)     = +0.00006000d0 / rhosec
      drx1(11)    = 0.00000000d0 / rhosec
      dry1(11)    = 0.00000000d0 / rhosec
      drz1(11)    = +0.00002000d0 / rhosec
      scale1(11)  = -1.58000d-9
      dscale1(11) = -0.01000d-9
      refepc1(11) = 2000.0d0

*** From ITRF94 to PACP00
*** use PA/ITRF00 rotation rates from Beavan et al., (2002)
      tx1(12)     = 0.9035d0
      ty1(12)     = -2.0202d0
      tz1(12)     = -0.5417d0
      dtx1(12)    = 0.00000d0
      dty1(12)    = 0.00060d0
      dtz1(12)    = 0.0014d0
      rx1(12)     = 0.027741d0 / rhosec
      ry1(12)     = 0.013469d0 / rhosec
      rz1(12)     = 0.002712d0 / rhosec
      drx1(12)    = -.000384d0 / rhosec
      dry1(12)    = 0.001007d0 / rhosec
      drz1(12)    = -.002166d0 / rhosec
      scale1(12)  = -1.55000d-9
      dscale1(12) = -0.010000d-9
      refepc1(12) = 1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
      tx1(13) = 0.9035d0
      ty1(13) = -2.0202d0
      tz1(13) = -0.5417d0
      dtx1(13) = -0.00000d0
      dty1(13) = 0.00060d0
      dtz1(13) = 0.00140d0
      rx1(13) = .028971d0 / rhosec
      ry1(13) = .01042d0 / rhosec
      rz1(13) = 0.008928d0 / rhosec
      drx1(13) = -0.00002d0 / rhosec
      dry1(13) =  0.000105d0 / rhosec
      drz1(13) = -0.000327d0 / rhosec
      scale1(13) = -1.55000d-9
      dscale1(13) = -0.01000d-9
      refepc1(13) = 1997.00d0

*** From ITRF94 to ITRF2005
*** assumes that         ITRF94   = ITRF96
*** uses IERS values for ITRF96   -> ITRF97
*** uses IERS values for ITRF97   -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
      tx1(14)     = -0.00680d0
      ty1(14)     = -0.00350d0
      tz1(14)     =  0.0285d0
      dtx1(14)    =  0.00020d0
      dty1(14)    =  0.00050d0
      dtz1(14)    =  0.00320d0
      rx1(14)     = 0.00000000d0 / rhosec
      ry1(14)     = 0.00000000d0 / rhosec
      rz1(14)     = +0.00006000d0 / rhosec
      drx1(14)    = 0.00000000d0 / rhosec
      dry1(14)    = 0.00000000d0 / rhosec
      drz1(14)    = +0.00002000d0 / rhosec
      scale1(14)  = -1.98000d-9
      dscale1(14) = -0.09000d-9
      refepc1(14) = 2000.0d0

*** From ITRF94 to ITRF2008 (also IGS08 and IGB08)
*** assumes that ITRF94 = ITRF96
*** uses IERS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)
      tx1(15)     = -0.00480d0
      ty1(15)     = -0.00260d0
      tz1(15)     =  0.03320d0
      dtx1(15)    = -0.00010d0
      dty1(15)    =  0.00050d0
      dtz1(15)    =  0.00320d0
      rx1(15)     =  0.00000000d0 / rhosec
      ry1(15)     =  0.00000000d0 / rhosec
      rz1(15)     = +0.00006000d0 / rhosec
      drx1(15)    = -0.00000000d0 / rhosec
      dry1(15)    =  0.00000000d0 / rhosec
      drz1(15)    = +0.00002000d0 / rhosec
      scale1(15)  = -2.92d-9
      dscale1(15) = -0.09d-9
      refepc1(15) = 2000.0d0

*** From ITRF94 to ITRF2014 
*** assumes that ITRF94 = ITRF96
*** uses IERS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008 (and IGS08 and IGB08)
*** uses IERS values for ITRF2008 -> ITRF2014

c     tx1(16)     = -0.00640d0                       !Differ from ITRF2008
      tx1(16)     = -0.00740d0
c     ty1(16)     = -0.00450d0                       !Differ from ITRF2008
      ty1(16)     =  0.00050d0
c     tz1(16)     =  0.03080d0                       !Differ from ITRF2008
      tz1(16)     =  0.06280d0

      dtx1(16)    = -0.00010d0                       !Like ITRF2008
      dty1(16)    =  0.00050d0                       !Like ITRF2008
      dtz1(16)    =  0.00330d0                       !Differ from ITRF2008

      rx1(16)     =  0.00000000d0 / rhosec           !Like ITRF2008
      ry1(16)     =  0.00000000d0 / rhosec           !Like ITRF2008
c     rz1(16)     = +0.00006000d0 / rhosec           !Like ITRF2008
      rz1(16)     = +0.00026000d0 / rhosec
      drx1(16)    = -0.00000000d0 / rhosec           !Like ITRF2008
      dry1(16)    =  0.00000000d0 / rhosec           !Like ITRF2008
      drz1(16)    = +0.00002000d0 / rhosec           !Like ITRF2008

c     scale1(16)  = -2.90d-9                         !Differ from ITRF2008, but only insignificantly
      scale1(16)  = -3.80d-9 
      dscale1(16) = -0.12d-9                         !Differ from ITRF2008, but only insignificantly

      refepc1(16) = 2010.0d0                         !Differ from ITRF2008

      return
      end
*************************************************************
      subroutine frit94(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts ITRF94 cartesian coordinates to cartesian
*** coordinates in the specified reference frame for the
*** given date
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (x1, y1, z1) --> input ITRF94 coordiates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc(iopt)
      tranx = tx(iopt) + dtx(iopt)*dtime
      trany = ty(iopt) + dty(iopt)*dtime
      tranz = tz(iopt) + dtz(iopt)*dtime
      rotnx  = rx(iopt) + drx(iopt)*dtime
      rotny  = ry(iopt) + dry(iopt)*dtime
      rotnz  = rz(iopt) + drz(iopt)*dtime
      ds     = 1.d0 + scale(iopt) + dscale(iopt)*dtime

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1
c     write (*,*) "FROM TIIT94 ",dtime
c     write (*,*) "FROM TIIT94 ",tx(iopt),ty(iopt),tz(iopt)
c     write (*,*) "FROM TIIT94 ",dtx(iopt),dty(iopt),dtz(iopt)
c     write (*,*) "FROM TIIT94 ",rx(iopt),ry(iopt),rz(iopt)
c     write (*,*) "FROM TIIT94 ",drx(iopt),dry(iopt),drz(iopt)
c     write (*,*) "FROM TIIT94 ",scale(iopt),dscale(iopt)
c     write (*,*) "FROM TIIT94 ",x2,y2,z2

      return
      end
*************************************************************************************
      subroutine frit94_IERS (x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts ITRF94 cartesian coordinates to cartesian
*** coordinates in the specified reference frame for the
*** given date
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (x1, y1, z1) --> input ITRF94 coordiates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)


      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc1(iopt)
      tranx = tx1(iopt) + dtx1(iopt)*dtime
      trany = ty1(iopt) + dty1(iopt)*dtime
      tranz = tz1(iopt) + dtz1(iopt)*dtime
      rotnx  = rx1(iopt) + drx1(iopt)*dtime
      rotny  = ry1(iopt) + dry1(iopt)*dtime
      rotnz  = rz1(iopt) + drz1(iopt)*dtime
      ds     = 1.d0 + scale1(iopt) + dscale1(iopt)*dtime

      x2 = tranx + ds*x1    + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1    + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1
c     write (*,*) "FROM TIIT94_IERS ",dtime
c     write (*,*) "FROM TIIT94_IERS ",tx1(iopt),ty1(iopt),tz1(iopt)
c     write (*,*) "FROM TIIT94_IERS ",dtx1(iopt),dty1(iopt),dtz1(iopt)
c     write (*,*) "FROM TIIT94_IERS ",rx1(iopt),ry1(iopt),rz1(iopt)
c     write (*,*) "FROM TIIT94_IERS ",drx1(iopt),dry1(iopt),drz1(iopt)
c     write (*,*) "FROM TIIT94_IERS ",scale1(iopt),dscale1(iopt)
c     write (*,*) "FROM TIIT94_IERS ",x2,y2,z2

      return
      end

***********************************************************
      subroutine toit94(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts cartesian coordinates in a specified reference
*** to ITRF94 cartesian coordinates for the given date
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (x1, y1, z1) --> input coordiates (meters)
*** (x2, y2, z2) --> output  ITRF94 coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of input reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc(iopt)
      tranx = -(tx(iopt) + dtx(iopt)*dtime)
      trany = -(ty(iopt) + dty(iopt)*dtime)
      tranz = -(tz(iopt) + dtz(iopt)*dtime)
      rotnx  = -(rx(iopt) + drx(iopt)*dtime)
      rotny  = -(ry(iopt) + dry(iopt)*dtime)
      rotnz  = -(rz(iopt) + drz(iopt)*dtime)
      ds     = 1.d0 - (scale(iopt) + dscale(iopt)*dtime)

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1
c     write (*,*) "TO TIIT94 ",dtime
c     write (*,*) "TO TIIT94 ",tx(iopt),ty(iopt),tz(iopt)
c     write (*,*) "TO TIIT94 ",dtx(iopt),dty(iopt),dtz(iopt)
c     write (*,*) "TO TIIT94 ",rx(iopt),ry(iopt),rz(iopt)
c     write (*,*) "TO TIIT94 ",drx(iopt),dry(iopt),drz(iopt)
c     write (*,*) "TO TIIT94 ",scale(iopt),dscale(iopt)
c     write (*,*) "TO TIIT94 ",x2,y2,z2

      return
      end

*****************************************************************
      subroutine toit94_IERS(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts cartesian coordinates in a specified reference
*** to ITRF94 cartesian coordinates for the given date
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (x1, y1, z1) --> input coordiates (meters)
*** (x2, y2, z2) --> output  ITRF94 coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of input reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)


      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc1(iopt)
      tranx = -(tx1(iopt) + dtx1(iopt)*dtime)
      trany = -(ty1(iopt) + dty1(iopt)*dtime)
      tranz = -(tz1(iopt) + dtz1(iopt)*dtime)
      rotnx  = -(rx1(iopt) + drx1(iopt)*dtime)
      rotny  = -(ry1(iopt) + dry1(iopt)*dtime)
      rotnz  = -(rz1(iopt) + drz1(iopt)*dtime)
      ds     = 1.d0 - (scale1(iopt) + dscale1(iopt)*dtime)

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1
c     write (*,*) "TO TIIT94_IERS ",dtime
c     write (*,*) "TO TIIT94_IERS ",tx1(iopt),ty1(iopt),tz1(iopt)
c     write (*,*) "TO TIIT94_IERS ",dtx1(iopt),dty1(iopt),dtz1(iopt)
c     write (*,*) "TO TIIT94_IERS ",rx1(iopt),ry1(iopt),rz1(iopt)
c     write (*,*) "TO TIIT94_IERS ",drx1(iopt),dry1(iopt),drz1(iopt)
c     write (*,*) "TO TIIT94_IERS ",scale1(iopt),dscale1(iopt)
c     write (*,*) "TO TIIT94_IERS ",x2,y2,z2

      return
      end
*****************************************************************
      subroutine MENU1(kopt, mframe)

** Write out options for reference frames

      implicit integer*4 (i-n)
      character    mframe*24, nframe*24
      dimension    nframe(24)
      dimension    iframe(24)
      common /files/ luin, luout, i1, i2, i3, i4, i5, i6

      iframe(1) = 1
      nframe(1) = 'NAD_83(2011/CORS96/2007)'       
      iframe(2) = 12
      nframe(2) = 'NAD_83(PA11/PACP00)     '
      iframe(3) = 13
      nframe(3) = 'NAD_83(MA11/MARP00)     '
      iframe(4) = 10
      nframe(4) = 'WGS_72                  '
      iframe(5) = 1
      nframe(5) = 'WGS_84(transit)         '
c     iframe(6) = 6 (This was incorrect in all versions of HTDP)
      iframe(6) = 5
      nframe(6) = 'WGS_84(G730)            '
      iframe(7) = 8
      nframe(7) = 'WGS_84(G873)            '
      iframe(8) = 11
      nframe(8) = 'WGS_84(G1150)           '
      iframe(9) = 15
      nframe(9) = 'WGS_84(G1674)           '
c     iframe(10)= 4
c     nframe(10)= 'PNEOS_90 or NEOS_90     '
      iframe(10)= 15
      nframe(10)= 'WGS_84(G1762)           '
      iframe(11)= 5
      nframe(11)= 'SIO/MIT_92              '
      iframe(12)= 2
      nframe(12)= 'ITRF88                  '
      iframe(13)= 3
      nframe(13)= 'ITRF89                  '
      iframe(14)= 4
      nframe(14)= 'ITRF90                  '
      iframe(15)= 5
      nframe(15)= 'ITRF91                  '
      iframe(16)= 6
      nframe(16)= 'ITRF92                  '
      iframe(17)= 7
      nframe(17)= 'ITRF93                  '
      iframe(18)= 8
      nframe(18)= 'ITRF94                  '
      iframe(19)= 8
      nframe(19)= 'ITRF96                  '
      iframe(20)= 9
      nframe(20)= 'ITRF97 or IGS97         '
      iframe(21)= 11
      nframe(21)= 'ITRF2000 or IGS00/IGb00 '
      iframe(22)= 14
      nframe(22)= 'ITRF2005 or IGS05       '
      iframe(23)= 15
      nframe(23)= 'ITRF2008 or IGS08/IGB08 '
      iframe(24)= 16
      nframe(24)= 'ITRF2014 or IGS14       '

      write(luout, 100)
  100 format(
     1'  1...NAD_83(2011/CORS96/2007)  (North American plate fixed) '/
     1'  2...NAD_83(PA11/PACP00)       (Pacific plate fixed) '/
     1'  3...NAD_83(MA11/MARP00)       (Mariana plate fixed) '/
     1'                                                   '/
c    1'  4...WGS_72                              '/
     1'  5...WGS_84(transit) (NAD_83(2011) used)  15...ITRF91 '/
     1'  6...WGS_84(G730) (ITRF91 used)           16...ITRF92 '/
     1'  7...WGS_84(G873) (ITRF94 used)           17...ITRF93 '/
     1'  8...WGS_84(G1150) (ITRF2000 used)        18...ITRF94 '/
     1'  9...WGS_84(G1674) (ITRF2008 used)        19...ITRF96 '/
c    1' 10...PNEOS_90 or NEOS_90  (ITRF90 used)   20...ITRF97 or IGS97'/
     1' 10...WGS_84(G1762) (IGb08 used)           20...ITRF97 or IGS97'/
     1' 11...SIO/MIT_92 (ITRF91 used)     21...ITRF2000 or IGS00/IGb00'/
     1' 12...ITRF88                       22...ITRF2005 or IGS05 '/
     1' 13...ITRF89                       23...ITRF2008 or IGS08/IGb08'/
     1' 14...ITRF90 or (PNEOS90/NEOS90)   24...ITRF2014 or IGS14   '/ )
c    1' 14...ITRF90                              '/ )

      read (luin, *,err=50,iostat=ios) iopt
      if (ios /= 0) goto 50
      if ( 1 .le. iopt .and. iopt .le. 24) then
        mframe = nframe(iopt)
        kopt = iframe(iopt)
      else
        mframe = '                '
        kopt = iopt
      endif
      return

 50   write (*,*) 'Failed to read option in MENU1:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      end

***************************************************************
      SUBROUTINE UPDATE(HTDP_version)

** Update positions and/or observations to a specified date.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 16)
      parameter (nbbdim = 10000)
      CHARACTER    OLDBB*30,NEWBB*30, NAMEIF*30
      CHARACTER    NAME24*24
      CHARACTER    OPT*1,ANSWER*1,BBTYPE*1,VOPT*1
      CHARACTER    LATDIR*1, LONDIR*1, LATDR*1, LONDR*1
      character    frame1*24, frame2*24
      character    HTDP_version*10
      LOGICAL      TEST
      LOGICAL      Is_iopt_NAD83
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CAUTION/MONTH2,IDAY2,IYEAR2,DATE2,frame1

      WRITE(LUOUT,20)
   20 FORMAT(' ********************************************'/
     1   ' Please enter the time to which the updated'/     
     1   ' positions and/or observations are to correspond.')
   15 CALL GETMDY(MONTH2,IDAY2,IYEAR2,DATE2,MIN2,TEST)
      IF(TEST) then
         write(luout,*) ' Do you wish to re-enter the time? (y/n)'
         read (luin, '(A1)',err=500,iostat=ios) ANSWER
         if (ios /= 0) goto 500
         if (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 15
         RETURN
      ENDIF

** Choosing reference frame for positions
   35 WRITE(LUOUT,30)
   30 FORMAT(' **************************************************'/
     1   ' Select the reference frame to be used for specifying'/
     2   ' positions. '/)
      call MENU1( iopt, frame1)
      IF(IOPT .LT. 1 .OR. IOPT .GT. numref) THEN
           WRITE(LUOUT,40)
   40      FORMAT(' Improper selection--try again.  ')
           GO TO 35
      ENDIF
      Is_iopt_NAD83 = (iopt == 1)

  999 WRITE(LUOUT,1000)
 1000 FORMAT(' ******************************'/
     1   ' Select option:'/
     2   '     0...No more updates. Return to main menu.'/
     3   '     1...Update positions for individual points', 
     4             ' entered interactively.'/
     5   '     2...Update positions for blue book stations.'/
     6   '     3...Update values for blue book observations.'/
     7   '     4...Update both the positions for blue book',
     8              ' stations'/
     9   '          and the values for blue book',
     9              ' observations.'/
     9   '     5...Update positions contained in batch file '/
     9   '         of delimited records of the form: '/
     9   '         LAT,LON,EHT,"TEXT" ' /
     9   '         LAT = latitude in degrees (positive north/DBL PREC)'/
     9   '         LON = longitude in degrees (positive west/DBL PREC)'/
     9   '         EHT = ellipsoid height in meters (DBL PREC)'/
     9   '         TEXT = Descriptive text (CHARACTER*24) '/
     9   '         Example:  '/
     9   '         40.731671553,112.212671753,34.241,"SALT AIR" '/)
      READ(LUIN,'(A1)',err=501,iostat=ios) OPT
      if (ios /= 0) goto 501
      IF(OPT .eq. '0') THEN
         RETURN
      ELSEIF(OPT .eq. '1') THEN
         WRITE(LUOUT,1020)
 1020    FORMAT(' ************************************'/
     1       ' Enter the time',
     2       ' to which the input positions will correspond.  ')
 1025    CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
         IF(TEST) then
            write(luout,*) ' Do you wish to re-enter the time? (y/n)'
            read(luin,'(A1)',err=500,iostat=ios) ANSWER
            if (ios /= 0) goto 500
            if(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y') GO TO 1025
            RETURN
         ENDIF
         WRITE(LUOUT,1030)
 1030    FORMAT(' ENTER file to save updated positions.  ')
         READ(LUIN,'(A30)',err=502,iostat=ios) NEWBB
         if (ios /= 0) goto 502
         OPEN(I2,FILE=NEWBB,STATUS='UNKNOWN')
         CALL HEADER

         WRITE(I2,1031) frame1     
 1031    FORMAT(' UPDATED POSITIONS IN ', a24)

         WRITE(I2,1040) MONTH1,IDAY1,IYEAR1,MONTH2,IDAY2,IYEAR2,
     1         DATE1, DATE2
 1040    FORMAT(' FROM ',I2,'-',I2.2,'-',I4,
     1     ' TO ',I2,'-',I2.2,'-',I4,' (month-day-year)'/
     1     ' FROM ',F8.3, ' TO ',F8.3, ' (decimal years)'//
     2     16X,'OLD COORDINATE    NEW COORDINATE',
     3     4X,'VELOCITY      DISPLACEMENT',/)
 1050    CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1       LONDIR, NAME24,  X,Y,Z,YLAT,YLON,EHT)
         ELON = -YLON
         call GETVLY(YLAT,ELON,VX,VY,VZ,VN,VE,VU,VOPT,210)
         IF ( VOPT .EQ. '0') then
         call PREDV( ylat, ylon, eht, date1, iopt,
     1         jregn, vn, ve, vu)
         IF(JREGN .eq. 0) THEN
           WRITE(LUOUT,1060)
 1060      FORMAT(' **************************************'/
     1       ' Can not update the position of this point'/
     2       ' because it is outside of the modeled region.'/)
          GO TO 1075
        ELSE
          call TOVXYZ(ylat,elon,vn,ve,vu,vx,vy,vz)
        ENDIF
      ENDIF
         call NEWCOR(ylat,ylon,eht,min1,min2, 
     1               ylatt,ylont,ehtnew,dn,de,du,vn,ve,vu)
         CALL TODMSS(YLATT,LATDN,LATMN,SLATN,ISIGN)
         LATDR = 'N'
         IF (ISIGN .eq. -1) LATDR = 'S'
         CALL TODMSS(YLONT,LONDN,LONMN,SLONN,ISIGN)
         LONDR = 'W'
         IF (ISIGN .eq. -1) LONDR = 'E'
         ELONT = -YLONT
         CALL TOXYZ(YLATT,ELONT,EHTNEW,X1,Y1,Z1)
         DX = X1 - X
         DY = Y1 - Y
         DZ = Z1 - Z
         WRITE(LUOUT,1065)LATDN,LATMN,SLATN,LATDR,LONDN,LONMN,SLONN,
     1              LONDR, EHTNEW, X1,Y1,Z1
 1065    FORMAT(' ****************************************'/
     1         ' Updated latitude   = ',I3,I3,1X,F8.5,1X,A1 /
     2         ' Updated longitude  = ',I3,I3,1X,F8.5,1X,A1 /
     2         ' Updated Ellip. Ht. = ', F12.3 ,' meters' /
     3         ' Updated X          = ',F12.3  ,' meters' /
     4         ' Updated Y          = ',F12.3  ,' meters' /
     5         ' Updated Z          = ',F12.3  ,' meters' /
     3         ' ****************************************'/)
         WRITE(I2,1070)NAME24,
     1      LATD,LATM,SLAT,LATDIR,LATDN,LATMN,SLATN,LATDR,vn,DN,
     1      LOND,LONM,SLON,LONDIR,LONDN,LONMN,SLONN,LONDR,ve,DE,
     1             EHT,  EHTNEW, vu,DU, X, X1, vx,DX,
     1             Y, Y1, vy,DY, Z, Z1, vz,DZ
 1070    FORMAT(1X,A24,/
     1   1X, 'LATITUDE  ',2(2X,I3,I3.2,1X,F8.5,1X,A1),
     1             F8.2, ' mm/yr', f8.3, ' m north'/
     1   1X, 'LONGITUDE ',2(2X,I3,I3.2,1X,F8.5,1x,A1),
     1             f8.2, ' mm/yr',f8.3,' m east'/
     1   1X, 'ELLIP. HT.',2(6X,F13.3),
     1             f8.2, ' mm/yr',F8.3,' m up'/
     1   1X, 'X         ',2(6X,F13.3),
     1             f8.2, ' mm/yr',F8.3,' m '/
     1   1X, 'Y         ',2(6X,F13.3),
     1             f8.2, ' mm/yr',f8.3, ' m'/
     1   1X, 'Z         ',2(6X,F13.3),
     1             f8.2, ' mm/yr',F8.3,' m'/)
 1075    WRITE(LUOUT,1080)
 1080    FORMAT(' Update more positions? (y/n)  ')
         READ(LUIN,'(A1)',err=500,iostat=ios) ANSWER
         if (ios /= 0) goto 500
         IF(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y') GO TO 1050
         CLOSE(I2,STATUS='KEEP')
         RETURN

*** Updating a blue book 

      ELSEIF(OPT .eq. '2' .or. OPT .eq. '3' .or. OPT .eq. '4') THEN

          if (nbbdim .eq. 10000) then
   90        WRITE(LUOUT,91)
   91        FORMAT(
     1       ' Identify type of blue book:'/
     2       '    1...Standard (4-digit SSN)'/
     3       '    2...Non-standard (5-digit SSN)  ')
             READ(LUIN,'(A1)',err=503,iostat=ios) BBTYPE
             if (ios /= 0) goto 503
             IF(BBTYPE .NE. '1' .AND. BBTYPE .NE. '2') GO TO 90
          else
             BBTYPE = '1'
          endif

          WRITE(LUOUT,100)
  100     FORMAT(' Enter name of the blue book file to be updated.'/)
          READ(LUIN,110,err=502,iostat=ios) OLDBB
          if (ios /= 0) goto 502
  110     FORMAT(A30)
          WRITE(LUOUT,120)
  120     FORMAT(' Enter name for the new blue book file that is to '/
     1       'contain the updated information.'/)
          READ(LUIN,110,err=502,iostat=ios) NEWBB
          if (ios /= 0) goto 502
          OPEN(I1,FILE = OLDBB , STATUS = 'OLD')

          OPEN(I4, ACCESS = 'DIRECT', RECL = 40,
     1         FORM = 'UNFORMATTED', STATUS = 'SCRATCH')

          IF(OPT .eq. '2' .or. OPT .eq. '4') THEN
             WRITE(LUOUT,122)
  122        FORMAT(' *********************************************'/
     1           ' Enter the time',
     1           ' to which the input positions correspond.'/)
  123        CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
             IF(TEST) then
              write(luout,*) ' Do you wish to re-enter the time? (y/n)'
              read(luin, '(A1)',err=500,iostat=ios) ANSWER
              if (ios /= 0) goto 500
              IF (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 123
              RETURN
             ENDIF
          ENDIF

*** Retrieve geodetic positions from old blue-book file

         IF(BBTYPE .EQ. '1') THEN
              CALL GETPO4(IOPT, DATE1)
         ELSEIF(BBTYPE .EQ. '2') THEN
              CALL GETPO5(IOPT, DATE1)
         ENDIF

*** Create new blue book file

         OPEN(I2,FILE = NEWBB, STATUS = 'UNKNOWN')

         write (i2,127) HTDP_version
  127    format (' ***CAUTION: This file was processed using HTDP',
     &           ' version ',a10, '***')
         write (i2,128) frame1       
  128    format (' ***CAUTION: Coordinates in this file are in ',    
     &           a24, '***')
         IF (OPT .EQ. '2' .OR. OPT .EQ. '4') THEN
            WRITE(I2, 129) MONTH2, IDAY2, IYEAR2, DATE2
  129       FORMAT(' ***CAUTION: Coordinates in this file have been ',
     *      'updated to ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
         ENDIF

c        IF (OPT .EQ. '3' .OR. OPT .EQ. '4') THEN
c          WRITE(I2, 130) MONTH2, IDAY2, IYEAR2, DATE2
c 130      FORMAT(' ***CAUTION: Observations in this file have been ',
c    *     'updated to ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
c        ENDIF

         IF(BBTYPE .EQ. '1') THEN
           CALL UPBB4(MIN1,MIN2,OPT,IOPT)
         ELSEIF(BBTYPE .EQ. '2') THEN
           CALL UPBB5(MIN1,MIN2,OPT,IOPT)
         ENDIF
         CLOSE(I1, STATUS = 'KEEP')
         CLOSE(I2, STATUS = 'KEEP')

*** Update G-FILE

         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
  605      WRITE(LUOUT,610)
  610      FORMAT(/' Is there a G-FILE to be updated? (y/n)  ')
           READ(LUIN,'(A1)',err=500,iostat=ios) ANSWER
           if (ios /= 0) goto 500
           IF(ANSWER .eq. 'N' .or. ANSWER .eq. 'n') THEN
             CONTINUE
           ELSEIF(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y')THEN
             WRITE(LUOUT,620)
  620        FORMAT(' Enter name of old G-FILE to be updated.  ')
             READ(LUIN,'(A30)',err=502,iostat=ios) OLDBB
             if (ios /= 0) goto 502
             WRITE(LUOUT,630)
  630        FORMAT(' Enter name for the new updated G-FILE.  ')
             READ(LUIN,'(A30)',err=502,iostat=ios) NEWBB
             if (ios /= 0) goto 502
             OPEN(I1,FILE=OLDBB,STATUS='OLD')
             OPEN(I2,FILE=NEWBB,STATUS='UNKNOWN')
             WRITE(I2,130) MONTH2, IDAY2, IYEAR2, DATE2
  130        FORMAT(' ***CAUTION: Observations in this file have been ',
     *       'updated to ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
        
  634        WRITE(LUOUT, 632)
  632        FORMAT(/' ***************************'/
     *              ' To what reference frame should the GPS'/
     *              ' vectors be transformed?'/
     *              '    -1...Do not transform GPS vectors.')
             CALL MENU1(kopt, frame2)
             IF(KOPT .LT. -1 .OR. KOPT .GT. numref) THEN
                WRITE(LUOUT, 40)
                GO TO 634
             ELSEIF (1 .LE. KOPT .AND. KOPT .LE. numref) then
               write( I2, 640) frame2
  640          format( ' ***CAUTION: All GPS interstation vectors',
     1                 ' have been transformed to ', a24, ' ***')
               write( I2, 641) HTDP_version
  641          format(' ***CAUTION: Observations were transformed using'
     1                ,' HTDP version ', a10, ' ***')
             ENDIF
         
             IF(BBTYPE .EQ. '1') THEN
                 CALL UPGFI4(DATE2, MIN2, IOPT, KOPT, 
     *                       MONTH2, IDAY2, IYEAR2)
             ELSEIF(BBTYPE .EQ. '2') THEN
                 CALL UPGFI5(DATE2, MIN2, IOPT, KOPT,
     *                       MONTH2, IDAY2, IYEAR2)
             ENDIF
             CLOSE(I1, STATUS = 'KEEP')
             CLOSE(I2, STATUS = 'KEEP')
           ELSE
             WRITE(LUOUT,700)
             GO TO 605
           ENDIF
         ENDIF
         CLOSE(I4, STATUS = 'DELETE')
         RETURN
      ELSEIF (OPT .eq. '5') then
         write (luout, 710)
  710    format(' Enter name of input file ')
         read( luin, 711,err=502,iostat=ios) NAMEIF
         if (ios /= 0) goto 502
  711    format( a30 )
         open(I1, FILE=NAMEIF, STATUS = 'OLD')
         write(luout, 1020)
  720    call GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
         if(TEST) then
            write(luout,*) ' Do you wish to re-enter the time? (y/n)'
            read(luin, '(A1)',err=500,iostat=ios) ANSWER
            if (ios /= 0) goto 500
            if(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y') GO TO 720
            RETURN
         endif
         write(luout, 1030)
         read(luin, '(a30)',err=502,iostat=ios)NEWBB
         if (ios /= 0) goto 502
         open(I2, FILE=NEWBB, STATUS='UNKNOWN')
         CALL HEADER
         write(I2, 1031) frame1
         write(I2, 1040) MONTH1,IDAY1,IYEAR1,
     1                   MONTH2,IDAY2,IYEAR2,DATE1,DATE2
  730    read(I1,*,END = 750,err=504,iostat=ios) XLAT, XLON, EHT, NAME24
         if (ios /= 0) goto 504
         YLAT = (XLAT*3600.d0) / rhosec
         Ylon = (XLON*3600.d0) / rhosec
         ELON = -YLON
         call TODMSS(ylat,latd,latm,slat,ISIGN)
         if (ISIGN .eq. 1) then
            latdir = 'N'
         else
            latdir = 'S'
         endif
         call TODMSS (ylon,lond,lonm,slon,isign)
         if (isign .eq. 1) then
            londir = 'W'
         else
            londir = 'E'
         endif
         call PREDV(ylat,ylon,eht,date1,iopt,jregn,vn,ve,vu)
         if(jregn .eq. 0) then
            write(I2, 735) name24, latd,latm,slat,latdir,
     1                 lond,lonm,slon,londir
  735       format(1x,a24,/
     1         1x, 'LATITUDE ', 2x,i3,i3.2,1x,f8.5,1x,a1,3x,
     1         'POINT LOCATED OUTSIDE OF MODELED REGION'/
     1         1x, 'LONGITUDE ', 1x,i3,i3.2,1x,f8.5,1x,a1 /)
         else
            call TOVXYZ(ylat,elon,vn,ve,vu,vx,vy,vz)
            call NEWCOR(ylat,ylon,eht,min1,min2,ylatt,ylont,ehtnew,
     1                  dn,de,du,vn,ve,vu)
            call TODMSS(YLATT,latdn,latmn,slatn,isign)
            latdr = 'N'
            if (isign .eq. -1) latdr = 'S'
            call TODMSS(ylont,londn,lonmn,slonn,isign)
            londr = 'W'
            if (isign .eq. -1) londr = 'E'
            call TOXYZ(ylat,elon,eht,x,y,z)
            elont = -ylont
            call TOXYZ(ylatt,elont,ehtnew,x1,y1,z1)
            dx = x1 - x
            dy = y1 - y
            dz = z1 - z
            write(I2, 1070) NAME24, latd,latm,slat,latdir,
     1         latdn,latmn,slatn,latdr,vn,dn,
     1         lond,lonm,slon,londir,londn,lonmn,slonn,londr,ve,de,
     1         eht,ehtnew,vu,du, x, x1,vx,dx,
     1         y, y1,vy,dy, z, z1, vz,dz
          endif
          go to 730
  750     close(I1, STATUS='KEEP')
          close(I2, status = 'KEEP')
          return
      ELSE
         WRITE(LUOUT,700)
  700    FORMAT(' Improper entry !'/)
         GO TO 999
      ENDIF

  500 write (*,'(/)') 
      write (*,*) "Failed to read answer in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  501 write (*,'(/)') 
      write (*,*) "Failed to read OPTION in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  502 write (*,'(/)') 
      write (*,*) "Failed to read file name in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  503 write (*,'(/)') 
      write (*,*) "Failed to read BBTYPE in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  504 write (*,'(/)') 
      write (*,*) "Failed to read input in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*******************************************************************
      SUBROUTINE GETPO4(IOPT, DATE)

*** Retrieve geodetic coordinates from the blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER    JN*1,JW*1
      CHARACTER    TYPE*4
      CHARACTER    CARD*80
      CHARACTER    PIDs*6

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      DO 100 I = 1, 10000
          LOC(I) = 0 
          HT(I) = 0.D0
  100 CONTINUE
          
      JREC = 0
  125 READ(I1,130,END=160,err=200,iostat=ios) CARD
      if (ios /= 0) goto 200
  130 FORMAT(A80)
      TYPE = CARD(7:10)
      IF(TYPE .EQ. '*80*') THEN
         JREC = JREC + 1
         READ(CARD,140,err=201,iostat=ios) ISN,LATD,LATM,SLAT,JN,LOND,
     &                                     LONM,SLON,JW,OH
         if (ios /= 0) goto 201
  140    FORMAT(BZ,10X,I4,T45,2I2,F7.5,A1,I3,I2,F7.5,A1,F6.2)
         LOC(ISN) = JREC
C        HT(ISN) = HT(ISN) + OH
      IF ( HT(ISN) .EQ. 0.D0 ) HT(ISN) = OH
c        write (*,*) SLAT,SLON,HT(ISN)                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
         YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
         IF(JN .EQ. 'S') YLAT = -YLAT
         IF(JW .EQ. 'E') YLON = -YLON
c        write (*,*) ylat,ylon            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call PREDV( ylat, ylon, ht(isn), date, iopt,
     1      IDG, vn, ve, vu)
c        write (*,*) 'IDG ',IDG          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(IDG .eq. 0) THEN
           WRITE(LUOUT,150) CARD
  150      FORMAT('  **WARNING**  The following point is outside of ',
     1     'the modeled region.'/
     2     ' It will be assumed that this point has not moved.'/A80/)
         ENDIF
         WRITE(I4, REC = JREC) YLAT, YLON, VN, VE, VU      
C     ELSEIF(TYPE .eq. '*84*') THEN
C        READ(CARD,155) ISN,GH
C 155    FORMAT(BZ,10X,I4,T70,F6.2)
C        HT(ISN) = HT(ISN) + GH
      ELSEIF(TYPE .eq. '*86*') THEN
         IF(CARD(46:52) .ne. '       ') THEN
            READ(CARD,156,err=202,iostat=ios) ISN, EHT
            if (ios /= 0) goto 202
  156       FORMAT(BZ,10X,I4,T46,F7.3)
            HT(ISN) = EHT
         ELSE
            READ(CARD,157,err=202,iostat=ios) ISN, OHT, GHT
            if (ios /= 0) goto 202
  157       FORMAT(BZ,10X,I4,T17,F7.3,T36,F7.3)
            HT(ISN) = OHT + GHT
         ENDIF
      ENDIF
      GO TO 125
  160 REWIND I1
      RETURN

  200 write (*,'(/)') 
      write (*,*) "Failed 1st reading in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  201 write (*,'(/)') 
      write (*,*) "Failed reading *80* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  202 write (*,'(/)') 
      write (*,*) "Failed reading *86* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
******************************************************************
      SUBROUTINE GETPO5(IOPT, DATE)

*** Retrieve geodetic position from blue book with 5-digit SSN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER   JN*1,JW*1
      CHARACTER   TYPE*3
      CHARACTER   CARD*80
      CHARACTER   PIDs*6

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      DO 100 I = 1, 10000
          LOC(I) = 0
          HT(I) = 0.D0
  100 CONTINUE
          
      JREC = 0

  125 READ(I1,130,END=160,err=200,iostat=ios) CARD
      if (ios /= 0) goto 200
  130 FORMAT(A80)
      TYPE = CARD(7:9)
      IF(TYPE .EQ. '*80') THEN
         JREC = JREC + 1
         READ(CARD,140,err=201,iostat=ios) ISN,LATD,LATM,SLAT,JN,LOND,
     &                                     LONM,SLON,JW,OH
         if (ios /= 0) goto 200
  140    FORMAT(BZ, 9X,I5,T45,2I2,F7.5,A1,I3,I2,F7.5,A1,F6.2)
         LOC(ISN) = JREC
C        HT(ISN) = HT(ISN) + OH
         IF ( HT(ISN) .EQ. 0.D0 ) HT(ISN) = OH
         YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
         YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
         IF(JN .EQ. 'S') YLAT = -YLAT
         IF(JW .EQ. 'E') YLON = -YLON
         call PREDV( ylat, ylon, ht(isn), date, iopt,
     1       idg, vn, ve, vu)
         IF(IDG .eq. 0) THEN
           WRITE(LUOUT,150) CARD
  150      FORMAT('  **WARNING**  The following point is outside of ',
     1     'the modeled region.'/
     2     ' It will be assumed that this point has not moved.'/A80/)
         ENDIF
         WRITE(I4, REC = JREC) YLAT, YLON, VN, VE, VU      
C     ELSEIF(TYPE .eq. '*84') THEN
C        READ(CARD,155) ISN,GH
C 155    FORMAT(BZ, 9X,I5,T70,F6.2)
C        HT(ISN) = HT(ISN) + GH
      ELSEIF(TYPE .eq. '*86') THEN
         IF(CARD(46:52) .ne. '       ') THEN
            READ(CARD,156,err=202,iostat=ios) ISN, EHT
            if (ios /= 0) goto 202
  156       FORMAT(BZ, 9X,I5,T46,F7.3)
            HT(ISN) = EHT
         ELSE
            READ(CARD,157,err=202,iostat=ios) ISN, OHT, GHT
            if (ios /= 0) goto 202
  157       FORMAT(BZ, 9X,I5,T17,F7.3,T36,F7.3)
            HT(ISN) = OHT + GHT
         ENDIF
      ENDIF
      GO TO 125
  160 REWIND I1
      RETURN

  200 write (*,'(/)') 
      write (*,*) "Failed 1st reading in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  201 write (*,'(/)') 
      write (*,*) "Failed reading *80* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  202 write (*,'(/)') 
      write (*,*) "Failed reading *86* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
******************************************************************
      SUBROUTINE UPBB4(MIN1,MIN2,OPT,IOPT)

*** Update blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      parameter (nbbdim = 10000)

      CHARACTER   CARD*80
      CHARACTER   DATE*6
      CHARACTER   TYPE*4
      CHARACTER   OPT*1,JN*1,JW*1
      CHARACTER   LATDIR*1, LONDIR*1
      CHARACTER   PIDs*6                
      LOGICAL TEST, TEST1
      COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim) ,PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C***    To update classical observations an *12* record
C***    is needed to identify the correct century, as
C***    classical observation records contain only a 2-digit year
      IREC12 = 0

  170 READ(I1,175,END=600,err=700,iostat=ios) CARD
      if (ios /= 0) goto 700
  175 FORMAT(A80)
      TYPE = CARD(7:10)
      IF(TYPE .EQ. '*A1*' .OR.
     1      TYPE .EQ. '*AA*' .OR.
     1      TYPE .EQ. '*10*' .OR.
     1      TYPE .EQ. '*11*' .OR.
     1      TYPE .EQ. '*13*' .OR.
     1      TYPE .EQ. '*21*' .OR.
     1      TYPE .EQ. '*25*' .OR.
     1      TYPE .EQ. '*26*' .OR.
     1      TYPE .EQ. '*27*' .OR.
     1      TYPE .EQ. '*28*' .OR.
     1      TYPE .EQ. '*29*' ) THEN
              CONTINUE
      ELSEIF(TYPE .EQ. '*31*' .OR.
     1      TYPE .EQ. '*40*' .OR.
     1      TYPE .EQ. '*41*' .OR.
     1      TYPE .EQ. '*42*' .OR.
     1      TYPE .EQ. '*45*' .OR.
     1      TYPE .EQ. '*46*' .OR.
     1      TYPE .EQ. '*47*' .OR.
     1      TYPE .EQ. '*70*' .OR.
     1      TYPE .EQ. '*81*' .OR.
     1      TYPE .EQ. '*82*' .OR.
     1      TYPE .EQ. '*83*' .OR.
     1      TYPE .EQ. '*84*' .OR.
     1      TYPE .EQ. '*85*' .OR.
     1      TYPE .EQ. '*86*' .OR.
     1      TYPE .EQ. '*90*') THEN
	      CONTINUE

      ELSEIF (TYPE .EQ. '*12*') THEN
         IF( CARD(11:14) .EQ. '    ' .OR.
     1       CARD(17:20) .EQ. '    ') THEN
             WRITE ( LUOUT, 176 )
  176        FORMAT( ' ERROR: The *12* record does not contain'/
     1               '   appropriate dates.')
         ELSE
             READ ( CARD, 177,err=701,iostat=ios) IYEAR1, IYEAR2
             if (ios /= 0) goto 701
  177        FORMAT ( 10X, I4, 2X, I4)
             IF ( (IYEAR2 - IYEAR1) .LT. 0 .OR.
     1            (IYEAR2 - IYEAR1) .GT. 99 ) THEN
                WRITE (LUOUT, 178)
  178           FORMAT(' ERROR: The dates in the *12* record are in'/
     1            ' error. Either they span more than 99 years or '/
     1            ' the end date preceedes the start date.'/
     1            ' If they span more than 99 years, then the Bluebook'/
     1            ' will need to divided into two or more bluebooks'/
     1            ' with the observations in each spanning no more'/
     1            ' than 99 years.')
             ELSE
                IREC12 = 1
             ENDIF
         ENDIF

*** Regular distance

      ELSEIF(TYPE .EQ. '*50*' .OR.
     1       TYPE .EQ. '*51*' .OR.
     1       TYPE .EQ. '*52*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4')THEN
           READ(CARD,180,err=702,iostat=ios) ISN,DATE,JSN,OBS
           if (ios /= 0) goto 702
  180      FORMAT(BZ,10X,I4,T35,A6,T46,I4,T64,F9.4)
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
             CALL TRFDAT(CARD,DATE,IREC12,IYEAR1,IYEAR2,MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)     
             IOBS = IDNINT((OBS+DS)*10000.D0)
             WRITE(CARD(64:72),190) IOBS
  190        FORMAT(I9)
           ENDIF
         ENDIF

*** Azimuth

      ELSEIF(TYPE .EQ. '*60*' .OR. TYPE .EQ. '*61*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
           READ(CARD,200,err=703,iostat=ios)ISN,DATE,JSN,IDEG,MIN,SEC
           if (ios /= 0) goto 703
  200      FORMAT(BZ,10X,I4,T40,A6,T51,I4,T64,I3,I2,F3.1)
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
             CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)    
             OBS = (DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
             OBS = OBS + DA
             IF(OBS .GE. TWOPI) THEN
               OBS = OBS - TWOPI
             ELSEIF(OBS .LT. 0.D0) THEN
               OBS = OBS + TWOPI
             ENDIF
             CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
             ISEC = IDNINT(SEC*10.D0)
             WRITE(CARD(64:71),210) IDEG,MIN,ISEC
  210        FORMAT(I3,I2.2,I3.3)
           ENDIF
         ENDIF

*** Direction observation

      ELSEIF(TYPE .EQ. '*20*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
            READ(CARD,220,err=704,iostat=ios) ISN0,LIST0,DATE,JSN
            if (ios /= 0) goto 704
  220       FORMAT(BZ,10X,I4,I2,T40,A6,T51,I4)
            CALL CHECK(ISN0,JSN,CARD,TEST) 
            CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
            IF(TEST) THEN
              CALL DSDA(LOC(ISN0),LOC(JSN),MINO,MIN2,DS,DA0)     
            ELSE
              DA0 = 0.D0
            ENDIF
         ENDIF
      ELSEIF(TYPE .EQ. '*22*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
           READ(CARD,230,err=705,iostat=ios) ISN,LIST,JSN,IDEG,MIN,SEC
           if (ios /= 0) goto 705
  230      FORMAT(BZ,10X,I4,I2,T51,I4,T64,I3,I2,F4.2)
           IF(LIST.NE.LIST0 .OR. ISN.NE.ISN0) THEN
             WRITE(LUOUT,240) CARD
  240        FORMAT(' Blue-book file has incorrect structure.'/
     1              ' A *22* record disagrees with its corresponding'/
     2              ' *20* record.  The record reads:'/A80)
             STOP
           ENDIF
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)    
             OBS=(DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
             OBS = OBS + DA - DA0
             IF(OBS .GT. TWOPI) THEN
               OBS = OBS - TWOPI
             ELSEIF(OBS .LT. 0.D0) THEN
               OBS = OBS + TWOPI
             ENDIF
             CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
             ISEC = IDNINT(SEC*100.D0)
             WRITE(CARD(64:72),250) IDEG,MIN,ISEC
  250         FORMAT(I3.3,I2.2,I4.4)
           ENDIF
         ENDIF

*** Long distance

      ELSEIF(TYPE .EQ. '*53*' .OR.
     1       TYPE .EQ. '*54*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
           READ(CARD,260,err=706,iostat=ios) ISN,DATE,JSN,OBS
           if (ios /= 0) goto 706
  260      FORMAT(BZ,10X,I4,T35,A6,T46,I4,T64,F10.3)
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
             CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)    
             IOBS = IDNINT((OBS + DS)*1000.D0)
             WRITE(CARD(64:73),270) IOBS
  270        FORMAT(I10)
           ENDIF
         ENDIF
  
*** Horizontal angle

      ELSEIF(TYPE .EQ. '*30*' .OR.
     1       TYPE .EQ. '*32*') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
           READ(CARD,280,err=707,iostat=ios) ISN,JSN,IDEG,MIN,SEC,KSN
           if (ios /= 0) goto 707
  280      FORMAT(BZ,10X,I4,T51,I4,T64,I3,I2,F3.1,I4)
           IF(TYPE.EQ.'*30*') THEN
             DATE = CARD(40:45)
             CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
           ENDIF
           CALL CHECK(ISN,JSN,CARD,TEST)
           CALL CHECK(ISN,KSN,CARD,TEST1)
           IF(TEST .and. TEST1) THEN
              CALL DSDA(LOC(ISN), LOC(JSN), 
     1           MINO, MIN2, DS, DA0)   
              CALL DSDA(LOC(ISN), LOC(KSN), 
     1           MINO, MIN2, DS, DA)     
              OBS=(DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
              OBS = OBS + DA - DA0
              IF(OBS .GE. TWOPI) THEN
                OBS = OBS - TWOPI
              ELSEIF(OBS .LT. 0.D0) THEN
                OBS = OBS + TWOPI
              ENDIF
              CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
              ISEC = IDNINT(SEC*10.D0)
              WRITE(CARD(64:71),290)IDEG,MIN,SEC
  290         FORMAT(I3.3,I2.2,I3.3)
            ENDIF
         ENDIF

*** position record

      ELSEIF(TYPE .EQ. '*80*') THEN
            IF(OPT .eq. '2' .or. OPT .eq. '4') THEN
              READ(CARD,300,err=708,iostat=ios) ISN,LATD,LATM,
     &                              SLAT,JN,LOND,LONM,SLON,JW
              if (ios /= 0) goto 708
  300         FORMAT(BZ,10X,I4,T45,2I2,F7.5,A1,I3,I2,F7.5,A1)
              YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
              YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
              IF(JN .EQ. 'S') YLAT = -YLAT
              IF(JW .EQ. 'E') YLON = -YLON
C             READ(I4, REC = LOC(ISN)) RLAT, RLON, IDG
              READ(I4, REC = LOC(ISN),err=709,iostat=ios) RLAT, 
     &                                        RLON, VN, VE, VU
              if (ios /= 0) goto 709
              call NEWCOR( ylat, ylon, ht(isn), min1, min2,
     1           ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
C             CALL NWCORD(YLAT,YLON,HT(ISN),
C    1          RLAT,RLON,IDG,MIN1,MIN2,
C    1          YLATT,YLONT,EHTNEW,DN,DE,DU,IOPT)
              CALL TODMSS(YLATT,LATDN,LATMN,SLATN,ISIGN)
              LATDIR = 'N'
              IF (ISIGN .eq. -1) LATDIR = 'S'
              CALL TODMSS(YLONT,LONDN,LONMN,SLONN,ISIGN)
              LONDIR = 'W'
              IF (ISIGN .eq. -1) LONDIR = 'E'
              LATS = IDNINT(SLATN*100000.D0)
              LONS = IDNINT(SLONN*100000.D0)
              WRITE(CARD(45:56),310) LATDN,LATMN,LATS,LATDIR
  310         FORMAT(I2,I2.2,I7.7,A1)
              WRITE(CARD(57:69),320) LONDN,LONMN,LONS,LONDIR
  320         FORMAT(I3,I2.2,I7.7,A1)
            ENDIF

*** Unrecognized blue book record

      ELSE
c           WRITE(LUOUT,500) TYPE
  500       FORMAT(' This software does not recognize an ',A4,/
     1       ' record.  The record will be copied to the new'/
     2       ' file without change.')
      ENDIF
      WRITE(I2,175) CARD
      GO TO 170
  600 CONTINUE
      RETURN

  700 write (*,'(/)') 
      write (*,*) "Failed to read card in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  701 write (*,'(/)') 
      write (*,*) "Failed to read card *12* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  702 write (*,'(/)') 
      write (*,*) "Failed to read *50,51,52* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  703 write (*,'(/)') 
      write (*,*) "Failed to read *60,61* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  704 write (*,'(/)') 
      write (*,*) "Failed to read *20* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  705 write (*,'(/)') 
      write (*,*) "Failed to read *22* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  706 write (*,'(/)') 
      write (*,*) "Failed to read *53,54* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  707 write (*,'(/)') 
      write (*,*) "Failed to read *30,32* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  708 write (*,'(/)') 
      write (*,*) "Failed to read *80* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  709 write (*,'(/)') 
      write (*,*) "Failed to read I4 of *80* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
***************************************************************
      SUBROUTINE UPBB5(MIN1,MIN2,OPT,IOPT)

*** Update 5-digit blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER     CARD*212
      CHARACTER     SUBCRD*80
C     CHARACTER     DATE*6
      CHARACTER     DATE8*8
      CHARACTER     TYPE*3
      CHARACTER     OPT*1,SIGN*1,JN*1,JW*1
      CHARACTER     LATDIR*1,LONDIR*1
      CHARACTER     PIDs*6         
      LOGICAL       TEST
      COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim) ,PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

  170 READ(I1,175,END=600,err=700,iostat=ios) CARD
      if (ios /= 0) goto 700
  175 FORMAT(A212)
      SUBCRD = CARD(1:80)
      TYPE = CARD(7:9) 
      IF(TYPE .EQ. '*A1' .OR.
     1      TYPE .EQ. '*AA' .OR.
     1      TYPE .EQ. '*10' .OR.
     1      TYPE .EQ. '*11' .OR.
     1      TYPE .EQ. '*12' .OR.
     1      TYPE .EQ. '*13' .OR.
     1      TYPE .EQ. '*21' .OR.
     1      TYPE .EQ. '*25' .OR.
     1      TYPE .EQ. '*26' .OR.
     1      TYPE .EQ. '*27' .OR.
     1      TYPE .EQ. '*28' .OR.
     1      TYPE .EQ. '*29' ) THEN
              CONTINUE
      ELSEIF(TYPE .EQ. '*31' .OR.
     1      TYPE .EQ. '*40' .OR.
     1      TYPE .EQ. '*41' .OR.
     1      TYPE .EQ. '*42' .OR.
     1      TYPE .EQ. '*45' .OR.
     1      TYPE .EQ. '*46' .OR.
     1      TYPE .EQ. '*47' .OR.
     1      TYPE .EQ. '*70' .OR.
     1      TYPE .EQ. '*81' .OR.
     1      TYPE .EQ. '*82' .OR.
     1      TYPE .EQ. '*83' .OR.
     1      TYPE .EQ. '*84' .OR.
     1      TYPE .EQ. '*85' .OR.
     1      TYPE .EQ. '*86' .OR.
     1      TYPE .EQ. '*90') THEN
          CONTINUE

*** Regular distance

      ELSEIF(TYPE .EQ. '*50' .OR.
     1       TYPE .EQ. '*51' .OR.
     1       TYPE .EQ. '*52') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4')THEN
C          READ(CARD,180) ISN,DATE,JSN,OBS
C 180      FORMAT(BZ, 9X,I5,T35,A6,T46,I5,T64,F9.4)
           READ(CARD,180,err=702,iostat=ios) ISN,JSN,OBS,DATE8
           if (ios /= 0) goto 702
  180      FORMAT(BZ, 9X,I5,T46,I5,T64,F9.4,T101,A8)
           CALL CHECK(ISN,JSN,SUBCRD,TEST)
           IF(TEST) THEN
             CALL TNFDAT(DATE8,MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)    
             IOBS = IDNINT((OBS+DS)*10000.D0)
             WRITE(CARD(64:72),190) IOBS
  190        FORMAT(I9)
             WRITE(CARD(174:181),191) DS
  191        FORMAT(F8.4)
           ENDIF
         ENDIF

*** Azimuth

      ELSEIF(TYPE .EQ. '*60' .OR. TYPE .EQ. '*61') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
C          READ(CARD,200) ISN,DATE,JSN,IDEG,MIN,SEC
C 200      FORMAT(BZ, 9X,I5,T40,A6,T51,I5,T64,I3,I2,F3.1)
           READ(CARD,200,err=703,iostat=ios) ISN,JSN,IDEG,MIN,
     &                        SEC,DATE8
           if (ios /= 0) goto 703
  200      FORMAT(BZ, 9X,I5,T51,I5,T64,I3,I2,F3.1,T101,A8)
           CALL CHECK(ISN,JSN,SUBCRD,TEST)
           IF(TEST) THEN
             CALL TNFDAT(DATE8,MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)     
             OBS = (DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
             OBS = OBS + DA
             IF(OBS .GE. TWOPI) THEN
               OBS = OBS - TWOPI
             ELSEIF(OBS .LT. 0.D0) THEN
               OBS = OBS + TWOPI
             ENDIF
             CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
             ISEC = IDNINT(SEC*10.D0)
             WRITE(CARD(64:71),210) IDEG,MIN,ISEC
  210        FORMAT(I3,I2.2,I3.3)
             IF(DA .GE. 0.0D0) THEN
               SIGN = ' '
             ELSE
               SIGN = '-'
             ENDIF
             ANGLE = DABS(DA)
             CALL TODMSS(ANGLE,IDEG,MIN,SEC,ISIGN)
             ISEC1 = SEC
             ISEC2 = IDNINT((SEC - DBLE(ISEC1)) * 100.D0)
             WRITE(CARD(152:162),211)SIGN,IDEG,MIN,ISEC1,ISEC2
  211        FORMAT(A1,I3.3,I2.2,I2.2,'.',I2.2)
           ENDIF
         ENDIF

*** Direction observation

      ELSEIF(TYPE .EQ. '*20') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
C           READ(CARD,220) ISN0,LIST0,DATE,JSN
C 220       FORMAT(BZ, 9X,I5,I2,T40,A6,T51,I5)
            READ(CARD,220,err=704,iostat=ios) ISN0,LIST0,JSN,DATE8
            if (ios /= 0) goto 704
  220       FORMAT(BZ, 9X,I5,I2,T51,I5,T101,A8)
            CALL CHECK(ISN0,JSN,SUBCRD,TEST) 
            CALL TNFDAT(DATE8,MINO)
            IF(TEST) THEN
              CALL DSDA(LOC(ISN0),LOC(JSN),MINO,MIN2,DS,DA0)     
            ELSE
              DA0 = 0.D0
            ENDIF
            CARD(132:142) = ' 0000000.00'
         ENDIF
      ELSEIF(TYPE .EQ. '*22') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
           READ(CARD,230,err=705,iostat=ios) ISN,LIST,JSN,IDEG,
     &                                       MIN,SEC
           if (ios /= 0) goto 705
  230      FORMAT(BZ, 9X,I5,I2,T51,I5,T64,I3,I2,F4.2)
           IF(LIST.NE.LIST0 .OR. ISN.NE.ISN0) THEN
             WRITE(LUOUT,240) CARD
  240           FORMAT(' Blue-book file has incorrect structure.'/
     1          ' A *22* record disagrees with its corresponding'/
     2          ' *20* record.  The record reads:'/A212)
             STOP
           ENDIF
           CALL CHECK(ISN,JSN,SUBCRD,TEST)
           IF(TEST) THEN
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)    
           ELSE
             DA = 0.D0
           ENDIF
           OBS=(DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
           OBS = OBS + DA - DA0
           IF(OBS .GT. TWOPI) THEN
             OBS = OBS - TWOPI
           ELSEIF(OBS .LT. 0.D0) THEN
             OBS = OBS + TWOPI
           ENDIF
           CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
           ISEC = IDNINT(SEC*100.D0)
           WRITE(CARD(64:72),250) IDEG,MIN,ISEC
  250      FORMAT(I3.3,I2.2,I4.4)
           IF((DA-DA0) .GE. 0.D0) THEN
             SIGN = ' '
           ELSE
             SIGN = '-'
           ENDIF
           ANGLE = DABS(DA - DA0)
           CALL TODMSS(ANGLE,IDEG,MIN,SEC,ISIGN)
           ISEC1 = SEC
           ISEC2 = IDNINT((SEC-DBLE(ISEC1))*100.D0)
           WRITE(CARD(132:142),251) SIGN,IDEG,MIN,ISEC1,ISEC2
  251      FORMAT(A1,I3.3,I2.2,I2.2,'.',I2.2)
         ENDIF

*** Long distance

      ELSEIF(TYPE .EQ. '*53' .OR.
     1       TYPE .EQ. '*54') THEN
         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
C          READ(CARD,260) ISN,DATE,JSN,OBS
C 260      FORMAT(BZ, 9X,I5,T35,A6,T46,I5,T64,F10.3)
           READ(CARD,260,err=706,iostat=ios) ISN,JSN,OBS,DATE8
           if (ios /= 0) goto 706
  260      FORMAT(BZ, 9X,I5,T46,I5,T64,F10.3,T101,A8)
           CALL CHECK(ISN,JSN,SUBCRD,TEST)
           IF(TEST) THEN
             CALL TNFDAT(DATE8,MINO)
             CALL DSDA(LOC(ISN),LOC(JSN),MINO,MIN2,DS,DA)     
             IOBS = IDNINT((OBS + DS)*1000.D0)
             WRITE(CARD(64:73),270) IOBS
  270        FORMAT(I10)
             CORR = IDNINT(DS*10000.D0)
             WRITE(CARD(174:181),271) CORR
  271        FORMAT(F8.4)
           ENDIF
         ENDIF
  
*** position record

      ELSEIF(TYPE .EQ. '*80') THEN
            IF(OPT .eq. '2' .or. OPT .eq. '4') THEN
              READ(CARD,300,err=708,iostat=ios) ISN,LATD,LATM,
     &                      SLAT,JN,LOND,LONM,SLON,JW
              if (ios /= 0) goto 708
  300         FORMAT(BZ,9X,I5,T45,2I2,F7.5,A1,I3,I2,F7.5,A1)
              YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
              YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
              IF(JN .EQ. 'S') YLAT = -YLAT
              IF(JW .EQ. 'E') YLON = -YLON
C             READ(I4, REC = LOC(ISN)) RLAT, RLON, IDG
C             CALL NWCORD(YLAT,YLON,HT(ISN),
C    1          RLAT,RLON,IDG,MIN1,MIN2,
C    1          YLATT,YLONT,EHTNEW,DN,DE,DU,IOPT)
          READ(I4, REC = LOC(ISN),err=709,iostat=ios) RLAT, RLON, 
     &                      VN, VE, VU
              if (ios /= 0) goto 709
          call NEWCOR( ylat, ylon, ht(isn), min1, min2,
     1           ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
              CALL TODMSS(YLATT,LATDN,LATMN,SLATN,ISIGN)
              LATDIR = 'N'
          IF (ISIGN .eq. -1) LATDIR = 'S'
              CALL TODMSS(YLONT,LONDN,LONMN,SLONN,ISIGN)
              LONDIR = 'W'
          IF (ISIGN .eq. -1) LONDIR = 'E'
              LATS = IDNINT(SLATN*100000.D0)
              LONS = IDNINT(SLONN*100000.D0)
              WRITE(CARD(45:56),310) LATDN,LATMN,LATS,LATDIR
  310         FORMAT(I2,I2.2,I7.7,A1)
              WRITE(CARD(57:69),320) LONDN,LONMN,LONS,LONDIR
  320         FORMAT(I3,I2.2,I7.7,A1)
            ENDIF

*** Unrecognized blue book record

      ELSE
            WRITE(LUOUT,500) TYPE
  500       FORMAT(' This software does not recognize an ',A3,/
     1       ' record.  The record will be copied to the new'/
     2       ' file without change.')
      ENDIF
      WRITE(I2,175) CARD
      GO TO 170
  600 CONTINUE
      RETURN

  700 write (*,'(/)') 
      write (*,*) "Failed to read card in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  701 write (*,'(/)') 
      write (*,*) "Failed to read card *12* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  702 write (*,'(/)') 
      write (*,*) "Failed to read *50,51,52* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  703 write (*,'(/)') 
      write (*,*) "Failed to read *60,61* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  704 write (*,'(/)') 
      write (*,*) "Failed to read *20* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  705 write (*,'(/)') 
      write (*,*) "Failed to read *22* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  706 write (*,'(/)') 
      write (*,*) "Failed to read *53,54* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  707 write (*,'(/)') 
      write (*,*) "Failed to read *30,32* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  708 write (*,'(/)') 
      write (*,*) "Failed to read *80* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  709 write (*,'(/)') 
      write (*,*) "Failed to read I4 of *80* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
***************************************************************
      SUBROUTINE CHECK(ISN,JSN,CARD,TEST)

*** Check if stations have *80* records

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER    CARD*80
      CHARACTER    PIDs*6
      LOGICAL TEST
      COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      TEST = .TRUE.

      IF(LOC(ISN) .eq. 0) THEN
          WRITE(LUOUT,100) ISN, CARD
  100     FORMAT(' No *80* record for SSN = ',I4/A80/)
          TEST = .FALSE.
      ENDIF

      IF(LOC(JSN) .eq. 0) THEN
          WRITE(LUOUT,100) JSN, CARD
          TEST = .FALSE.
      ENDIF

      RETURN
      END
C*************************************************
      SUBROUTINE UPGFI4(DATE2, MIN2, IOPT, KOPT,
     *     MONTH, IDAY, IYEAR)

*** Update G-FILE of blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
c     CHARACTER*80 CARD
      CHARACTER     CARD*120
      CHARACTER     TYPE*1
      CHARACTER     NRF*2
      CHARACTER     ZT*2
      CHARACTER     CHAR14*14
      LOGICAL       TEST
      LOGICAL Is_inp_NAD83,Is_out_NAD83,Is_out1_NAD83
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      ZT = 'ZT'

*** Obtain blue-book reference frame identifier
*** corresponding to KOPT  (output frame for vectors)
      IF (KOPT .NE. -1) THEN
         CALL RFCON1(KOPT, NBBREF)
         WRITE(NRF,10) NBBREF
   10    FORMAT(I2.2)
      ENDIF

   90 READ(I1,100,END=200,err=300,iostat=ios) CARD
      if (ios /= 0) goto 300
c 100 FORMAT(A80)
  100 FORMAT(A120)
      TYPE = CARD(1:1)
      IF(TYPE .eq. 'A') THEN
         IF (CARD(79:80) .EQ. ZT) THEN
            WRITE(LUOUT, 103)
  103       FORMAT(
     1     ' *****************************************************'/
     1     ' * ERROR: The input GFILE contains the letters, ZT,  *'/
     1     ' * in columns 79-80 of its A-record.  This indicates *'/
     1     ' * that the GPS vectors in this file have already    *'/
     1     ' * been updated to a common date.  This software     *'/
     1     ' * will not further modify the GPS vectors.          *'/
     1     ' *****************************************************'/)
            RETURN
         ENDIF
         WRITE(CARD(79:80),101) ZT
  101    FORMAT(A2)
         WRITE(CARD(4:11), 102) IYEAR, MONTH, IDAY
  102    FORMAT(I4, I2.2, I2.2)
         WRITE(CARD(12:19), 102) IYEAR, MONTH, IDAY
                  
      ELSEIF(TYPE .eq. 'B') THEN
         READ(CARD,110,err=301,iostat=ios)IYEAR1,MONTH1,IDAY1,
     1        IYEAR2,MONTH2,IDAY2,IBBREF
         if (ios /= 0) goto 301
  110    FORMAT(1X,I4,I2,I2,4X,I4,I2,I2,30X,I2)
C        CALL TOTIME(IYEAR1,MONTH1,IDAY1,MINO1)
C        CALL TOTIME(IYEAR1, 1, 1, MIN00)
C        DECYR1 = DBLE(IYEAR1) + DBLE(MINO1 - MIN00)/525600.D0
C        CALL TOTIME(IYEAR2,MONTH2,IDAY2,MINO2)
C        CALL TOTIME(IYEAR2, 1, 1, MIN00)
C        DECYR2 = DBLE(IYEAR2) + DBLE(MINO2 - MIN00)/525600.D0

         CALL IYMDMJ(IYEAR1,MONTH1,IDAY1,MJD1)
         MINO1 = MJD1 * 24 * 60
         CALL IYMDMJ(IYEAR1, 1, 1, MJD0)
         DECYR1 = DBLE(IYEAR1) + DBLE(MJD1 - MJD0)/365.D0
         CALL IYMDMJ(IYEAR2,MONTH2,IDAY2, MJD2)
         MINO2 = MJD2 * 24 * 60
         CALL IYMDMJ(IYEAR2, 1, 1, MJD0)
         DECYR2 = DBLE(IYEAR2) + DBLE(MJD2 - MJD0)/365.D0
         MINO = (MINO1 + MINO2) / 2
         DECYR = (DECYR1 + DECYR2) / 2.D0
         CALL RFCON(IBBREF, JREF)         !JREF is the current frame or frame of input
         IF (KOPT .NE. -1) THEN
            CARD(52:53) = NRF
         ENDIF
         Is_inp_NAD83 = (JREF == 1)
         Is_out_NAD83 = (IOPT == 1)       !IOPT is the frame of the positions
         Is_out1_NAD83 = (KOPT == 1)       !IOPT is the frame of the positions
      ELSEIF(TYPE .eq. 'C') THEN
         READ(CARD,120,err=302,iostat=ios)ISN,JSN,DX,DY,DZ
         if (ios /= 0) goto 302
  120    FORMAT(BZ,1X,2I4,F11.4,5X,F11.4,5X,F11.4)
   
         CALL CHECK(ISN, JSN, CARD, TEST)
         IF (TEST) THEN
            CALL DDXYZ(ISN, JSN,                                 
     1              MINO, MIN2, DDX, DDY, DDZ)
*** Convert vector from JREF to IOPT frame                        
            if (Is_inp_NAD83 .or. Is_out_NAD83) then
              call TRAVEC (DX, DY, DZ, DECYR, JREF, IOPT)
            else
              call TRAVEC_IERS (DX, DY, DZ, DECYR, JREF, IOPT)
            endif
            DX = DX + DDX
            DY = DY + DDY
            DZ = DZ + DDZ

*** Convert GPS vector from IOPT to KOPT reference frame
            IF (KOPT .NE. -1) THEN
              if (Is_out_NAD83 .or. Is_out1_NAD83) then
                CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, KOPT)
              else
                CALL TRAVEC_IERS(DX, DY, DZ, DATE2, IOPT, KOPT)
              endif
            ELSE
              if (Is_inp_NAD83 .or. Is_out_NAD83) then
                CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, JREF)
              else
                CALL TRAVEC_IERS(DX, DY, DZ, DATE2, IOPT, JREF)
              endif
            ENDIF

*** Rewrite GPS observational record
            CALL TOCHAR(DX,CHAR14)          
            CARD(10:20) = CHAR14(3:13)          

            CALL TOCHAR(DY,CHAR14)
            CARD(26:36) = CHAR14(3:13)    

            CALL TOCHAR(DZ,CHAR14)
            CARD(42:52) = CHAR14(3:13)    
         ENDIF

      ELSEIF(TYPE .eq. 'F') THEN
         READ(CARD,140,err=303,iostat=ios)ISN,JSN,DX,DY,DZ
         if (ios /= 0) goto 303
  140    FORMAT(BZ,1X,2I4,F13.4,5X,F13.4,5X,F13.4)

         CALL CHECK(ISN, JSN, CARD, TEST)
         IF (TEST) THEN
            CALL DDXYZ(ISN, JSN,                               
     1              MINO, MIN2, DDX, DDY, DDZ)

            if (Is_inp_NAD83 .or. Is_out_NAD83) then
              call TRAVEC( DX, DY, DZ, DECYR, JREF, IOPT) 
            else
              call TRAVEC_IERS( DX, DY, DZ, DECYR, JREF, IOPT) 
            endif

            DX = DX + DDX
            DY = DY + DDY
            DZ = DZ + DDZ

*** Convert GPS vector to output reference frame
            IF (KOPT .NE. -1) THEN
              if (Is_out_NAD83 .or. Is_out1_NAD83) then
                CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, KOPT)
              else
                CALL TRAVEC_IERS (DX, DY, DZ, DATE2, IOPT, KOPT)
              endif
            ELSE
              if (Is_out_NAD83 .or. Is_out1_NAD83) then
                CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, JREF)
              else
                CALL TRAVEC_IERS (DX, DY, DZ, DATE2, IOPT, JREF)
              endif
            ENDIF

*** Rewrite GPS observational record
            CALL TOCHAR(DX,CHAR14)
            CARD(10:22) = CHAR14(1:13)

            CALL TOCHAR(DY,CHAR14)
            CARD(28:40) = CHAR14(1:13)           

            CALL TOCHAR(DZ,CHAR14)
            CARD(46:58) = CHAR14(1:13)
         ENDIF
      ENDIF
      WRITE(I2,100) CARD
      GO TO 90
  200 CONTINUE
      RETURN

  300 write (*,'(/)') 
      write (*,*) "Failed to read gfile in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  301 write (*,'(/)') 
      write (*,*) "Failed to read B card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  302 write (*,'(/)') 
      write (*,*) "Failed to read C card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  303 write (*,'(/)') 
      write (*,*) "Failed to read F card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*********************************************************
      SUBROUTINE UPGFI5(DATE2, MIN2, IOPT, KOPT,
     *                  MONTH, IDAY, IYEAR)

*** Update G-FILE of 5-digit blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER     CARD*180
      CHARACTER     SUBCRD*80
      CHARACTER     TYPE*1
      CHARACTER     NRF*2
      CHARACTER     ZT*2
      CHARACTER     CHAR14*14
      LOGICAL       TEST
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      ZT = 'ZT'

*** Obtain blue-book reference frame identifier
*** corresponding to IOPT
      IF (KOPT .NE. -1) THEN
         CALL RFCON1(KOPT, NBBREF)
         WRITE(NRF, 10) NBBREF
   10    FORMAT(I2.2)
      ENDIF

   90 READ(I1,100,END=200,err=300,iostat=ios) CARD
      if (ios /= 0) goto 300
  100 FORMAT(A180)
      TYPE = CARD(1:1)

      IF (TYPE .eq. 'A') THEN
         IF (CARD(79:80) .EQ. ZT) THEN
            WRITE(LUOUT, 103)
  103       FORMAT(
     1     ' *****************************************************'/
     1     ' * ERROR: The input GFILE contains the letters, ZT,  *'/
     1     ' * in columns 79-80 of its A-record.  This indicates *'/
     1     ' * that the GPS vectors in this file have already    *'/
     1     ' * been updated to a common date.  This software     *'/
     1     ' * will not further modify the GPS vectors.          *'/
     1     ' *****************************************************'/)
            RETURN
         ENDIF
         WRITE(CARD(79:80), 101) ZT
  101    FORMAT(A2)
         WRITE(CARD(4:11), 102) IYEAR, MONTH, IDAY
  102    FORMAT(I4, I2.2, I2.2)
         WRITE(CARD(12:19), 102) IYEAR, MONTH, IDAY
                            
      ELSEIF(TYPE .eq. 'B') THEN
         READ(CARD,110,err=301,iostat=ios)IYEAR1,MONTH1,IDAY1,
     1        IYEAR2,MONTH2,IDAY2,IBBREF
         if (ios /= 0) goto 301
  110    FORMAT(1X,I4,I2,I2,4X,I4,I2,I2,30X,I2)
C        CALL TOTIME(IYEAR1,MONTH1,IDAY1,MINO1)
C	 CALL TOTIME(IYEAR1, 1, 1, MIN00)
C	 DECYR1 = DBLE(IYEAR1) + DBLE(MINO1 - MIN00)/525600.D0
C        CALL TOTIME(IYEAR2,MONTH2,IDAY2,MINO2)
C 	 CALL TOTIME(IYEAR2, 1, 1, MIN00)
         CALL IYMDMJ(IYEAR1,MONTH1,IDAY1,MJD1)
         MINO1 = MJD1 * 24 * 60
         CALL IYMDMJ(IYEAR1, 1, 1, MJD0)
         DECYR1 = DBLE(IYEAR1) + DBLE(MJD1 - MJD0)/365.D0
         CALL IYMDMJ(IYEAR2, MONTH2,IDAY2, MJD2)
         MINO2 = MJD2 * 24 * 60
         CALL IYMDMJ(IYEAR2, 1, 1, MJD0)
         DECYR2 = DBLE(IYEAR2) + DBLE(MJD2 - MJD0)/365.D0
C	 DECYR2 = DBLE(IYEAR2) + DBLE(MINO2 - MIN00)/525600.D0
         MINO = (MINO1 + MINO2) / 2
         DECYR = (DECYR1 + DECYR2) / 2.D0
         CALL RFCON(IBBREF, JREF)
         IF (KOPT .NE. -1) THEN
            CARD(52:53) = NRF
         ENDIF

      ELSEIF(TYPE .eq. 'F') THEN
         READ(CARD,140,err=303,iostat=ios)ISN,DX,DY,DZ,JSN
         if (ios /= 0) goto 303
  140    FORMAT(BZ,1X,I5,3X,F13.4,5X,F13.4,5X,F13.4,T81,I5)

         SUBCRD = CARD(1:80)
         CALL CHECK(ISN, JSN, SUBCRD, TEST)
         IF (TEST) THEN
            CALL DDXYZ(ISN, JSN,                                 
     1              MINO,MIN2,DDX, DDY, DDZ)

            call TRAVEC( DX, DY, DZ, DECYR, JREF, IOPT)

            DX = DX + DDX
            DY = DY + DDY
            DZ = DZ + DDZ

*** Convert GPS vector to outout reference frame
            IF (KOPT .NE. -1) THEN
              CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, KOPT)
            ELSE
              CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, JREF)
            ENDIF

*** Rewrite GPS observational record
            CALL TOCHAR(DX,CHAR14)
            CARD(10:22) = CHAR14(1:13)

            CALL TOCHAR(DY,CHAR14)
            CARD(28:40) = CHAR14(1:13)           

            CALL TOCHAR(DZ,CHAR14)
            CARD(46:58) = CHAR14(1:13)
 
            ISHIFT = IDNINT(DDX * 10000.D0)
            WRITE(CARD(115:121),150) ISHIFT
  150       FORMAT(I7)

            ISHIFT = IDNINT(DDY * 10000.D0)
            WRITE(CARD(125:131),150) ISHIFT

            ISHIFT = IDNINT(DDZ * 10000.D0)
            WRITE(CARD(135:141),150) ISHIFT
         ENDIF

      ENDIF
      WRITE(I2,100) CARD
      GO TO 90
  200 CONTINUE
      RETURN

  300 write (*,'(/)') 
      write (*,*) "Failed to read gfile in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  301 write (*,'(/)') 
      write (*,*) "Failed to read B card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  302 write (*,'(/)') 
      write (*,*) "Failed to read C card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  303 write (*,'(/)') 
      write (*,*) "Failed to read F card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*********************************************************
      SUBROUTINE TODMSS(val,id,im,s,isign)
 
*** convert position radians to deg,min,sec
*** range is [-twopi to +twopi]
 
      implicit double precision(a-h,o-z)
      IMPLICIT INTEGER*4 (I-N)
      common/CONST/A,F,E2,EP2,AF,PI,TWOPI,RHOSEC
 
    1 if(val.gt.twopi) then
        val=val-twopi
        go to 1
      endif
 
    2 if(val.lt.-twopi) then
        val=val+twopi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*RHOSEC/3600.D0)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
*** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
      endif
 
      return
      end
*********************************************************
      SUBROUTINE SETRF

*** Specify arrays that may be used to convert
*** a reference frame identifier in the blue book
*** to a reference frame identifier in HTDP and back

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 16 )
      COMMON /REFCON/ IRFCON(36), JRFCON(numref)

*** From blue book identifier to HTDP indentifier
*** WGS 72 Precise
c     IRFCON(1) = 10
      IRFCON(1) = 1
C HTDP no longer supports WGS 72. Hence, if a BlueBook
C file contains WGS 72 coordinates, HTDP treats these
C coordinates as if they were NAD 83(2011)coordinates.

*** WGS 84 (orig) Precise (set  equal to NAD 83)
      IRFCON(2) = 1

*** WGS 72 Broadcast
c     IRFCON(3) = 10
      IRFCON(3) = 1
C HTDP no longer supports WGS 72. Hence, if a BlueBook
C file contains WGS 72 coordinates, HTDP treats these
C coordinates as if they were NAD 83(2011)coordinates.

*** WGS 84 (orig) Broadcast (set equal to NAD 83)
      IRFCON(4) = 1

*** ITRF89
      IRFCON(5) = 3

*** PNEOS 90 or NEOS 91.25 (set equal to ITRF90)
      IRFCON(6) = 4

*** NEOS 90 (set equal to ITRF90)
      IRFCON(7) = 4

*** ITRF91
      IRFCON(8) = 5

*** SIO/MIT 92.57 (set equal to ITRF91)
      IRFCON(9) = 5

*** ITRF91
      IRFCON(10) = 5

*** ITRF92
      IRFCON(11) = 6

*** ITRF93
      IRFCON(12) = 7

*** WGS 84 (G730) Precise (set equal to ITRF91)
      IRFCON(13) = 5

*** WGS 84 (G730) Broadcast (set equal to ITRF91)
      IRFCON(14) = 5

*** ITRF94
      IRFCON(15) = 8

*** WGS 84 (G873) Precise  (set equal to ITRF94)
      IRFCON(16) = 8

*** WGS 84 (G873) Broadcast (set equal to ITRF94)
      IRFCON(17) = 8

*** ITRF96
      IRFCON(18) = 8

*** ITRF97
      IRFCON(19) = 9

*** IGS97
      IRFCON(20) = 9

*** ITRF00
      IRFCON(21) = 11

*** IGS00
      IRFCON(22) = 11

*** WGS 84 (G1150)
      IRFCON(23) = 11

*** IGb00
      IRFCON(24) = 11

*** ITRF2005
      IRFCON(25) = 14

*** IGS05
      IRFCON(26) = 14

*** IGS08
      IRFCON(27) = 15

*** IGB08
      IRFCON(28) = 15

*** ITRF2008
      IRFCON(29) = 15

*** WGS84 (G1674)
      IRFCON(30) = 15

*** WGS84 (G1762)
      IRFCON(31) = 15

*** ITRF2014
      IRFCON(32) = 16

*** IGB14
      IRFCON(33) = 16

*** NAD83 (2011/2007/CORS96/FBN/HARN)
      IRFCON(34) = 1

*** NAD83 (PA11)
      IRFCON(35) = 12

*** NAD83 (MA11)
      IRFCON(36) = 13

*** From HTDP identifier to blue book identifier
*** NAD 83 (set equal to WGS 84 (transit))
c     JRFCON(1) = 2
      JRFCON(1) = 34

*** ITRF88 (set equal to ITRF89)
      JRFCON(2) = 5

*** ITRF89
      JRFCON(3) = 5

*** ITRF90 (set equal to NEOS 90)
      JRFCON(4) = 7

*** ITRF91
      JRFCON(5) = 8

*** ITRF92
      JRFCON(6) = 11

*** ITRF93
      JRFCON(7) = 12

*** ITRF96 (= ITRF94)
      JRFCON(8) = 18

*** ITRF97
      JRFCON(9) = 19

*** WGS 72
      JRFCON(10) = 1

*** ITRF00
      JRFCON(11) = 21

*** NAD 83 (PACP00) or NAD 83 (PA11)
c     JRFCON(12) = 0
c     JRFCON(12) = 2
      JRFCON(12) = 35
C  Michael Dennis requested that HTDP identify
C  NAD 83 (PA11) coordinates as WGS 84(transit)
C  coordinates in the Bluebook context

*** NAD 83 (MARP00) or NAD 83 (MA11)
c     JRFCON(13) = 0
c     JRFCON(13) = 2 
      JRFCON(13) = 36
C  Michael Dennis requested that HTDP identify
C  NAD 83 (MA11) coordinates as WGS 84(transit)
C  coordinates in the Bluebook context

*** ITRF2005 or IGS05
      JRFCON(14) = 26

*** ITRF2008 or IGS08
      JRFCON(15) = 27

*** IGB08
      JRFCON(15) = 28

*** ITRF2014 or IGS14
      JRFCON(16) = 33

      RETURN
      END
***************************************************
      SUBROUTINE RFCON(IBBREF, JREF)

*** Convert reference frame identifier from
*** system used in the blue-book to the
*** system used in HTDP

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 16 )
      COMMON /REFCON/ IRFCON(36), JRFCON(numref)

      IF (1 .LE. IBBREF .AND. IBBREF .LE. 36) THEN
          JREF = IRFCON(IBBREF)
      ELSE
          WRITE(6, 10) IBBREF
   10     FORMAT(' Improper reference frame identifier (=',
     1      I4, ')' /
     1      ' appearing in B-record of the G-FILE')
          STOP
       ENDIF

       RETURN
       END
*******************************************************
      SUBROUTINE RFCON1(JREF, IBBREF)

*** Convert reference frame identifier from
*** system used in HTDP to the  system
*** used in the blue-book

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 16 )
      COMMON /REFCON/ IRFCON(36), JRFCON(numref)

      IF (JREF .EQ. 0) THEN
        I = 1
      ELSE
        I = JREF
      ENDIF

      IF(1. LE. I .AND. I .LE. numref) THEN
        IBBREF = JRFCON(I)
        IF ( IBBREF .eq. 0) THEN
          write(6,5) JREF
    5     format(' ERROR: The BlueBook does not recognize '/
     *           'this reference frame, HTDP ID = ',I2)
          stop
        ENDIF

      ELSE
          WRITE(6, 10) JREF
   10     FORMAT(' Improper reference frame identifier (=',
     *       I4, ')' / 'appearing in routine RFCON1')
          STOP
      ENDIF

      RETURN
      END
*****************************************************************

      subroutine TRAVEC(dxi, dyi, dzi, date, jopt1, jopt2)

*** Transform GPS or other vector from one reference frame to another
*** for the given date                                     
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (dxi, dyi, dzi) --> (input) components of input vector in meters
***                 --> (output) components of transformed vector in meters

*** date --> (input) time (decimal years) to which the input & output
***          vectors correspond

*** jopt1 --> (input) specifier of reference frame for input vector
*** jopt2 --> (input) specifier of reference frame for output vector

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

*** Transform input vector to ITRF94 reference frame
      if (jopt1 .eq. 0) then
         iopt = 1
      else
         iopt = jopt1
      endif

      dtime = date - refepc(iopt)
      rotnx  = -(rx(iopt) + drx(iopt)*dtime)
      rotny  = -(ry(iopt) + dry(iopt)*dtime)
      rotnz  = -(rz(iopt) + drz(iopt)*dtime)
      ds     = 1.d0 - (scale(iopt) + dscale(iopt)*dtime)

      dxt = + ds*dxi + rotnz*dyi - rotny*dzi
      dyt = - rotnz*dxi + ds*dyi + rotnx*dzi
      dzt = + rotny*dxi - rotnx*dyi + ds*dzi

*** Transform ITRF94 vector to new reference frame
      if (jopt2 .eq. 0) then
         iopt = 1
      else
         iopt = jopt2
      endif

      dtime = date - refepc(iopt)
      rotnx  = rx(iopt) + drx(iopt)*dtime
      rotny  = ry(iopt) + dry(iopt)*dtime
      rotnz  = rz(iopt) + drz(iopt)*dtime
      ds     = 1.d0 + scale(iopt) + dscale(iopt)*dtime

      dxi = + ds*dxt + rotnz*dyt - rotny*dzt
      dyi = - rotnz*dxt + ds*dyt + rotnx*dzt
      dzi = + rotny*dxt - rotnx*dyt + ds*dzt

      return
      end

***********************************************************

      subroutine TRAVEC_IERS(dxi, dyi, dzi, date, jopt1, jopt2)

*** Transform GPS or other vector from one reference frame to another
*** for the given date, using the IERS 96-97 transformation parameters                                     
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

*** (dxi, dyi, dzi) --> (input) components of input vector in meters
***                 --> (output) components of transformed vector in meters

*** date --> (input) time (decimal years) to which the input & output
***          vectors correspond

*** jopt1 --> (input) specifier of reference frame for input vector
*** jopt2 --> (input) specifier of reference frame for output vector

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)


*** Transform input vector to ITRF94 reference frame
      if (jopt1 .eq. 0) then
         iopt = 1
      else
         iopt = jopt1
      endif

      dtime = date - refepc1(iopt)
      rotnx  = -(rx1(iopt) + drx1(iopt)*dtime)
      rotny  = -(ry1(iopt) + dry1(iopt)*dtime)
      rotnz  = -(rz1(iopt) + drz1(iopt)*dtime)
      ds     = 1.d0 - (scale1(iopt) + dscale1(iopt)*dtime)

      dxt = + ds*dxi + rotnz*dyi - rotny*dzi
      dyt = - rotnz*dxi + ds*dyi + rotnx*dzi
      dzt = + rotny*dxi - rotnx*dyi + ds*dzi

*** Transform ITRF94 vector to new reference frame
      if (jopt2 .eq. 0) then
         iopt = 1
      else
         iopt = jopt2
      endif

      dtime = date - refepc1(iopt)
      rotnx  = rx1(iopt) + drx1(iopt)*dtime
      rotny  = ry1(iopt) + dry1(iopt)*dtime
      rotnz  = rz1(iopt) + drz1(iopt)*dtime
      ds     = 1.d0 + scale1(iopt) + dscale1(iopt)*dtime

      dxi = + ds*dxt + rotnz*dyt - rotny*dzt
      dyi = - rotnz*dxt + ds*dyt + rotnx*dzt
      dzi = + rotny*dxt - rotnx*dyt + ds*dzt

      return
      end

***********************************************************
      SUBROUTINE GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1     LONDIR, NAME, X,Y,Z,XLAT,XLON,EHT)
   
*** Interactively obtain name and coordinates for a point.
*** Output longitude will be positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER    NAME*24
      CHARACTER    COPT*1,LATDIR*1,LONDIR*1
      LOGICAL      FRMXYZ
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

      WRITE(LUOUT,100)
  100 FORMAT(' Enter name for point (24 character max).  ')
      READ(LUIN,105,err=200,iostat=ios) NAME
      if (ios /= 0) goto 200
  105 FORMAT(A24)

  110 WRITE(LUOUT,111)
  111 FORMAT(' How do you wish to specify positional coordinates:'/
     1       '     1...geodetic latitude, longitude, ellipsoid height'/
     2       '     2...Cartesian (X,Y,Z) coordinates.  ')
      READ(LUIN,'(A1)',err=201,iostat=ios) COPT
      if (ios /= 0) goto 201

      IF(COPT .EQ. '1') THEN
        WRITE(LUOUT,115)
  115   FORMAT(
     1    ' Enter latitude degrees-minutes-seconds in free format'/,
     2    ' with north being positive. For example,    35,17,28.3  '/
     2    ' For a point in the southern hemisphere, enter a minus sign'/
     3    ' before each value. For example, -35,-17,-28.3')
        READ(LUIN,*,err=202,iostat=ios) LATD, LATM, SLAT
        if (ios /= 0) goto 202
        WRITE(LUOUT,120)
  120   FORMAT(
     1    ' Enter longitude degrees-minutes-seconds in free format'/,
     2    ' with west being positive.  To express a longitude measured'/
     3    ' eastward, enter a minus sign before each value.')
        READ(LUIN,*,err=203,iostat=ios) LOND, LONM, SLON
          if (ios /= 0) goto 203
          WRITE(LUOUT, 125)
  125     FORMAT(
     1    ' Enter ellipsoid height in meters. (Note that'/,
     1    ' predicted motions are independent of this height.)  ')
          READ(LUIN,*,err=204,iostat=ios) EHT
          if (ios /= 0) goto 204
          XLAT =  (DBLE((LATD*60 + LATM)*60) + SLAT)/RHOSEC
          LATDIR = 'N'
          IF (XLAT .lt. 0.0D0) then
            LATD = - LATD
            LATM = - LATM
            SLAT = - SLAT
            LATDIR = 'S'
          ENDIF
          XLON = (DBLE((LOND*60 + LONM)*60) + SLON)/RHOSEC
          ELON = -XLON
          CALL TOXYZ(XLAT,ELON,EHT,X,Y,Z)
          LONDIR = 'W'
          IF (XLON .lt. 0.0D0) then
            LOND   = -LOND
            LONM   = -LONM
            SLON   = -SLON
            LONDIR = 'E'
          ENDIF

      ELSEIF(COPT .EQ. '2') THEN
          WRITE(LUOUT,130)
  130     FORMAT(' Enter X coordinate in meters.  ')
          READ(LUIN,*,err=205,iostat=ios) X
          if (ios /= 0) goto 205
          WRITE(LUOUT,140)
  140     FORMAT(' Enter Y coordinate in meters.  ')
          READ(LUIN,*,err=206,iostat=ios) Y
          if (ios /= 0) goto 206
          WRITE(LUOUT,150)
  150     FORMAT(' Enter Z coordinate in meters.  ')
          READ(LUIN,*,err=207,iostat=ios) Z
          if (ios /= 0) goto 207
          IF(.NOT.FRMXYZ(X,Y,Z,XLAT,XLON,EHT)) STOP 666
          XLON = -XLON
          IF(XLON .LT. 0.0D0) XLON = XLON + TWOPI
          CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
          LATDIR = 'N'
          IF (ISIGN .eq. -1) LATDIR = 'S'
          CALL TODMSS(XLON,LOND,LONM,SLON,ISIGN)
          LONDIR = 'W'
          IF (ISIGN .eq. -1) LONDIR = 'E'
      ELSE
          WRITE(LUOUT,*) ' Improper response -- try again.  '
          GO TO 110
      ENDIF

      RETURN

  200 write (*,'(/)') 
      write (*,*) "Failed to read point name: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  201 write (*,'(/)') 
      write (*,*) "Failed to read Coord. form option:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  202 write (*,'(/)') 
      write (*,*) "Failed to read latitude:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  203 write (*,'(/)') 
      write (*,*) "Failed to read longitude:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  204 write (*,'(/)') 
      write (*,*) "Failed to read ellipsoidal height:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  205 write (*,'(/)') 
      write (*,*) "Failed to read the X coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  206 write (*,'(/)') 
      write (*,*) "Failed to read the Y coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  207 write (*,'(/)') 
      write (*,*) "Failed to read the Z coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*****************************************************************
      SUBROUTINE GETVLY(GLAT, GLON, VX, VY, VZ,
     1      VNORTH, VEAST, VUP, VOPT, IFORM)

*** Interactively obtain velocity for a point

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER   VOPT*1
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

  200 CONTINUE
      IF (IFORM .eq. 210) then
         WRITE(LUOUT, 210)
  210    FORMAT(' What is the velocity of this site relative to '/ 
     1       ' the input reference frame:'/ 
     1       '     0...use velocity predicted by this software.'/ 
     1       '     1...user will specify the north-east-up'/ 
     1       '         components of site velocity.'/ 
     1       '     2...user will specify the global x-y-z'/
     1       '         components of site velocity.' )
      else
        write(luout, 211)
  211   FORMAT(' How do you wish to specify the velocity: '/
     1       '     1...north-east-up components.'/
     1       '     2...global x-y-z components.' )
      endif

      READ(LUIN, '(A1)',err=300,iostat=ios) VOPT
      if (ios /= 0) goto 300

      IF (VOPT .EQ. '0') THEN
        CONTINUE
      ELSEIF (VOPT .EQ. '1') THEN
        WRITE(LUOUT, 220)
  220   FORMAT( ' Enter north-south component of velocity in mm/yr'/
     1          ' with north being positive and south being negative')
        READ(LUIN, *,err=301,iostat=ios) VNORTH
        if (ios /= 0) goto 301
        WRITE(LUOUT, 221)
  221   FORMAT( ' Enter east-west component of velocity in mm/yr'/
     1           ' with east being positive and west being negative')
        READ(LUIN, *,err=302,iostat=ios) VEAST
        if (ios /= 0) goto 302
        WRITE(LUOUT, 222)
  222   FORMAT( ' Enter vertical component of velocity in mm/yr'/
     1          ' with up being positive and down being negative'/
     1           ' or enter 0.0 if unknown.')
        READ(LUIN, *,err=303,iostat=ios) VUP
        if (ios /= 0) goto 303
        CALL TOVXYZ( GLAT, GLON, VNORTH, VEAST, VUP, VX, VY, VZ)
      ELSEIF (VOPT .EQ. '2') THEN
        WRITE(LUOUT, *) ' Enter x-component of velocity in mm/yr.'
        READ(LUIN, *,err=304,iostat=ios) VX
        if (ios /= 0) goto 304
        WRITE(LUOUT, *) ' Enter y-component of velocity in mm/yr.'
        READ(LUIN, *,err=305,iostat=ios) VY
        if (ios /= 0) goto 305
        WRITE(LUOUT, *) ' Enter z-component of velocity in mm/yr.'
        READ(LUIN, *,err=306,iostat=ios) VZ
        if (ios /= 0) goto 306
        CALL TOVNEU(GLAT, GLON, VX, VY, VZ, VNORTH, VEAST, VUP)
      ELSE
        WRITE(LUOUT, *) 'Improper response -- try again. '
        GO TO 200
      ENDIF
      RETURN

 300  write (*,'(/)') 
      write (*,*) "Failed to read form of velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 301  write (*,'(/)') 
      write (*,*) "Failed to read North velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 302  write (*,'(/)')
      write (*,*) "Failed to read East velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 303  write (*,'(/)')
      write (*,*) "Failed to read Up velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 304  write (*,'(/)')
      write (*,*) "Failed to read X velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 305  write (*,'(/)')
      write (*,*) "Failed to read Y velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 306  write (*,'(/)')
      write (*,*) "Failed to read Z velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END

************************************************************************
      SUBROUTINE COMPSN(YLATT,YLONT,HTT,YLAT,YLON,HT,
     1                  MIN,VN, VE, VU)

*** Compute the position of a point at specified time
*** Upon input VN, VE, and VU are in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NDLOC = 2195)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /QPARM/ STRIKE(NDLOC), HL(NDLOC), EQLAT(NDLOC),
     1          EQLON(NDLOC), SSLIP(NDLOC), DSLIP(NDLOC),
     1          DIP(NDLOC), DEPTH(NDLOC), WIDTH(NDLOC),
     1               EQLATR(50),EQLONR(50),EQRAD(50),
     1               ITEQK(50),NLOC(50),NFP(50),NUMEQ
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

** Compute the contribution due to constant velocity
         DTIME = DBLE(MIN - ITREF) / 525960.D0
         CALL RADR8T(YLAT,VN,VE,VNR,VER)
         YLATT = YLAT + VNR*DTIME
         YLONT = YLON - VER*DTIME
         HTT   = HT + ((VU * DTIME) /1000.D0)
c        write (*,*) "FROM COMPSN ",YLATT*180.d0/pi,YLONT*180.d0/pi,
c    &               HTT
c        write (*,*) "FROM COMPSN ",ITREF,MIN,DTIME
       
** Compute the contribution due to earthquakes.
** It is assumed that the components of displacement,
** DNORTH,DWEST,DUP, do not vary from one reference
** frame to another given the accuracy of dislocation
** models.
      DO 10 I = 1, NUMEQ
          IF(ITEQK(I) .GT. ITREF) THEN 
               NTIME = 1
          ELSE
            NTIME = 0
          ENDIF
          IF(MIN .LT. ITEQK(I)) NTIME = NTIME - 1
          IF(NTIME .NE. 0) THEN 
             CALL RADII(EQLATR(I),RADMER,RADPAR)
             DDLAT = (YLAT - EQLATR(I))*RADMER
             DDLON = (YLON - EQLONR(I))*RADPAR
             DIST = DSQRT(DDLAT*DDLAT + DDLON*DDLON)
             IF(DIST .LE. EQRAD(I)) THEN
               ISTART = NLOC(I)
               IEND = NLOC(I) + NFP(I) - 1
               DO 5 JREC = ISTART,IEND
                 CALL DISLOC(YLAT,YLON,STRIKE(JREC),HL(JREC),
     &              EQLAT(JREC),EQLON(JREC),SSLIP(JREC),     
     &              DSLIP(JREC),DIP(JREC),DEPTH(JREC),
     &              WIDTH(JREC),DNORTH,DWEST,DUP)     
                 YLATT = YLATT + NTIME*DNORTH
                 YLONT = YLONT + NTIME*DWEST     
                 HTT   = HTT   + NTIME*DUP
    5          CONTINUE
             ENDIF
          ENDIF  
   10 CONTINUE

*** Compute contribution due to postseismic deformation
      CALL PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
c     write (*,*) "FROM COMPSN DISP ",DNORTH, DEAST, DUP
      CALL RADII (YLAT, RMER, RPAR)
      YLATT = YLATT + DNORTH/RMER
      YLONT = YLONT - DEAST/RPAR
      HTT  = HTT + DUP
      RETURN
      END
***************************************************************
      SUBROUTINE NEWCOR(YLAT,YLON,HTOLD,MIN1,MIN2,
     1         YLAT3,YLON3,HTNEW,DN,DE,DU,VN, VE, VU)

*** Predict coordinates at time MIN2 given coordinates at time MIN1.
*** Predict displacements from time MIN1 to time MIN2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

      HT = HTOLD

      CALL COMPSN(YLAT1,YLON1,HT1,YLAT,YLON,HT,
     1            MIN1,VN, VE, VU)

      CALL COMPSN(YLAT2,YLON2,HT2,YLAT,YLON,HT,
     1            MIN2, VN, VE, VU)

      YLAT3 = YLAT + YLAT2 - YLAT1
      YLON3 = YLON + YLON2 - YLON1
      HTNEW   = HT   + HT2   - HT1
c     write (*,*) "FROM NEWCOR ",YLAT3*180.d0/pi,YLON3*180.d0/pi,
c    &  HTNEW

      CALL RADII(YLAT,RADMER,RADPAR)

      DN =  RADMER * (YLAT2 - YLAT1)
      DE = -RADPAR * (YLON2 - YLON1)
      DU =  HT2 - HT1
c     write (*,*) "FROM NEWCOR displacement ",DN,DE,DU
c     write (*,*) "FROM NEWCOR ",YLAT*180.d0/pi,YLON*180.d0/pi,HT
c     write (*,*) "FROM NEWCOR ",YLAT1*180.d0/pi,YLON1*180.d0/pi,HT1
c     write (*,*) "FROM NEWCOR ",YLAT2*180.d0/pi,YLON2*180.d0/pi,HT2
c     write (*,*) "FROM NEWCOR ",YLAT3*180.d0/pi,YLON3*180.d0/pi,HTNEW

      RETURN
      END
******************************************************************
      subroutine PREDV (ylat,ylon,eht,date,iopt,jregn,vn,ve,vu) 

** Predict velocity in iopt reference frame       

** ylat       input - north latitude (radians)
** ylon       input - west longitude (radians)
** eht        input - ellipsoid height (meters)
** date       input - date (decimal years)
** iopt       input - reference frame
** jregn      output - deformation region
** vn         output - northward velocity in mm/yr
** ve         output - eastward velocity in mm/yr
** vu         output - upward velocity in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      logical  Is_iopt_NAD83

** Get reference latitude (RLAT) and reference longitude (RLON)

C  The following 2 lines were added on 07/22/2015 after Rich found this bug
         elon = -ylon
         call TOXYZ(ylat, elon, eht, x, y, z)

c        IF(IOPT .EQ. 0 .OR. IOPT .EQ. 1) THEN   !Velocity grids are No longer in NAD83
         IF(IOPT .EQ. 15) THEN                   !They are in ITRF2008
            RLAT = YLAT
            RLON = YLON
         ELSE
c            elon = -ylon                         !Was a bug, commented out on 07222015
c            call TOXYZ(ylat, elon, eht, x, y, z) !Was a bug, commented out on 07222015
c            CALL XTONAD(X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !No longer NAD83
             CALL XTO08 (X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !Positions should be in ITRF2008
         ENDIF

** Get deformation region

         CALL GETREG(RLAT,RLON,JREGN)
c        write (*,*) JREGN
         IF (JREGN .EQ. 0) THEN
           VN = 0.D0
           VE = 0.D0
           VU = 0.D0
           RETURN
         ENDIF
         CALL COMVEL( RLAT, RLON, JREGN, VN, VE, VU)       !Those velocities are in ITRF2008

** Convert  velocity to reference of iopt, if iopt != NAD83   !No, was in ITRF2008, not in NAD83  (since 09/12/2014)

         Is_iopt_NAD83 = (iopt == 1)
c        IF (IOPT .NE. 0 .AND. IOPT .NE. 1) THEN
         IF (IOPT .NE. 15) THEN
           CALL TOVXYZ( YLAT, ELON, VN, VE, VU, VX, VY, VZ)
           if (Is_iopt_NAD83) then
c            CALL VTRANF( X, Y, Z, VX, VY, VZ, 1, IOPT)
             CALL VTRANF( X, Y, Z, VX, VY, VZ, 15, IOPT)
           else
             CALL VTRANF_IERS( X, Y, Z, VX, VY, VZ, 15, IOPT)
           endif
           CALL TOVNEU( YLAT, ELON, VX, VY, VZ, VN, VE, VU)
         ENDIF
c        write (*,*) "From PREDV ",VN, VE, VU

         RETURN
         END

****************************************************************
      subroutine TRFVEL

*** Transform velocities from one reference frame to another

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 16)
      character    nameif*80,name24*80
      character    namef*30
      character    frame1*24, frame2*24
      character    option*1
      character    vopt*1, LATDIR*1,LONDIR*1
      character    record*120                 
      LOGICAL      Is_inp_NAD83,Is_out_NAD83
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

      write( luout, 100)
  100 format(
     1  ' Please enter the name of the file to contain '/
     1  ' the transformed velocities. ')
      read( luin, '(a30)',err=600,iostat=ios) namef
      if (ios /= 0) goto 600
      
      open( i2, file = namef, status = 'unknown')
      CALL HEADER

  105 write( luout, 110)
  110 format( /'*******************************'/
     1   ' Enter the reference frame of the input velocities.')
      call MENU1(iopt1, frame1)
      if (iopt1 .lt. 1 .or. iopt1 .gt. numref) then
      write( luout, *) ' Improper selection -- try again.'
      go to 105
      endif
      Is_inp_NAD83 = (iopt1 ==  1)

  115 write( luout, 120)
  120 format( /' Enter the reference frame for the output velocities.')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
      write( luout, *) 'Improper selection -- try again.'
      go to 115
      endif
      Is_out_NAD83 = (iopt2 ==  1)

      write( i2, 125) frame1, frame2
  125 format( ' TRANSFORMING VELOCITIES FROM ', A24, ' TO ', a24//
     1   16X, ' INPUT VELOCITIES      OUTPUT VELOCITIES'/)

  130 write( luout, 140)
  140 format( /'**********************************'/
     1  ' Velocities will be transformed at each specified point.'/
     1  ' Please indicate how you wish to input points.'/
     1  '    0...No more points.  Return to main menu.'/
     1  '    1...Individual points entered interactively.'/
     1  '    2...Transform velocities contained in batch file '/
     1  '        of delimited records of the form: '/
     1  '        LAT,LON,VN,VE,VU,TEXT ' /
     1  '        LAT = latitude in degrees (positive north/DBL PREC)'/
     1  '        LON = longitude in degrees (positive west/DBL PREC)'/
     1  '        VN = northward velocity in mm/yr (DBL PREC) '/
     1  '        VE = eastwars velocity in mm/yr (DBL PREC) '/
     1  '        VU = upward velocity in mm/yr (DBL PREC) '/
     1  '        TEXT = descriptive text (CHARACTER*24) '/
     1  '        Example: '/
     1  '        40.731671553,112.212671753, 3.7,3.8,-2.4,SALT AIR '/)
      read( luin,'(a1)',err=601,iostat=ios) option
      if (ios /= 0) goto 601
      if (option .eq. '0') then
        go to 500
      elseif (option .eq. '1') then
        call GETPNT( latd, latm, slat, LATDIR, lond, lonm, slon,
     1       LONDIR, name24, x, y, z, ylat, ylon, eht)
        elon = - ylon
        call GETVLY( ylat, elon, vx, vy, vz, vn, ve, vu, vopt, 211)
        vx1 = vx
        vy1 = vy
        vz1 = vz
        if (Is_inp_NAD83 .or. Is_out_NAD83) then
          call VTRANF( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
        else
          call VTRANF_IERS( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
        endif
        call TOVNEU( ylat, elon, vx1, vy1, vz1, vn1, ve1, vu1)
        xlat = (ylat*rhosec)/3600.d0
        xlon = (ylon*rhosec)/3600.d0
        call PRNTVL(vn, ve, vu, vx, vy, vz, vn1, ve1, vu1,
     1               vx1, vy1, vz1, name24, 1,xlat, xlon)

        go to 130
      elseif (option .eq. '2') then
         eht = 0.d0
         write (luout, 200)
  200    format(/' Enter name of input file: ')
         read(luin, '(a)',err=602,iostat=ios) nameif
         if (ios /= 0) goto 602

         open (i1, file = nameif, status = 'old')
  210    read(i1,'(a)',end = 220,err=603,iostat=ios) record          
c 210    read(i1, *, end = 220,err=603,iostat=ios) xlat, xlon,vn,ve,vu, name24
         call interprate_velocity_record (record,xlat,xlon,vn,ve,
     &                                    vu,name24)
         if (ios /= 0) goto 603
         ylat = (xlat*3600.d0) / rhosec
         ylon = (xlon*3600.d0) / rhosec
         elon = -ylon
         call TOXYZ (ylat, elon, eht, x, y, z)
         call TOVXYZ (ylat, elon, vn,ve,vu,vx,vy,vz)
         vx1 = vx
         vy1 = vy
         vz1 = vz
         if (Is_inp_NAD83 .or. Is_out_NAD83) then
           call VTRANF ( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
         else
           call VTRANF_IERS ( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
         endif
         call TOVNEU (ylat, elon, vx1, vy1, vz1, vn1, ve1, vu1)
         call PRNTVL (vn, ve,vu,vx,vy,vz,vn1,ve1,vu1,vx1,vy1,vz1,
     1      name24, 0, xlat, xlon)
         go to 210
  220    close (i1, status = 'keep')
      endif
  500 continue
      close (i2, status = 'keep')
      return

  600 write (*,'(/)')
      write (*,*) "Failed to read file name in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)')
      write (*,*) "Failed to read option in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)')
      write (*,*) "Failed to read input file name in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)')
      write(*,*)"Failed to read input file format in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      end
*******************************************************************
      SUBROUTINE VTRANF(X,Y,Z,VX,VY,VZ, IOPT1, IOPT2)

*** Convert velocity from reference frame of IOPT1 to 
*** reference frame of IOPT2.
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 16)
      common /tranpa/ tx(numref), ty(numref), tz(numref),
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref),
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      IF(IOPT1 .le. numref .and. IOPT2. le. numref
     &   .and. IOPT1 .gt. 0 .and. IOPT2 .gt. 0 ) THEN

*** Convert from mm/yr to m/yr
         VX = VX /1000.d0
         VY = VY / 1000.d0
         VZ = VZ / 1000.d0

*** From IOPT1 to ITRF94 
*** (following equations use approximations assuming
*** that rotations and scale change are small)
         WX = -drx(iopt1)           
         WY = -dry(iopt1)      
         WZ = -drz(iopt1)      
         DS = -dscale(iopt1)
         VX = VX - dtx(iopt1) + DS*X + WZ*Y - WY*Z
         VY = VY - dty(iopt1) - WZ*X  +DS*Y + WX*Z
         VZ = VZ - dtz(iopt1) + WY*X - WX*Y + DS*Z

*** From ITRF94 to IOPT2 reference frame
*** (following equations use approximations assuming
***  that rotations and scale change are small)
         WX = drx(iopt2)
         WY = dry(iopt2)
         WZ = drz(iopt2)
         DS = dscale(iopt2)
         VX = VX + dtx(iopt2) + DS*X + WZ*Y - WY*Z
         VY = VY + dty(iopt2) - WZ*X + DS*Y + WX*Z 
         VZ = VZ + dtz(iopt2) + WY*X - WX*Y + DS*Z

*** FROM m/yr to mm/yr
         VX = VX * 1000.d0
         VY = VY * 1000.d0
         VZ = VZ * 1000.d0

      ELSE
         write(luout,*) ' Improper reference frame in routine vtranf'
         stop
      ENDIF
c     write (*,*) "Inside VTRANF  ",Vx,Vy,Vz

      RETURN
      END
******************************************************
      SUBROUTINE VTRANF_IERS(X,Y,Z,VX,VY,VZ, IOPT1, IOPT2)

*** Convert velocity from reference frame of IOPT1 to 
*** reference frame of IOPT2.
****************
C  Important note:
C  The parameters in common block tranpa are computed using the IGS values of ITRF96==>ITRF97
C  The parameters in common block tranpa1 are computed using the IERS values of ITRF96==>ITRF97

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 16)
      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      IF(IOPT1 .le. numref .and. IOPT2. le. numref
     &   .and. IOPT1 .gt. 0 .and. IOPT2 .gt. 0 ) THEN

*** Convert from mm/yr to m/yr
         VX = VX /1000.d0
         VY = VY / 1000.d0
         VZ = VZ / 1000.d0

*** From IOPT1 to ITRF94 
*** (following equations use approximations assuming
*** that rotations and scale change are small)
         WX = -drx1(iopt1)           
         WY = -dry1(iopt1)      
         WZ = -drz1(iopt1)      
         DS = -dscale1(iopt1)
         VX = VX - dtx1(iopt1) + DS*X + WZ*Y - WY*Z
         VY = VY - dty1(iopt1) - WZ*X  +DS*Y + WX*Z
         VZ = VZ - dtz1(iopt1) + WY*X - WX*Y + DS*Z

*** From ITRF94 to IOPT2 reference frame
*** (following equations use approximations assuming
***  that rotations and scale change are small)
         WX = drx1(iopt2)
         WY = dry1(iopt2)
         WZ = drz1(iopt2)
         DS = dscale1(iopt2)
         VX = VX + dtx1(iopt2) + DS*X + WZ*Y - WY*Z
         VY = VY + dty1(iopt2) - WZ*X + DS*Y + WX*Z 
         VZ = VZ + dtz1(iopt2) + WY*X - WX*Y + DS*Z

*** FROM m/yr to mm/yr
         VX = VX * 1000.d0
         VY = VY * 1000.d0
         VZ = VZ * 1000.d0

      ELSE
         write(luout,*) ' Improper reference frame in routine vtranf'
         stop
      ENDIF

      RETURN
      END
******************************************************
      subroutine PRNTVL(VN, VE, VU, VX, VY, VZ, VN1, VE1, VU1,
     1                  VX1, VY1, VZ1, NAME24, IPRINT,XLAT,XLON)

*** Print transformed velocities

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      character    name24*80
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      if (iprint .eq. 1) then
      write( luout, 100) vn1, ve1, vu1, vx1, vy1, vz1
  100 format( ' ****************************************'/
     1  ' New northward velocity = ', f8.2, ' mm/yr' /
     1  ' New eastward velocity  = ', f8.2, ' mm/yr'/
     1  ' New upward velocity    = ', f8.2, ' mm/yr'/
     1  ' New x velocity         = ', f8.2, ' mm/yr'/
     1  ' New y velocity         = ', f8.2, ' mm/yr'/
     1  ' New z velocity         = ', f8.2, ' mm/yr'/)
      endif 

      write( i2, 200) name24, xlat, xlon,
     1    vn, vn1, ve, ve1, vu, vu1,
     1    vx, vx1, vy, vy1, vz, vz1
  200 format(1x, a24 /
     1   1x, 'latitude = ',F14.9,2x,'longitude = ',F14.9 /
     1   5x, 'northward velocity ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'eastward velocity  ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'upward velocity    ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'x velocity         ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'y velocity         ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'z velocity         ', f8.2, 6x, f8.2, ' mm/yr' /)

      return 
      end
*******************************************************************
       SUBROUTINE HEADER

       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       IMPLICIT INTEGER*4 (I-N)
       character  HTDP_version*10
       COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
       COMMON /VERSION/ HTDP_version                        

       WRITE(I2, 10) HTDP_version
   10  FORMAT(' HTDP (VERSION ',a,') OUTPUT' / )
       RETURN
       END
*********************************************
      SUBROUTINE GETMDY(MONTH, IDAY, IYEAR, DATE, MINS, TEST)

*** Read month-day-year and convert to decimal years
*** and Julian time in minutes      
***    MONTH      output - number from 1 to 12
***    IDAY       output - number from 1 to 31
***    IYEAR      output - must be after 1906
***    DATE       output - corresponding time in decimal years
***    MINS       output - corresponding julian time in minutes
***    TEST       output - if (true) then there is an error

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION M(12)
      LOGICAL TEST
      CHARACTER   TOPT*1
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      M(1) = 31
      M(2) = 28
      M(3) = 31
      M(4) = 30
      M(5) = 31
      M(6) = 30
      M(7) = 31
      M(8) = 31
      M(9) = 30
      M(10) = 31
      M(11) = 30
      M(12) = 31

   3  WRITE (LUOUT,1)
   1  FORMAT (' How do you wish to enter the time?'/
     1        '     1. month-day-year in free format'/
     1        '        for example, 5,12,1979'/
     1        '        represents May 12, 1979'/
     1        '     2. decimal year'/
     1        '        for example the entry 1979.359'/
     1        '        represents UTC midnight at the begining'/
     1        '        of May 12, 1979.'/) 

      READ (LUIN,5,err=100,iostat=ios) TOPT
      if (ios /= 0) goto 100
   5  format (A1)

      IF (TOPT .EQ. '1') Then
          write (luout,*) ' Enter month-day-year '
          READ(LUIN,*,err=101,iostat=ios) MONTH,IDAY,IYEAR
          if (ios /= 0)  goto 101
    
          IF(IYEAR .le. 1906) THEN
               WRITE(LUOUT,10)
   10          FORMAT(' The model is not valid for dates prior ',
     1           'to 1906.'/)
               TEST = .TRUE.
               RETURN
          ENDIF

          IF(MONTH .le. 0 .or. MONTH .gt. 12) THEN
               WRITE(LUOUT,20)
   20          FORMAT(' Improper month specified.'/)
               TEST = .TRUE.
               RETURN
          ENDIF

          IF(IDAY .le. 0 .or. IDAY .gt. 31) THEN
               WRITE(LUOUT,30)
   30          FORMAT(' Improper day specified.'/)
               TEST = .TRUE.
               RETURN
          ENDIF

          CALL IYMDMJ(IYEAR, MONTH, IDAY, MJD)
          CALL IYMDMJ(IYEAR, 1, 1, MJD0)
          IYEAR1 = IYEAR + 1
          CALL IYMDMJ(IYEAR1, 1, 1, MJD1)
          DAY = DBLE(MJD - MJD0)
          DENOM = DBLE(MJD1 - MJD0)
          DATE = DBLE(IYEAR) + (DAY / DENOM)
          MINS = MJD * 24 * 60
          TEST = .FALSE.
          RETURN

      ELSEIF (TOPT .EQ. '2') then
          write(luout,*) ' Enter decimal year '
          READ (LUIN, *,err=102,iostat=ios) DATE
          if (ios /= 0)  goto 102
          
          IF (DATE .lt. 1906.0d0) then
             write (luout, 10)
             TEST = .TRUE.
             RETURN
          ENDIF
**** add small increment to circumvent round-off error
c         DATE = DATE + 0.0004D0
****
          IYEAR = DATE
          CALL IYMDMJ(IYEAR, 1, 1, MJD0)
          IYEAR1 = IYEAR + 1
          CALL IYMDMJ(IYEAR1, 1, 1, MJD1)
          LEAP = 0
          IF ((MJD1 - MJD0) .eq. 366) LEAP = 1
          REMDAY = (DATE - IYEAR)* (MJD1 - MJD0)
          IBEGIN = 0
          ITOTAL = 31
          IF (REMDAY .LT. ITOTAL) then
              MONTH = 1
              IDAY = REMDAY - IBEGIN + 1
              CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
              MINS = MJD * 24 * 60
              TEST = .FALSE.
              RETURN
          ENDIF
          IBEGIN = ITOTAL
          ITOTAL = ITOTAL + LEAP
          DO I = 2, 12
             ITOTAL = ITOTAL + M(I)
             IF (REMDAY .LT. ITOTAL) then
                  MONTH = I
                  IDAY = REMDAY - IBEGIN + 1
                  TEST = .FALSE.
                  CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
                  MINS = MJD * 24 * 60
                  RETURN
             ENDIF
             IBEGIN = ITOTAL
           ENDDO
           Write(LUOUT, 60)
   60      Format (' Error could not convert Decimal years to '
     1             ,'month-day-year')
           TEST = .TRUE.
           RETURN
      ELSE
         write (luout, 70)
   70    format (' Improper entry')
         Go TO 3
      ENDIF

 100   write (*,'(/)')
       write (*,*) "Failed to read option in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

 101   write (*,'(/)')
       write (*,*) "Wrong MDY input in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

 102   write (*,'(/)')
       write (*,*) "Wrong decimal year in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

       END
************************************
      SUBROUTINE IYMDMJ( IYR, IMON, IDAY, MJD )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       IYMDMJ
C VERSION:    Sep. 17, 2010
C WRITTEN BY: R. SNAY (after M. SCHENEWERK)
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE 
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IDAY              DAY
C IMON              MONTH
C IYR               YEAR
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C MJD               MODIFIED JULIAN DATE 
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C A                 TEMPORARY STORAGE
C B                 TEMPORARY STORAGE
C C                 TEMPORARY STORAGE
C D                 TEMPORARY STORAGE
C IMOP              TEMPORARY STORAGE
C IYRP              TEMPORARY STORAGE
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:     DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:       
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
C                              I.E. 1992 RATHER THAN 92.  
C
C********1*********2*********3*********4*********5*********6*********7**
C::LAST MODIFICATION
C::8909.06, MSS, DOC STANDARD IMPLIMENTED
C::9004.17, MSS, CHANGE ORDER YY MM DD
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      INTEGER*4     A, B, C, D

      IYRP = IYR
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IF( IMON .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMON + 12
      ELSE
        IMOP= IMON
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP*0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0*IYRP
      D=  30.6001D0*(IMOP + 1)
      MJD =  (B + C + D + IDAY - 679006) 
C      
      RETURN
      END
*****************************************************
      SUBROUTINE PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
********
*   Compute total postseismic displacement for all earthquakes
*
* INPUT
*   YLAT       latitude of point in radians, positive north
*   YLON       longitude of point in radians, positive west
*   MIN        modified julian date of reference epoch for new coordinates
*               in minutes
*
*   DNORTH     Total northward postseismic displacement at point during
*              period from ITREF to MIN in meters
*   DEAST      TOTAL eastward postseismic displacement
*   DUP        Total upward postseismic displacement 
*******

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL INSIDE

      parameter (NUMPSG = 1)
      COMMON /CONST/ A, F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /PGRID/ PS(18000)
      DIMENSION ITEQ(NUMPSG)      
      DIMENSION TAU(NUMPSG)   
      DIMENSION WEI(2,2)
      DIMENSION AMP(2,2,3)

*** Relaxation constant (in years) for 2002 Denali earthquake
      TAU(1) = 5.0D0

*** Modofied Julian Date (in minutes) for the 2002 Denali earthquake
      IYEAR = 2002
      IMO = 11
      IDAY = 3
      CALL IYMDMJ(IYEAR,IMO,IDAY, MJD)
      ITEQ(1) = MJD*60*24

      DNORTH = 0.0D0
      DEAST = 0.0D0
      DUP = 0.0D0

      DO K = 1, NUMPSG
*** Check if the point is inside the grid
         POSX = YLON*180.d0/PI
         POSX = 360.d0 - POSX
         IF (POSX .GT. 360.D0) POSX = POSX - 360.D0
         POSY = YLAT*180.D0/PI
         CALL GRDCHK(POSX, POSY, PSGLX(K), PSGUX(K),
     1            PSGLY(K), PSGUY(K), INSIDE)

         IF (INSIDE ) THEN
*** Get the indices for the lower left-hand corner of the grid
         CALL PSGWEI(POSX,POSY,K,I,J,WEI)
*** Get the displacement amplitude at the four corners
         CALL GRDAMP(K,I,J,AMP,PS)

         ANORTH = WEI(1,1)*AMP(1,1,1) + WEI(1,2)*AMP(1,2,1)
     1          + WEI(2,1)*AMP(2,1,1) + WEI(2,2)*AMP(2,2,1)     
         AEAST  = WEI(1,1)*AMP(1,1,2) + WEI(1,2)*AMP(1,2,2)
     1          + WEI(2,1)*AMP(2,1,2) + WEI(2,2)*AMP(2,2,2)
         AUP    = WEI(1,1)*AMP(1,1,3) + WEI(1,2)*AMP(1,2,3)
     1          + WEI(2,1)*AMP(2,1,3) + WEI(2,2)*AMP(2,2,3)

c        write (6, 30) anorth, aeast, aup
c  30    format (1x, 3f15.5)

*** Convert amplitudes from mm to meters
         ANORTH = ANORTH / 1000.D0
         AEAST =  AEAST  / 1000.D0
         AUP =    AUP    / 1000.D0

         IF (MIN .GT. ITEQ(K)) THEN
            DTIME = DBLE(MIN - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH + ANORTH*FACTOR
            DEAST = DEAST + AEAST*FACTOR
            DUP = DUP + AUP*FACTOR
         ENDIF
         IF (ITREF .GT. ITEQ(K)) THEN
            DTIME = DBLE(ITREF - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH - ANORTH*FACTOR
            DEAST = DEAST - AEAST*FACTOR
            DUP = DUP - AUP*FACTOR
         ENDIF
         ENDIF
      ENDDO
      RETURN
      END

****************************************************
      SUBROUTINE GRDCHK (POSX, POSY, GRDLX, GRDUX,
     1             GRDLY, GRDUY, INSIDE)

C
C ROUTINE CHECKS IF THE POINT HAVING COORDINATES (POSX, POSY)
C IS WITHIN THE REGION SPANNED BY THE GRID
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL INSIDE

      INSIDE = .TRUE.

      IF (POSX .LT. GRDLX .OR. POSX .GT. GRDUX) THEN
         INSIDE = .FALSE.
      ENDIF
      IF (POSY .LT. GRDLY .OR. POSY .GT. GRDUY) THEN
         INSIDE = .FALSE.
      ENDIF

      RETURN
      END
*******************************************************************
      SUBROUTINE PSGWEI (POSX, POSY, K, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C POSX         LONGITUDE OF POINT IN DEGREES, POSITIVE EAST
C POSY         LATITUDE OF POINT IN DEGREES, POSITIVE NORTH
C K            ID OF EARTHQUAKE GRID                   
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
C              FOR THE CORNER VECTORS
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   PSDISP
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /PSGRID/, /CONST/
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9511.09, RAS, MODIFIED FOR HTDP
C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
C********1*********2*********3*********4*********5*********6*********7**
    
C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMPSG = 1)
      DIMENSION WEI(2,2)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG), 
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (PSGUX(K) - PSGLX(K)) / ICNTPX(K)
      STEPY = (PSGUY(K) - PSGLY(K)) / ICNTPY(K)
      I = IDINT((POSX - PSGLX(K))/STEPX) + 1
      J = IDINT((POSY - PSGLY(K))/STEPY) + 1
c     write(6,1001) K, I, J
c1001 format(1x, 'quake = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = PSGLX(K) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = PSGLY(K) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END

C*********************************************************************
      SUBROUTINE GRDAMP (K, I, J, AMP, PS)
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE AMPLITUDES OF THE FOUR
C              GRID NODES OFGRID K WHERE I,J ARE THE INDICES OF
C              THE LOWER LEFT HAND CORNER
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C
C K            ID OF EARTHQUAKE CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C PS           THE ARRAY CONTAINING ALL THE GRIDDED AMPLITUDES
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C AMP          A TWO BY TWO ARRAY CONTAINING THE 3D AMPLITUDES
C              FOR THE CORNERS OF THE GRID
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   PSDISP
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::2011.08.17, RAS, ORIGINAL CREATION FOR TRANS4D
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION AMP(2,2,3), PS(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IPSGRD(K, I + II, J + IJ, IVEC)
               AMP(II + 1, IJ + 1, IVEC) = PS(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

      RETURN
      END

C***************************************************
      INTEGER FUNCTION IPSGRD(IGRID, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMPSG = 1)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)

      IPSGRD = NBASEP(IGRID) +
     1      3 * ((J - 1) * (ICNTPX(IGRID) + 1) +  (I - 1)) + IVEC

      RETURN
      END
C---------------------------------------------------------------------
      subroutine extract_name (name,i)

C  In f90, use trim(name) to truncate all spaces after the name.
C  But "trim" were not available in f77. This is what this subroutine is for

      implicit none
      integer*4 i
      character name*80,scratch*80

      do i=80,1,-1
        if (name(i:i) /= ' ') exit 
      enddo
      scratch = name(1:i)
      name    = scratch

      return
      end
C-----------------------------------------------------------------------------------
         subroutine interprate_XYZ_record (record,x,y,z,name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,z
         character   name*80,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,zzzz*80

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get z or h     

         zzzz = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             zzzz(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (zzzz,*) z
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end
C-----------------------------------------------------------------------------------
         subroutine interprate_XYZVxVyVz_record (record,x,y,z,Vx,Vy,Vz,
     &                                           name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,z,Vx,Vy,Vz
         character   name*80,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,zzzz*80

         record1 = trim(adjustl(record))
         length  = len_trim(record1)
         chars   = ' '

         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get z or h     

         zzzz = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             zzzz(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (zzzz,*) z
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get Vx 

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) Vx
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get Vy 

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) Vy
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get z or h     

         zzzz = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             zzzz(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (zzzz,*) Vz
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
C  Get name

         name = trim(record1)

C  Done

         return 
         end
C-----------------------------------------------------------------------------------
         subroutine interprate_latlon_record (record,x,y,name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,z
         character   name*24,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,zzzz*80

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end
C-----------------------------------------------------------------------------------
         subroutine interprate_velocity_record (record,x,y,vn,ve,
     &                                          vu,name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,vn,ve,vu
         character   name*80,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,vnnn*80,veee*80,vuuu*80

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))        !adjust left then trim all trailing spaces
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get vn         

         vnnn = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             vnnn(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (vnnn,*) vn
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get ve         

         veee = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             veee(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (veee,*) ve
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get vu         

         vuuu = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             vuuu(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (vuuu,*) vu
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end
C-----------------------------------------------------------------------------------
        logical function am_I_in_or_near_CONUS (fi,la)

        implicit none

        integer*4         nrows,ncols,xcell,ycell,rec_num
        integer*2         sea,code,code_UL,code_UR,code_LL,code_LR
        real*8            fi,la,fi_max,la_min,dlamda,dfi,fi_min
        real*8            fi1,la1,fi2,la2,la_max

        character         grid*80                              
c       logical           am_I_in_or_near_CONUS,CONUS_on_land
        logical           CONUS_on_land

C  Constants

        parameter (fi_max = 50.d0,
     &             fi_min = 23.5d0,
     &             la_min = 235.d0, 
     &             la_max = 295.d0, 
     &             sea    = 0)

C  Initialize this function

        am_I_in_or_near_CONUS = .false.

C  If clearly outside Conus

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_CONUS = .true.   
          return
        endif

        return
        end
C-----------------------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_AK(fi,la)

C  Given your latitude and longitude, this little routine tells you if you  
C  are in Alaska (AK) or within 10 km of its boundaries.                 

C  "AK_direct_cell_by_cell_1min" is an unformatted direct access file of nowrs*ncols records
C  each record is of length 2 bytes (integer*2), which contains a geographic code defined as follows:

C  1) The AK code is 278 
C  2) The codes for Canada are 265 to 276
C  4) The code for sea is 0

C  Geographic conditions to be in or near AK:          

C  1) The point is on land in AK    

C            a) code is 278

C  2) In the Canada near the AK:                      

C            a) code = 265 to 276
C            c) Point is no more than 10km from AK                          

C  3) In the sea near AK:                        

C            a) code = 0
C            c) Point is no more than 10km from AK                          

        implicit none

        integer*4         nrows,ncols,xcell,ycell,rec_num
        integer*2         sea,code,code_UL,code_UR,code_LL,code_LR
        real*8            fi,la,fi_max,dlamda,dfi,la_max
        real*8            fi1,la1,fi2,la2,la_min,fi_min

        character         grid*27,path*46                
c       logical           am_I_in_or_near_AK,AK_on_land
        logical           AK_on_land

C  Constants

        parameter (nrows     = 1291, 
     &             ncols     = 3511,     
     &             dfi       = 1.d0/60.d0,
     &             dlamda    = 1.d0/60.d0,
     &             fi_max    = 71.5d0,
     &             fi_min    = 50.d0, 
     &             la_min    = 172.d0, 
     &             la_max    = 230.5d0,
     &             sea       = 0)

C  Initialize this function

        am_I_in_or_near_AK = .false.

C  When it is clearly not AK

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_AK = .true.    
          return
        endif

C  On 1/22/2016 I disabled the following code following Giovanni's suggestion    
C  and for simplicity and portability

        return
        end
C----------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_HI(fi,la)

C  Given your latitude and longitude, this little routine tells you if you  
C  are in the Hawaian Islands (Hawaii, Kahoolawe,Kauai,Lanai,Maui,Molokai, Nihau and Oahu)
C  or within 5 km of their coastlines.                 

C  "HI_Hawaii.gmt", "HI_Kahoolawe.gmt", "HI_Kauai.gmt", "HI_Lanai.gmt", "HI_Maui.gmt", "HI_Molokai.gmt",
C  "HI_Nihau.gmt" and "HI_Oahu.gmt" are ASCII files that contain coastline polygons for the islands.
C  Each polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.

        implicit none

        integer*4    i,NPC,n_Hawaii,n_Kahoolawe,n_Kauai,n_Lanai,n_Maui
        integer*4    n_Molokai,n_Nihau,n_Oahu
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_HI,Hawaii,Kahoolawe,Kauai,Lanai
        logical      Hawaii,Kahoolawe,Kauai,Lanai
        logical      Maui,Molokai,Nihau,Oahu
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_Hawaii    = 329,
     &             n_Kahoolawe = 143,
     &             n_Kauai     = 107,
     &             n_Lanai     =  52,
     &             n_Maui      = 205,
     &             n_Molokai   = 114,
     &             n_Nihau     =  68,
     &             n_Oahu      = 123,
     &             fi_max    = 24.d0,
     &             fi_min    = 18.d0, 
     &             la_min    = 199.d0, 
     &             la_max    = 207.d0)

        real*8            X1(n_Hawaii)    ,Y1(n_Hawaii) !I am not sure the cheap compilers in subversion can handle dynamic memory
        real*8            X2(n_Kahoolawe) ,Y2(n_Kahoolawe)
        real*8            X3(n_Kauai)     ,Y3(n_Kauai)
        real*8            X4(n_Lanai)     ,Y4(n_Lanai)
        real*8            X5(n_Maui)      ,Y5(n_Maui)
        real*8            X6(n_Molokai)   ,Y6(n_Molokai)
        real*8            X7(n_Nihau)     ,Y7(n_Nihau)
        real*8            X8(n_Oahu)      ,Y8(n_Oahu)

C  Initialize this function

        am_I_in_or_near_HI = .false.

C  When it is clearly not HI

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_HI = .true.    
          return
        endif

C  On 1/22/2016 I disabled the following code following Giovanni's suggestion
C  and to simplify.    

        return
        end
C-------------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_PR(fi,la)

C  Given your latitude and longitude, this little routine tells you if you  
C  are in Puerto Rico Islands (the main Island, Culebra, Desecheo, Mona or Viequez) 
C  or within 10 km of their coastlines.                 

C  "PR_PuertoRico.gmt", "PR_Culebra", "PR_Mona.gmt", "PR_Viequez.gmt", and PR_"Desecheo.gmt" 
C  are ASCII files that contain coastline polygons for each one of the 5 PR islands.
C  Each polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_PR,n_Culebra,n_Mona,n_Desecheo,n_Viequez
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_PR,PR
        logical      PR
        logical      Culebra,Mona,Desecheo,Viequez
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_PR       = 228,
     &             n_Culebra  = 185,
     &             n_Mona     =  52,
     &             n_Desecheo =  34,
     &             n_Viequez  = 309,
     &             fi_max     =  19.d0,
     &             fi_min     =  17.d0,
     &             la_max     = 295.d0,
     &             la_min     = 292.d0)

        real*8            X1(n_PR)      ,Y1(n_PR)   !I am not sure the cheap compilers in subversion can handle dynamic memory
        real*8            X2(n_Culebra) ,Y2(n_Culebra)
        real*8            X3(n_Mona)    ,Y3(n_Mona)
        real*8            X4(n_Desecheo),Y4(n_Desecheo)
        real*8            X5(n_Viequez) ,Y5(n_Viequez)

C  Initialize this function

        am_I_in_or_near_PR = .false.

C  When it is clearly not PR

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_PR = .true.    
          return
        endif

C  On 1/22/2016 I disabled the following code following Giovanni's suggestion
C  and for simplicity.

        return
        end
C-------------------------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_VQ(fi,la)

C  Given your latitude and longitude, this little routine tells you if you are in
C  the Virgin  Islands (St. Croix, St John and St Thomas) or within 10 km of their coastlines.                 

C  "VQ_StCroix.gmt", "VQ_StJohn.gmt" and  "VQ_StThomas.gmt" are ASCII files that
C  contain coastline polygons for each one of the 3 VQ islands.
C  Each polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_VQ,n_StCroix,n_StJohn,n_StThomas
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_VQ
        logical      StCroix,StJohn,StThomas       
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_StCroix  = 142,
     &             n_StJohn   = 356,
     &             n_StThomas = 304,
     &             fi_max     = 18.4d0,
     &             fi_min     = 15.5d0,
     &             la_max     = 295.6d0,
     &             la_min     = 295.d0)

        real*8            X1(n_StCroix) ,Y1(n_StCroix)  !I am not sure the cheap compilers in subversion can handle dynamic memory
        real*8            X2(n_StJohn ) ,Y2(n_StJohn )
        real*8            X3(n_StThomas),Y3(n_StThomas)

C  Initialize this function

        am_I_in_or_near_VQ = .false.

C  When it is clearly not VQ

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_VQ = .true.    
          return
        endif

C  On 1/22/2016 the following code was disabled following Giovanni's suggestion
C  and for simplicity

        return
        end
C-------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_CQ (fi,la)

C  Given your latitude and longitude, this little routine tells you if you are in
C  the North Mariana Islands (Saipan, Tinian and Rota) or within 10 km of their coastlines.                 

C  "CQ_Rota.gmt", "CQ_Saipan.gmt" and "CQ_Tinian.gmt" are ASCII files that
C  contain coastline polygons for each one of the 3 CQ islands.
C  Each polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_CQ,n_Saipan ,n_Tinian,n_Rota     
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_CQ
        logical      Saipan,Tinian,Rota            
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_Saipan   = 122,
     &             n_Tinian   =  70,
     &             n_Rota     =  80,
     &             fi_max     =  15.6d0,
     &             fi_min     =  13.8d0,
     &             la_max     = 146.d0,
     &             la_min     = 144.8d0)

        real*8            X1(n_Saipan ),Y1(n_Saipan ) !I am not sure the cheap compilers in subversion can handle dynamic memory
        real*8            X2(n_Tinian ),Y2(n_Tinian )
        real*8            X3(n_Rota   ),Y3(n_Rota   )

C  Initialize this function

        am_I_in_or_near_CQ = .false.

C  When it is clearly not CQ

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_CQ = .true.    
          return
        endif

C  On 1/22/2016 the following code was disabled following Giovanni's suggestion
C  and for simplicity

        return
        end
C-------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_AS(fi,la)

C  Given your latitude and longitude, this little routine tells you if you are in
C  the American Samoa Islands (Ofu, Rose, Tau and Tutila) or within 10 km of their coastlines.                 

C  "AS_Ofu.gmt", "AS_Rose.gmt" and  " AS_Tau.gmt" and "AS_Tutila.gmt" are ASCII files that
C  contain coastline polygons for each one of the 4 AS islands.
C  Each polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_AS,n_Ofu,n_Rose,n_Tau,n_Tutila
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_AS
        logical      Ofu,Rose,Tau,Tutila           
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_Ofu      =  64,
     &             n_Rose     =   8,
     &             n_Tau      =  76,
     &             n_Tutila   = 239,
     &             fi_max     = -13.7d0,
     &             fi_min     = -14.7d0,
     &             la_max     = 190.8d0,
     &             la_min     = 189.d0)

        real*8            X1(n_Ofu )  ,Y1(n_Ofu )      !not sure the cheap compilers in subversion can handle dynamic memory
        real*8            X2(n_Rose)  ,Y2(n_Rose)
        real*8            X3(n_Tau )  ,Y3(n_Tau )
        real*8            X4(n_Tutila),Y4(n_Tutila)

C  Initialize this function

        am_I_in_or_near_AS = .false.

C  When it is clearly not AS

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_AS = .true.    
          return
        endif

C  On 1/22/2016 the following code was disabled following Giovanni's suggestion
C  and for simplicity

        return
        end
C-------------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_Guam(fi,la)

C  Given your latitude and longitude, this little routine tells you if you are in
C  the Guam Island or within 5 km of its coastlines.                 

C  "GQ_Guam.gmt", is an ASCII file that contains coastline polygons Guam.                             
C  This polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_GQ
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min
c       logical      am_I_in_or_near_Guam
        logical      GQ                            
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_GQ   = 300,
     &             fi_max = 13.7d0, 
     &             fi_min = 13.d0,
     &             la_max = 145.d0,
     &             la_min = 144.5d0)

        real*8            X1(n_GQ) ,Y1(n_GQ)              !I am not sure the cheap compilers in subversion can handle dynamic memory

C  Initialize this function

        am_I_in_or_near_Guam = .false.

C  When it is clearly not Guam

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_Guam = .true.    
          return
        endif

C  On 1/22/2016 the following code was disabled following Giovanni's suggestion
C  and for simplicity

        return
        end
C----------------------------------------------------------------------------------------------
        logical function am_I_in_or_near_KW(fi,la)

C  Given your latitude and longitude, this little routine tells you if you are in
C  the a little square around Kwajalein of the Marshal islands       

C  "Kwajalein.gmt", is an ASCII file that contains a rectangle around Kwajalein.                             
C  This polygon is closed, starts and ends at the same point. 
C  File format: (East) longitude and latitude in decimal degrees.


        implicit none

        integer*4    i,NPC,n_KW
        
        real*8       fi,la,fi1,la1,fi2,la2,fi_max,fi_min,la_max,la_min

c       logical      am_I_in_or_near_KW   
        logical      KW                            
        logical      UL,UR,LL,LR

C  Constants

        parameter (n_KW   = 5,
     &             fi_max = 10.5d0,
     &             fi_min =  5.d0 ,
     &             la_max = 174.d0,
     &             la_min = 165.7d0)

        real*8            X1(n_KW) ,Y1(n_KW)              !I am not sure the cheap compilers in subversion can handle dynamic memory

C  Initialize this function

        am_I_in_or_near_KW = .false.

C  When it is clearly not KW   

        if (fi>fi_max.or.fi<fi_min.or.la>la_max.or.la<la_min) then
          return
        else
          am_I_in_or_near_KW = .true.    
          return
        endif

C  On 1/22/2016 the following code was disabled following Giovanni's suggestion
C  and for simplicity

        return
        end
