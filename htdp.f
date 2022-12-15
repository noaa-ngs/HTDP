      PROGRAM HTDP

**************************************************************
*  NAME:       HTDP (Horizontal Time-Dependent Positioning)
*
*  WRITTEN BY: Michael Dennis, Jarir Saleh, Richard Snay, 
*              and Chris Pearson
*
*  PURPOSE:    Transform coordinates across time
*              and between reference frames
*
**************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      character    HTDP_version*8
      character    Version_date*20
      CHARACTER    OPTION*1

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /VERSION/ HTDP_version

C  You must change HTDP version and date here, if necessary

      HTDP_version = '3.5.0'
      Version_date = 'November 29, 2022'

*** Introduce variables for file IDs

      LUIN = 5
*       interactive input
      LUOUT = 6
*       interactive output
      I1 = 11
*       input of velocity grid in GETVEL      
*       input of earthquake parameters in GETEQ
*       input of Bluebook BFILE in DLACE, VELOC, UPDATE and TRFPOS
*       input of Bluebook GFILE in UPDATE
      I2 = 12
*       output of predicted displacements in DLACE
*       output of predicted velocities in VELOC
*       output of updated Bluebook BFILE in UPDATE
*       output of updated Bluebook GFILE in UPDATE
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
*       output of transformed Bluebook BFILE in TRFPOS

*** Obtain parameters defining crustal motion model
      CALL MODEL
       
*** Initialize transformation parameters between reference frames
      CALL SETTP 

*** Initialize conversion table between reference frame identifiers
      CALL SETRF

      WRITE(LUOUT,5) HTDP_version, Version_date
    5 FORMAT(
     1 '********************************************************'/
     1 '   HTDP (Horizontal Time-Dependent Positioning)         '/
     1 '   SOFTWARE VERSION: ', a8                               /
     1 '   VERSION DATE:     ', a20                              /)
      WRITE(LUOUT,501) 
  501 FORMAT(
     1 '   AUTHORS: Michael Dennis, Jarir Saleh, Richard Snay,  '/
     1 '            and Chris Pearson                          '//
     1 '   Web: https://geodesy.noaa.gov/TOOLS/Htdp/Htdp.shtml  '/
     1 '   Email: ngs.cors.htdp@noaa.gov                        '/
     1 '********************************************************')
      WRITE(LUOUT,10)
   10 FORMAT( 
     1 ' This software incorporates numerical models that',/
     3 ' characterize continuous crustal motion as well as ',/
     3 ' the episodic motion associated with earthquakes.'/)
      WRITE(LUOUT,11)
   11 FORMAT(
     5 ' The User Guide contains additional information and a set'/
     5 ' of exercises to familiarize users with the software.')

   25 WRITE(LUOUT,26)
   26 FORMAT('********************************************************'/
     1 ' MAIN MENU:',/
     6 '    0... Exit software.',/
     7 '    1... Estimate horizontal displacements between two dates.'/
     8 '    2... Estimate horizontal velocities.'/
     9 '    3... Transform observations to a specified reference '
     &           ,'frame and/or date.'/
     & '    4... Transform positions between reference frames '
     &           ,'and/or dates.'/
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

      A = 6378137.D0
      F = 1.D0 / 298.257222101D0
      E2 = F*(2.D0 - F)   !Calculate E2 from F rather than hard-coded as done previously
      AF = A / (1.D0 - F)
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
*** Region 1: San Andreas fault in central California         
*** Region 2: southern California
*** Region 3: northern California
*** Region 4: Pacific Northwest
*** Region 5: western CONUS
*** Region 6: CONUS
*** Region 7: St. Elias, Alaska
*** Region 8: south-central Alaska
*** Region 9: southeast Alaska
*** Region 10: all mainland Alaska

* Placeholder regions inserted to make compatible with new initbd.f
* that contains tectonic plates added in December 2020 (v3.3.0).
*** Region 11: placeholder region
*** Region 12: placeholder region
*** Region 13: placeholder region

*** Tectonic plates
*** All plate boundaries are from Bird 2003 with minor topological corrections.
*** Main four plates included through v3.2.9:
***   Region 14: North America (NA)
***   Region 15: Pacific (PA)
***   Region 16: Caribbean (CA)
***   Region 17: Mariana (MA)
***
*** Plates added or updated after v3.2.9
*** Primary plates:
***   Region 18: Africa (AF)
***   Region 19: Antarctica (AN)
***   Region 20: Australia (AU)
***   Region 21: Eurasia (EU)
***   Region 22: South America (SA)

*** Secondary plates:
***   Region 23: Amur (AM)
***   Region 24: Arabia (AR)
***   Region 25: Cocos (CO)
***   Region 26: India (IN)
***   Region 27: Juan de Fuca (JF)
***   Region 28: Nazca (NZ)
***   Region 29: Philippine Sea (PS)
***   Region 30: Scotia (SC)
***   Region 31: Somalia (SO)
***   Region 32: Sunda (SU)

*** Tertiary plates:
***   Region 33: Aegean Sea (AS)
***   Region 34: Altiplano (AP)
***   Region 35: Anatolia (AT)
***   Region 36: Balmoral Reef (BR)
***   Region 37: Banda Sea (BS)
***   Region 38: Birds Head (BH)
***   Region 39: Burma (BU)
***   Region 40: Caroline (CL)
***   Region 41: Conway Reef (CR)
***   Region 42: Easter (EA)
***   Region 43: Futuna (FT)
***   Region 44: Galapagos (GP)
***   Region 45: Juan Fernandez (JZ)
***   Region 46: Kermadec (KE)
***   Region 47: Manus (MN)
***   Region 48: Maoke (MO)
***   Region 49: Molucca Sea (MS)
***   Region 50: New Hebrides (NH)
***   Region 51: Niuafo'ou (NI)
***   Region 52: North Andes (ND)
***   Region 53: North Bismarck (NB)
***   Region 54: Okhotsk (OK)
***   Region 55: Okinawa (ON)
***   Region 56: Panama (PM)
***   Region 57: Rivera (RI)
***   Region 58: Sandwich (SW)
***   Region 59: Shetland (SL)
***   Region 60: Solomon Sea (SS)
***   Region 61: South Bismarck (SB)
***   Region 62: Timor (TI)
***   Region 63: Tonga (TO)
***   Region 64: Woodlark (WL)
***   Region 65: Yangtze (YA)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 65)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /BNDRY/ X(12000), Y(12000), NPOINT(70)
 
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
      parameter (NMREGN = 65)
      COMMON /BNDRY/ X(12000), Y(12000), NPOINT(70)
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
      SUBROUTINE POLYIN (X0,Y0,X,Y,N,NPC)
* This version is from TRANS4D v0.3.1 and replaces version in HTDP v3.2.9.
C     DETERMINES IF A POINT AT (X0,Y0) IS INSIDE OR OUTSIDE 
C     OF A CLOSED FIGURE DESCRIBED BY A SEQUENCE OF CONNECTED
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
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      DIMENSION X(N),Y(N)

      IS=0
      NPC=0
C
C     FIND STARTING POINT WHERE X(I).NE.X0
      IP=0
   10 IP=IP+1
      IF(X(IP)-X0) 15,12,16
   12 IF(IP.LE.N) GO TO 10
      WRITE(LUOUT,6001)
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
      IF(IL) 30,50,40
   30 IF(X(I)-X0) 90,32,34
   32 IS=-1
      GO TO 60
   34 IL=1
      GO TO 80
   40 IF(X(I)-X0) 42,44,90
   42 IL=-1
      GO TO 80
   44 IS=1
      GO TO 60
   50 IF(X(I)-X0) 52,55,54
   52 IL=-1
      IF(IS) 90,140,80
   54 IL=1
      IF(IS) 80,140,90
   55 IF(Y(I)-Y0) 57,120,58
   57 IF(YL-Y0) 90,120,120
   58 IF(YL-Y0) 120,120,90
C
   60 IL=0
      IF(Y(I)-Y0) 90,120,90
   80 IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)) 90,120,85
   85 NPC=NPC+1
   90 XL=X(I)
      YL=Y(I)
  100 CONTINUE
      NPC=MOD(NPC,2)
      RETURN
  120 NPC=2
      RETURN
 140  WRITE(LUOUT,6002)
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
C Compute the ITRF2008 velocity at a point in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 10)
      parameter (NMREGN = 65)
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
        
*** Subtract the number of placeholder regions
        IPLATE = IPLATE - 3
        
        ELON = - YLON
        HT = 0.D0
        CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
        CALL PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
        VX = VX * 1000.D0
        VY = VY * 1000.D0
        VZ = VZ * 1000.D0
        CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)

      ELSEIF(JREGN .GE. 1 .AND. JREGN .LE. NUMGRD) THEN

C*** Get indices for the lower left hand corner of the grid
C*** and get the weights for the four corners
        CALL GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C*** Get the velocity vectors at the four corners
        CALL GRDVEC (JREGN, I, J, VEL, B)

        VN = WEI(1,1) * VEL(1,1,1) + WEI(1,2) * VEL(1,2,1)
     *     + WEI(2,1) * VEL(2,1,1) + WEI(2,2) * VEL(2,2,1)

        VE = WEI(1,1) * VEL(1,1,2) + WEI(1,2) * VEL(1,2,2)
     *     + WEI(2,1) * VEL(2,1,2) + WEI(2,2) * VEL(2,2,2)
  
        VU = WEI(1,1) * VEL(1,1,3) + WEI(1,2) * VEL(1,2,3)
     *     + WEI(2,1) * VEL(2,1,3) + WEI(2,2) * VEL(2,2,3)        

C*** If the point in one of the four Alaskan regions,
C*** then set its vertical velocity to 0.0
        IF(JREGN .GE. 7 .AND. JREGN .LE. 10) THEN
           VU = 0.D0
        ENDIF

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
*** with coordinates X, Y, Z (in meters)
*** The resulting velocities--VX, VY, and VZ--will be in meters/yr
*** References: 
***   Altamimi et al. 2012. "ITRF2008 plate motion model" J. Geophys. Res.
***   Altamimi et al. 2017. "ITRF2014 plate motion model." Geophys. J. Int, 209.
***   Bird 2003. "An updated digital model of plate boundaries." Geochem. Geophys. Geosyst., 4(3), 1027.
***   DeMets et al. 2010. "Geologically current plate motions." Geophys. J. Int., 181:1-80.
***   Kreemer et al. 2014. "A geodetic plate motion and global strain rate model." Geochem. Geophys. Geosyst., 15, 3849-3889.
***   Mora-Páez et al. 2018. "Crustal deformation in the northern Andes – A new GPS velocity field." J. S. Amer. Earth Sci, 2018.
***   Snay 2003. "Introducing two spatial reference frames for regions of the Pacific Ocean." Surv. Land. Info. Sci., 63(1):5-12.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (NUMPLATE = 52)
      
      DIMENSION WX(NUMPLATE), WY(NUMPLATE), WZ(NUMPLATE)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      
* Main four plates included through version 3.2.9 (revised order):
** Rates for first three plates referenced to ITRF2008, fourth (MA) referenced to ITRF2000.
*** IPLATE = 1 (NA) North America (Altamimi et al. 2012)
***          2 (PA) Pacific (Altamimi et al. 2012)
***          3 (CA) Caribbean (Altamimi et al. 2012)
***          4 (MA) Mariana (Snay 2003)
***
* Plates added or edited in version 3.4.0:
** Rates for all plates referenced to ITRF2014.
***   Primary plates
***          5 (AF) Africa (Altamimi et al. 2017)
***          6 (AN) Antarctica (Altamimi et al. 2017)
***          7 (AU) Australia (Altamimi et al. 2017)
***          8 (EU) Eurasia (Altamimi et al. 2017)
***          9 (SA) South America (Altamimi et al. 2017)
***
***   Secondary plates
***         10 (AM) Amur (Kreemer et al. 2014)
***         11 (AR) Arabia (Altamimi et al. 2017)
***         12 (CO) Cocos (DeMets 2010)
***         13 (IN) India (Altamimi et al. 2017)
***         14 (JF) Juan de Fuca (DeMets 2010)
***         15 (NZ) Nazca (Altamimi et al. 2017)
***         16 (PS) Philippine Sea (Kreemer et al. 2014)
***         17 (SC) Scotia (DeMets 2010)
***         18 (SO) Somalia (Altamimi et al. 2017)
***         19 (SU) Sunda (Kreemer et al. 2014)
***
***   Tertiary plates
***         20 (AS) Aegean Sea (Kreemer et al. 2014)
***         21 (AP) Altiplano (Lamb 2000, in Bird 2003)
***         22 (AT) Anatolia (McClusky et al. 2000, in Bird 2003)
***         23 (BR) Balmoral Reef (Bird 2003)
***         24 (BS) Banda Sea (Rangin et al. 1999, in Bird 2003)
***         25 (BH) Birds Head (Bird 2003)
***         26 (BU) Burma (Circum-Pacific Map Project 1986, in Bird 2003)
***         27 (CL) Caroline (Seno et al. 1993, in Bird 2003)
***         28 (CR) Conway Reef (Bird 2003)
***         29 (EA) Easter (Engeln and Stein 1984, in Bird 2003)
***         30 (FT) Futuna (Bird 2003)
***         31 (GP) Galapagos (Lonsdale 1988, in Bird 2003)
***         32 (JZ) Juan Fernandez (Anderson-Fontana et al. 1986, in Bird 2003)
***         33 (KE) Kermadec (Bird 2003)
***         34 (MN) Manus (Martinez and Taylor 1996, in Bird 2003)
***         35 (MO) Maoke (Bird 2003)
***         36 (MS) Molucca Sea (Rangin et al. 1999, in Bird 2003)
***         37 (NH) New Hebrides (Bird 2003)
***         38 (NI) Niuafo'ou (Zellmer and Taylor 2001, in Bird 2003)
***         39 (ND) North Andes (Mora-Páez et al. 2018)
***         40 (NB) North Bismarck (Kreemer et al. 2014)
***         41 (OK) Okhotsk (Kreemer et al. 2014)
***         42 (ON) Okinawa (Kreemer et al. 2014)
***         43 (PM) Panama (Kreemer et al. 2014)
***         44 (RI) Rivera (DeMets 2010)
***         45 (SW) Sandwich (DeMets 2010)
***         46 (SL) Shetland (Kreemer et al. 2014)
***         47 (SS) Solomon Sea (Bird 2003)
***         48 (SB) South Bismarck (Kreemer et al. 2014)
***         49 (TI) Timor (Bird 2003)
***         50 (TO) Tonga (Kreemer et al. 2014)
***         51 (WL) Woodlark (Kreemer et al. 2014)
***         52 (YA) Yangtze (Kreemer et al. 2014)

* Plate rotation rates in radians per year. 
*** Plates #1-4 were defined in v3.2.9, with #1-3 referenced to ITRF2008 and
*** #4 referenced to ITRF2000. All other plates (#5-52) referenced to ITRF2014.

      DATA WX / 0.17000D-9, -1.99300D-9, 0.23800D-9, -0.09700D-9,     !North America-Mariana
     &          0.47997D-9, -1.20234D-9, 7.32069D-9, -0.41209D-9,     !Africa-Eurasia
     &         -1.30900D-9, -0.68964D-9, 5.59475D-9, -10.38028D-9,    !South America-Cocos
     &          5.59475D-9, 6.63589D-9, -1.61443D-9, 9.22141D-9,      !India-Philippine Sea
     &         -0.52297D-9, -0.58662D-9, -0.34003D-9, 1.26854D-9,     !Scotia-Aegean Sea
     &          0.05864D-9, 13.71303D-9, -2.85343D-9, -21.10731D-9,   !Altiplano-Banda Sea
     &         -1.79893D-9, 9.52310D-9, 1.73361D-9, -63.15807D-9,     !Birds Head-Conway Reef
     &         68.15282D-9, -85.23370D-9, 14.27334D-9, 106.03040D-9,  !Easter-Juan Fernandez
     &         31.33555D-9, -779.82645D-9, -0.46179D-9, 36.24013D-9,  !Kermadec-Molucca Sea
     &         42.92984D-9, -57.32448D-9, -1.96421D-9, -13.09213D-9,  !New Hebrides-North Bismarck
     &         -0.27495D-9, -15.50125D-9, 2.08835D-9, -21.93297D-9,   !Okhotsk-Rivera
     &         16.58716D-9, -8.64703D-9, -19.17916D-9, 97.26078D-9,   !Sandwich-South Bismarck
     &        -11.38292D-9, 137.84116D-9, -22.45705D-9, -1.01317D-9 / !Timor-Yangtze

      DATA WY /-3.20900D-9, 5.02300D-9, -5.27500D-9, 0.50900D-9,      !North America-Mariana
     &         -2.97676D-9, -1.57080D-9, 5.73050D-9, -2.57436D-9,     !Africa-Eurasia
     &         -1.45929D-9, -2.18428D-9, -0.65935D-9, -14.90060D-9,   !South America-Cocos
     &         -0.02424D-9, 11.76141D-9, -7.48552D-9, -4.96294D-9,    !India-Philippine Sea
     &         -1.79239D-9, -3.84942D-9, -3.69407D-9, 2.71543D-9,     !Scotia-Aegean Sea
     &         -8.07657D-9, 7.54290D-9, 2.80815D-9, 35.16250D-9,      !Altiplano-Banda Sea
     &         10.23283D-9, -39.44962D-9, 1.28481D-9, 10.29152D-9,    !Birds Head-Conway Reef
     &        165.61029D-9, 2.61244D-9, 94.43957D-9, 304.53677D-9,    !Easter-Juan Fernandez
     &          3.26279D-9, 445.94761D-9, 12.81479D-9, -53.21492D-9,  !Kermadec-Molucca Sea
     &         -4.47050D-9, -5.81369D-9, -1.51836D-9, 12.88081D-9,    !New Hebrides-North Bismarck
     &         -3.05687D-9, 10.47202D-9, -23.03737D-9, -70.43203D-9,  !Okhotsk-Rivera
     &        -11.88078D-9, 8.85561D-9, 22.26207D-9, -61.73877D-9,    !Sandwich-South Bismarck
     &         28.13885D-9, 10.44750D-9, 26.05666D-9, -2.26006D-9 /   !Timor-Yangtze

      DATA WZ /-0.48500D-9, -10.50100D-9, 3.21900D-9, -1.68200D-9,    !North America-Mariana
     &          3.55368D-9, 3.27249D-9, 5.89049D-9, 3.73307D-9,       !Africa-Eurasia
     &         -0.67874D-9, 4.19822D-9, 7.00071D-9, 9.13337D-9,       !South America-Cocos
     &          7.04919D-9, -10.62984D-9, 7.86853D-9, -11.55438D-9,   !India-Philippine Sea
     &          0.63488D-9, 4.28575D-9, 4.56571D-9, 3.07496D-9,       !Scotia-Aegean Sea
     &         -1.65936D-9, 13.29303D-9, -8.00888D-9, -0.28835D-9,    !Altiplano-Banda Sea
     &         -9.36606D-9, -3.31898D-9, -9.56706D-9, -24.27100D-9,   !Birds Head-Conway Reef
     &         83.81255D-9, -25.43833D-9, 4.51959D-9, 220.01253D-9,   !Easter-Juan Fernandez
     &         25.92570D-9, -57.95220D-9, 2.92132D-9, 3.16382D-9,     !Kermadec-Molucca Sea
     &          0.08496D-9, -3.72208D-9, 0.40012D-9, -10.73837D-9,    !New Hebrides-North Bismarck
     &          1.56236D-9, 14.79722D-9, 6.72874D-9, 27.07095D-9,     !Okhotsk-Rivera
     &        -12.18588D-9, 27.07392D-9, -1.89243D-9, 13.80350D-9,    !Sandwich-South Bismarck
     &         -1.68457D-9, 68.45806D-9, -1.16558D-9, 5.04044D-9 /    !Timor-Yangtze

      IF (IPLATE .LE. 0 .OR. IPLATE .GT. NUMPLATE) THEN
          WRITE (LUOUT, 1) IPLATE
    1     FORMAT(' Improper plate ID in PLATVL = ', I6)
          STOP
      ENDIF

      VX = -WZ(IPLATE) * Y + WY(IPLATE) * Z
      VY =  WZ(IPLATE) * X - WX(IPLATE) * Z
      VZ = -WY(IPLATE) * X + WX(IPLATE) * Y

*** The parameters for plate #4 (Mariana plate) refer to ITRF2000 (Snay 2003).
*** The following code converts  VX, VY, VZ from ITRF2000 to ITRF2008 velocities
*** for this plate.
      IF (IPLATE .EQ. 4) THEN
         VX = VX*1000.d0
         VY = VY*1000.d0
         VZ = VZ*1000.d0
         CALL VTRANF(X, Y, Z, VX, VY, VZ, 11, 15)
         VX = VX/1000.d0
         VY = VY/1000.d0
         VZ = VZ/1000.d0

*** The rotation rates for plates added after v3.2.9 refer to ITRF2014
*** (plates 5 through 52). The following code converts the rates to 
*** ITRF2008 velocities for these plates.
      ELSEIF (IPLATE .GE. 5) THEN
         VX = VX*1000.d0
         VY = VY*1000.d0
         VZ = VZ*1000.d0
         CALL VTRANF(X, Y, Z, VX, VY, VZ, 16, 15)
         VX = VX/1000.d0
         VY = VY/1000.d0
         VZ = VZ/1000.d0

*** The following translations rates are added per Altamimi et al. (2012)
*** for the first three plates referenced to ITRF2008.
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

*********************************************************
      SUBROUTINE TRFDAT(CARD,DATE,IREC12,IYEAR1,IYEAR2,MINS)

C Convert Bluebook date to time in minutes

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER    CARD*80
      CHARACTER   DATE*6
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

      IF (IREC12 .EQ. 0) THEN
         WRITE (LUOUT, 5)
    5    FORMAT(' ABORT: The Bluebook needs a valid *12* record.')
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
     1   ' The Bluebook file contains an observation that'/
     1   ' predates 1906.  The TDP model may not be valid'/
     1   ' and the computed corrections may be erroneous.')
      ENDIF
   
      IF (IDAY .EQ. 0) IDAY = 15

      IF (MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF

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

C Convert Bluebook date to time in years

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      character DATE*8
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      READ(DATE,10) IYEAR,MONTH,IDAY
   10 FORMAT(I4,I2,I2)
      IF(IYEAR .LE. 1906) THEN
        WRITE(LUOUT,20)
   20     FORMAT(' ***WARNING***'/
     1    ' The Bluebook file contains an observation that'/
     2    ' predates 1906.  The TDP model is not valid and the'/
     3    ' computed correction may be erroneous.')
      ENDIF
      IF(IDAY .EQ. 0) IDAY = 15
      IF(MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF

      CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
      MINS = MJD * 24 * 60
      RETURN
      END

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
     1       ' as deg-min-sec in free format (positive north).'/
     2       ' For example 35 17 28.3 or 35,17,28.3 '/)
      READ(LUIN,*,err=51,iostat=ios) LATD,LATM,SLAT
      if (ios /= 0) goto 51
      XLAT = (DBLE(3600*LATD + 60*LATM)+SLAT)/RHOSEC
      WRITE(LUOUT,130)
  130 FORMAT(' Specify the longitude for the origin of the line'/
     1       ' as deg-min-sec in free format (west positive).')
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
     1            MIN1, VNI, VEI, VUI)
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
*  Computes displacements at the point X1,X2 on the
*  Earth's surface due to 1.0 meter of right-lateral 
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

C********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,VERT,U1SS,U2SS,
     &                  U3SS,U1DS,U2DS,U3DS)

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
         F3 = 0.5D0*RATIO*(ETA/(R + DBAR) + YBAR*Q/
     &        ((R + DBAR)*(R + DBAR)) - DLOG(R + ETA))
         F1 = -0.5D0*RATIO*PSI*Q/((R + DBAR)*(R + DBAR))
      ELSE
         IF(DABS(PSI) .LE. 0.1D0) then
            F5 = 0.d0
         ELSE
            F5 = 2.D0*RATIO*DATAN((ETA*(X+Q*CDIP)+X*(R+X)*SDIP)/
     &           (PSI*(R+X)*CDIP))/CDIP
         ENDIF
         F4 = RATIO*(DLOG(R+DBAR)-SDIP*DLOG(R+ETA))/CDIP
         F3 = RATIO*(YBAR/(CDIP*(R+DBAR)) - DLOG(R+ETA)) + SDIP*F4/CDIP
         F1 = -RATIO*(PSI/(CDIP*(R+DBAR))) - SDIP*F5/CDIP
      ENDIF
      F2 = -RATIO*DLOG(R+ETA) - F3

      U1SS = -(PSI*Q/(R*(R+ETA)) + TERM + F1*SDIP)/TWOPI
      U2SS = -(YBAR*Q/(R*(R+ETA)) + Q*CDIP/(R+ETA) + F2*SDIP)/TWOPI
      U3SS = -(DBAR*Q/(R*(R+ETA)) + Q*SDIP/(R+ETA) + F4*SDIP)/TWOPI
      U1DS = -(Q/R - F3*SDIP*CDIP)/TWOPI

* Following part modified based on changes provided by Jarir Saleh 
* on 11/19/2020 to avoid floating point exceptions when a point is
* located on the extension of a fault plane.  
* Two conditional statements added to avoid floating point exception
* using variables TERM1 and TERM2: 
      TERM1 = YBAR*Q
      IF (ABS(TERM1) < 1D-3) THEN
        U2DS1 = 0.D0
      ELSE
        U2DS1 = -YBAR*Q/R/(R+PSI)/TWOPI
      ENDIF      

      U2DS2 = -CDIP*TERM/TWOPI
      U2DS3 = F1*SDIP*CDIP/TWOPI
      U2DS  = U2DS1 + U2DS2 + U2DS3
      
      TERM2 = DBAR*Q
      IF (ABS(TERM2) < 1D-3) THEN
        U3DS1 = 0.D0
      ELSE
        U3DS1 = -DBAR*Q/R/(R+PSI)/TWOPI
      ENDIF
      
      U3DS2 = -SDIP*TERM/TWOPI
      U3DS3 = F5*SDIP*CDIP/TWOPI
      U3DS  = U3DS1 + U3DS2 + U3DS3
* End of changes provided by Jarir

      RETURN
      END
      
C********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDWEI
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     RETURNS THE INDICES OF THE LOWER-LEFT
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
C PURPOSE:     RETRIEVES THE APPROXIMATE VALUES OF THE
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
      parameter (numref = 19)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER   CARD*80
      CHARACTER   record*120
      CHARACTER   NAMEF*80, NAME*80, NAMEBB*80, NAMEIF*80
      CHARACTER   NAME24*24
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
         write(luout,*) 'Do you wish to re-enter T1? (y/n)'
         read (luin,21,err=600,iostat=ios) ANSWER
         if (ios /= 0) goto 600
   21    format( A1 )
         IF (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 20
         RETURN
      ENDIF
      WRITE(LUOUT,*) 'Please enter T2 '
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
   50 FORMAT(A80)
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
      WRITE(LUOUT,80)                  !Remove label 75 (v3.3.0)
   80 FORMAT(' ********************************'/
     1   ' Displacements will be predicted at each point whose',/
     2   ' horizontal position is specified.',/
     3   ' Please indicate how you wish to supply positions.'/ )
   85 WRITE(LUOUT,86)
   86 FORMAT(
     1   '    0. No more points. Return to main menu.'/
     2   '    1. Individual points entered interactively.'/
     3   '    2. Points on a specified grid.'/
     4   '    3. The position records in a specified Bluebook file.'/
     5   '    4. Points on a specified line.  ' /
     6   '    5. File of delimited records of form LAT,LON,TEXT: ' /
     7   '       LAT = latitude in degrees (positive north)' /
     8   '       LON = longitude in degrees (positive west)' /
     9   '       TEXT = descriptive text (maximum 24 characters)' /
     1   '       Example:  40.731671553,112.212671753,SALT AIR' /)
     
      READ(LUIN,'(A1)',err=602,iostat=ios) OPTION
      if (ios /= 0) goto 602

      IF(OPTION .EQ. '0') THEN
        GO TO 510               
      ELSEIF(OPTION .EQ. '1') THEN
        CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1        LONDIR,NAME24, X,Y,Z,YLAT,YLON,EHT)
        ELON = -YLON
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
  150   FORMAT(' *****************************************'/
     1         ' Northward displacement = ',F7.3,' meters'/
     1         ' Eastward displacement  = ',F7.3,' meters'/
     1         ' Upward displacement    = ',F7.3,' meters'/
     1         ' *****************************************'//
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
  300   FORMAT(' Enter name of Bluebook file  ')
        READ(LUIN,310,err=603,iostat=ios) NAMEBB
        if (ios /= 0) goto 603
  310   FORMAT(A80)
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
     2  ' Bluebook file.  If you wish to calculate additional'/
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
  450   format(' Enter name of input file ')
        read(luin, 451,err=606,iostat=ios) NAMEIF
        if (ios /= 0) goto 606
  451   format(A80)
 
        open(I1,FILE=NAMEIF,STATUS='OLD')
        LINE = 0  !Inititalize line number of input file.
  455   read(I1,'(a)',END=460,err=607,iostat=ios) record
        if (ios /= 0) goto 607
        call interprate_latlon_record (record,XLAT,XLON,name24)
        LINE = LINE + 1   !Increment line number of input file.
        
*** If latitude magnitudes > 90 degrees or longitude magnitude > 360 degrees,
*** write error message and terminate program (added for v3.3.0).
        IF ((DABS(XLAT) .GT. 90).OR.(DABS(XLON) .GT. 360)) THEN
          WRITE(*,*)'***********************************************'
          WRITE(*,*)'Invalid latitude or longitude in input file'
          WRITE(*,4551)'on line', LINE, ':'
          WRITE(*,4552) XLAT, XLON, name24 
          WRITE(*,*)'Please check your input file and try again.'
          WRITE(*,*)'***********************************************'
 4551     FORMAT (1X, A, I8, A/)
 4552     FORMAT (1X, 2F17.10, 2X, A/)
          STOP
        ENDIF
        
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
      write (*,*) 'Wrong input file name in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  607 write (*,'(/)') 
      write (*,*) 'Wrong record from input in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
********************************************************************
      SUBROUTINE VELOC

*** Compute velocities at specified locations

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 19)
      parameter (nrsrch = 0)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER   CARD*80
      CHARACTER   NAMEF*80,NAME*80,NAMEBB*80,PVFILE*80
      CHARACTER   NAME24*24
      CHARACTER   BLAB*17
      CHARACTER   NAMEG*10
      CHARACTER   TYPE*4
      CHARACTER   OPTION*1,JN*1,JW*1,PVOUT*1
      CHARACTER   LATDIR*1, LONDIR*1
      character   frame1*24
      character   record*120

*** Commented out following temporary code (v3.3.0).
*** temporary code to plot results ****************
C       character TEMPNA*8   !Added string length to match string below.
C
C       TEMPNA = '        '
C       DUMMY = 0.0D0
*** end of temporary code *************************

      BLAB = 'OUTSIDE OF REGION'

      WRITE(LUOUT,10)
   10 FORMAT(' Please enter name for the file to contain the'/
     1       ' predicted velocities.  ')
      READ(LUIN,20,err=700,iostat=ios) NAMEF
      if (ios /= 0) goto 700
   20 FORMAT(A80)
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
            READ(LUIN,'(A80)',err=702,iostat=ios) PVFILE
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
 1002   FORMAT('VELOCITIES IN MM/YR RELATIVE TO ', a24 )
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
 1060   FORMAT(/A24,/
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
     4 '    3... The position records in a specified Bluebook file.'/
     5 '    4... Points on a specified line.  '/
     6 '    5... File of delimited records of form LAT,LON,TEXT:'/
     7 '         LAT = latitude in degrees (positive north) '/
     8 '         LON = longitude in degrees (positive west) '/
     9 '         TEXT = Descriptive text (maximum 24 characters) '/
     1 '         Example:  40.731671553,112.212671753,SALT AIR '/)

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
       LON = MINLON + J*JDS
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

*** Commented out following temporary code (v3.3.0).
*** temporary code to plot results
C         ZLAT = DBLE(LATD) + DBLE(LATM)/60.D0 + SLAT/3600.D0
C         ZLON = DBLE(LOND) + DBLE(LONM)/60.D0 + SLON/3600.D0
C         ZLON = 360.D0 - ZLON
C         TOTVEL = DSQRT(VN*VN + VE*VE)
C
*** code to create vectors to be plotted
C            IF (TOTVEL .GE. 5.0D0 .AND. TOTVEL .LE. 50.D0) THEN
C              TVE = VE/10.D0
C              TVN = VN/10.D0
C              WRITE(I2,129) ZLON,ZLAT,TVE,TVN,DUMMY,DUMMY,DUMMY,
C    1                   TEMPNA
C 129      FORMAT(F10.6, 2x, F9.6, 2X, 5(F7.2, 2x), A8)
C            ENDIF
C
*** code to create contour plots
C            WRITE(I2,129) ZLON, ZLAT, TOTVEL
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
  200   FORMAT(' Enter name of Bluebook file.  ')
        READ(LUIN,210,err=705,iostat=ios) NAMEBB
        if (ios /= 0) goto 705
  210   FORMAT(A80)
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
     2     ' Bluebook file.  If you wish to calculate additional'/
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
  410      FORMAT(A80)

           OPEN(I1,FILE=NAMEBB,STATUS='OLD')
           LINE = 0  !Inititalize line number of input file.
  420      READ(I1,'(a)',END=450,err=709,iostat=ios) record
           call interprate_latlon_record (record,XLAT,XLON,NAME24)
           if (ios /= 0) goto 709
           LINE = LINE + 1   !Increment line number of input file.
        
*** If latitude magnitudes > 90 degrees or longitude magnitude > 360 degrees,
*** write error message and terminate program (added for v3.3.0).
           IF ((DABS(XLAT) .GT. 90).OR.(DABS(XLON) .GT. 360)) THEN
             WRITE(*,*)'***********************************************'
             WRITE(*,*)'Invalid latitude or longitude in input file'
             WRITE(*,4201)'on line', LINE, ':'
             WRITE(*,4202) XLAT, XLON, name24 
             WRITE(*,*)'Please check your input file and try again.'
             WRITE(*,*)'***********************************************'
 4201        FORMAT (1X, A, I8, A/)
 4202        FORMAT (1X, 2F17.10, 2X, A/)
             STOP
           ENDIF

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
     2     ' file.  If you wish to calculate additional velocities, '/
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
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

      Is_iopt_NAD83 = (IOPT == 1)

*** Get reference latitude RLAT and reference longitude RLON in NAD 83
c     IF(IOPT .EQ. 0 .OR. IOPT .EQ. 1) THEN   !Not NAD83 anymore
      IF(IOPT .EQ. 15) THEN
          RLAT = YLAT
***     Added conditional statement to get rid of out-of-region error when
***     YLON is negative for reference frame (v3.3.0):
          IF(YLON .GE. 0.D0) THEN
            RLON = YLON
          ELSE
            RLON = YLON + TWOPI
          ENDIF
      ELSE
          ELON = -YLON
          CALL TOXYZ(YLAT,ELON,EHT,X,Y,Z)
          DATE = 2010.0d0
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

      RETURN
      END
 
******************************************************
      subroutine TRFPOS

*** Transform position between reference frames

      implicit double precision (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)
      parameter (nbbdim = 10000)
      parameter (rad2deg = 180.d0/3.14159265358979d0)

      character    card*80,namebb*80,nameif*80,name24*80
      character    record*120                    
      character    NAMEF*80
      character    frame1*24, frame2*24
      character    jn*1,jw*1,LATDIR*1,LONDIR*1
      character    option*1, answer*1, vopt*1
      character    PID*6,PIDs*6
      character    HTDP_version*8
      LOGICAL      FRMXYZ
      LOGICAL      TEST
      LOGICAL      Is_inp_NAD83, Is_out_NAD83
c      LOGICAL      Is_inp_NAD83PAC,Is_out_NAD83PAC    !Comment out v3.4.0
c      LOGICAL      Is_inp_NAD83MAR,Is_out_NAD83MAR    !Comment out v3.4.0

      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /ARRAYS/ HT(nbbdim), LOC(nbbdim),PIDs(nbbdim)
      COMMON /VERSION/ HTDP_version

C  Initialize NAD 83 input/poutput logical variables as false:
      Is_inp_NAD83 = .FALSE.
      Is_out_NAD83 = .FALSE.

      write(luout,80)
   80 format(
     1  ' Please enter the name for the file to contain'/
     2  ' the transformed positions.  ')
      read(luin,'(A80)',err=600,iostat=ios) namef
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

      IF(iopt1 .EQ. 1) Is_inp_NAD83 = .TRUE.

c      Is_inp_NAD83    = (iopt1 == 1)     !Comment out v3.4.0
c      Is_inp_NAD83PAC = (iopt1 == 12)    !Comment out v3.4.0
c      Is_inp_NAD83MAR = (iopt1 == 13)    !Comment out v3.4.0

  105 write(luout,110)
  110 format (/' Enter the reference frame for the output positions')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
         write(luout,*) ' Improper selection -- try again.  '
         go to 105
      endif

      IF(iopt2 .EQ. 1) Is_out_NAD83 = .TRUE.

c      Is_out_NAD83    = (iopt2 == 1)     !Comment out v3.4.0
c      Is_out_NAD83PAC = (iopt2 == 12)    !Comment out v3.4.0
c      Is_out_NAD83MAR = (iopt2 == 13)    !Comment out v3.4.0

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
     1  '    2... The position records in a specified Bluebook file.  '/
     1  '    3... Transform positions contained in file of delimited'/
     1  '         records of the form LAT,LON,EHT,TEXT:' /
     1  '         LAT = latitude in degrees (positive north)'/
     1  '         LON = longitude in degrees (positive west)'/
     1  '         EHT = ellipsoid height in meters'/
     1  '         TEXT = descriptive text (maximum 24 characters) '/
     1  '         Example: 40.731671553,112.212671753,34.241,SALT'/
     1  '    4... Transform positions contained in file of delimited'/
     1  '         records of the form X,Y,Z,TEXT, where:' /
     1  '         X, Y, Z = Cartesian coordinates in meters'/
     1  '         TEXT = descriptive text (maximum 24 characters) '/
     1  '         Example: -86682.104,-5394026.861,3391189.647,SALT'/
     1  '    5... Transform positions using velocities in a file of '/
     1  '         delimited records of the form X,Y,Z,Vx,Vy,Vz,TEXT' /
     1  '         X, Y, Z = Cartesian coordinates in meters'/
     1  '         Vx, Vy, Vz = Cartesian velocities in m/year'/
     1  '         TEXT = descriptive text (maximum 24 characters) '/
     1  '         Example:   '/
     1'   -86682.104,-5394026.861,3391189.647,0.012,0.006,-0.001,SALT'/)

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

          call GETVLY(ylat1,elon1,vx,vy,vz,vn,ve,vu,vopt,210)

          if (vopt .eq. '0') then
            call PREDV(ylat1,ylon1,eht1,date1,iopt1,jregn,vn,ve,vu)

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
      elseif (option .eq. '2') then     !Bluebook input and output
         vxsave = 0.d0
         vysave = 0.d0
         vzsave = 0.d0
         vnsave = 0.d0
         vesave = 0.d0
         vusave = 0.d0
         write(luout, 300)
  300    format(/ 
     1   ' Enter name of the Bluebook file:  ')
         read(luin, '(a)',err=603,iostat=ios) namebb
         if (ios /= 0) goto 603
         open (i1, file = namebb, status = 'old')

*** Obtaining the ellipsoid heights from the Bluebook file

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

*** Commented out due to unused label (v3.3.0).
*1309     format (/'The file, transformed_',a,
*     &   ', contains the transformed input file.'/)

         write (i2,1310) HTDP_version
 1310    format (' ***CAUTION: This file was processed using HTDP',
     &           ' version ',a8, '***')
         write (i2,1311) frame2       
 1311    format (' ***CAUTION: Coordinates in this file are in ',    
     &           a24, '***')
         WRITE (i2, 1312) MONTH2, IDAY2, IYEAR2, DATE2
 1312    FORMAT(' ***CAUTION: Coordinates in this file have been ',
     *   'updated to ',I2,'-',I2.2,'-',I4, '=(',F8.3,') ***'/)

*** Reread Bluebook file to get latitudes and longitudes
*** and to perform transformations
*** and to write transformed lat and lon to the *80* records of the output bfile
*** and to write transformed ellip. h to the *86* records of the output bfile

  310    read(i1, '(a80)', end = 390,err=604,iostat=ios) card
         if (ios /= 0) goto 604
         if (card(7:10) .eq. '*80*') then
           read(card,320,err=605,iostat=ios) isn,name24(1:30),latd,
     &                           latm,slat,jn,lond, lonm, slon, jw

           if (ios /= 0) goto 605
  320      format(10x, i4, a30, i2, i2, f7.5, a1,
     &            i3, i2, f7.5, a1)

*** Commented out due to unused label (v3.3.0).
*  321      format(i2, i2, f7.5, a1,3x, i3, i2, f7.5, a1)

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
         write (i2,1310) HTDP_version
         write (i2,1311) frame2       
         WRITE (I2, 1312) MONTH2, IDAY2, IYEAR2, DATE2

         LINE = 0  !Inititalize line number of input file.
  410    read (i1,'(a)',end=450,err=607,iostat=ios) record
         if (ios /= 0) goto 607
         call interprate_XYZ_record (record,xlat,xlon,eht,name24)
         LINE = LINE + 1   !Increment line number of input file.
        
*** If latitude magnitudes > 90 degrees or longitude magnitude > 360 degrees,
*** write error message and terminate program (added for v3.3.0).
         IF ((DABS(xlat) .GT. 90).OR.(DABS(xlon) .GT. 360)) THEN
           WRITE(*,*)'***********************************************'
           WRITE(*,*)'Invalid latitude or longitude in input file'
           WRITE(*,4101)'on line', LINE, ':'
           WRITE(*,4102) xlat, xlon, eht, name24 
           WRITE(*,*)'Please check your input file and try again.'
           WRITE(*,*)'***********************************************'
 4101      FORMAT (1X, A, I8, A/)
 4102      FORMAT (1X, 2F17.10, F11.3, 2X, A)
           STOP
         ENDIF
         
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

         write (i2,449) outlat,outlon,ehtnew,name24(1:iii)

c        write (i6,449) outlat,outlon,ehtnew,trim(name24)
  449    format (2f16.10,f10.3,4x,a)
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
 4001    format (' Enter name of input file: ')
         read (luin, '(a)',err=600,iostat=ios) nameif
         if (ios /= 0) goto 600
         open (i1, file = nameif, status = 'old')
c        open (i6,file='transformed_'//nameif)

C  write some comments in the transformed files

         call extract_name (nameif,iii)
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

         write (i2,1449) xt,yt,zt,name24(1:iii)

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

         write (i2,1449) xt,yt,zt,name24(1:iii)

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
C  The parameters in common block tranpa are computed using the IGS convention ITRF96 <> ITRF97.
C  The parameters in common block tranpa1 are computed using the IERS convention ITRF96 = ITRF97.
C  The latter parameters were added to HTDP in 09/2014. They will be used to transform between
C  ITRF systems. They will not be used if the transformation involves NAD83(2011).
C  However, the NAD 83 Pacific and Mariana frames use the IERS convention ITRF96 = ITRF97.
C  Removed equivalence between NAD83(2011)and WGS84 original (Transit) in v3.4.0.

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)

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

C  Parameters computed with the IGS convention ITRF96 <> ITRF97

*** From ITRF94 to NAD 83
*** In v3.4.0, changed rotations from radians to mas to make consistent
*** with published parameters.
      tx(1) =     0.9910d0        
      ty(1) =    -1.9072d0      
      tz(1) =    -0.5129d0      
      dtx(1) =    0.d0        
      dty(1) =    0.d0       
      dtz(1) =    0.d0      
      rx(1) =     0.02579d0 / rhosec     !Previous  1.25033d-7 rad (prior to v3.4.0)
      ry(1) =     0.00965d0 / rhosec     !Previous  0.46785d-7 rad (prior to v3.4.0)
      rz(1) =     0.01166d0 / rhosec     !Previous  0.56529d-7 rad (prior to v3.4.0)
      drx(1) =    0.0000532d0 / rhosec   !Previous  0.00258d-7 rad (prior to v3.4.0)
      dry(1) =   -0.0007423d0 / rhosec   !Previous -0.03599d-7 rad (prior to v3.4.0)
      drz(1) =   -0.0000316d0 / rhosec   !Previous -0.00153d-7 rad (prior to v3.4.0)
      scale(1) =  0.d0
      dscale(1) = 0.0d0
      refepc(1) = 1997.0d0

*** From ITRF94 to ITRF88
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx(2) =     0.018d0
      ty(2) =     0.000d0
      tz(2) =    -0.092d0
      dtx(2) =    0.0d0    
      dty(2) =    0.0d0
      dtz(2) =    0.0d0
      rx(2) =    -0.0001d0 / rhosec
      ry(2) =     0.0d0 / rhosec
      rz(2) =     0.0d0 / rhosec
      drx(2) =    0.0d0 / rhosec
      dry(2) =    0.0d0 / rhosec
      drz(2) =    0.0d0 / rhosec
      scale(2) =  0.749d-8    !Previous 0.74d-8 (prior to v3.3.0)
      dscale(2) = 0.0d0
      refepc(2) = 1988.0d0

*** From ITRF94 to ITRF89
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx(3) =     0.023d0
      ty(3) =     0.036d0
      tz(3) =    -0.068d0
      dtx(3) =    0.0d0    
      dty(3) =    0.0d0
      dtz(3) =    0.0d0
      rx(3) =     0.0d0 / rhosec
      ry(3) =     0.0d0 / rhosec
      rz(3) =     0.0d0 / rhosec
      drx(3) =    0.0d0 / rhosec
      dry(3) =    0.0d0 / rhosec
      drz(3) =    0.0d0 / rhosec
      scale(3) =  0.439d-8    !Previous 0.43d-8 (prior to v3.3.0)
      dscale(3) = 0.0d0
      refepc(3) = 1988.0d0

*** From ITRF94 to ITRF90
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx(4) =     0.018d0
      ty(4) =     0.012d0
      tz(4) =    -0.030d0
      dtx(4) =    0.0d0    
      dty(4) =    0.0d0
      dtz(4) =    0.0d0
      rx(4) =     0.0d0 / rhosec
      ry(4) =     0.0d0 / rhosec
      rz(4) =     0.0d0 / rhosec
      drx(4) =    0.0d0 / rhosec
      dry(4) =    0.0d0 / rhosec
      drz(4) =    0.0d0 / rhosec
      scale(4) =  0.099d-8    !Previous 0.09d-8 (prior to v3.3.0)
      dscale(4) = 0.0d0
      refepc(4) = 1988.0d0

*** From ITRF94 to ITRF91
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx(5) =     0.020d0
      ty(5) =     0.016d0
      tz(5) =    -0.014d0
      dtx(5) =    0.0d0    
      dty(5) =    0.0d0
      dtz(5) =    0.0d0
      rx(5) =     0.0d0 / rhosec
      ry(5) =     0.0d0 / rhosec
      rz(5) =     0.0d0 / rhosec
      drx(5) =    0.0d0 / rhosec
      dry(5) =    0.0d0 / rhosec
      drz(5) =    0.0d0 / rhosec
      scale(5) =  0.069d-8    !Previous 0.06d-8 (prior to v3.3.0)
      dscale(5) = 0.0d0
      refepc(5) = 1988.0d0

*** From ITRF94 to ITRF92
***   Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx(6) =     0.008d0
      ty(6) =     0.002d0
      tz(6) =    -0.008d0
      dtx(6) =    0.0d0    
      dty(6) =    0.0d0
      dtz(6) =    0.0d0
      rx(6) =     0.0d0 / rhosec
      ry(6) =     0.0d0 / rhosec
      rz(6) =     0.0d0 / rhosec
      drx(6) =    0.0d0 / rhosec
      dry(6) =    0.0d0 / rhosec
      drz(6) =    0.0d0 / rhosec
      scale(6) = -0.071d-8    !Previous -0.08d-8 (prior to v3.3.0)
      dscale(6) = 0.0d0
      refepc(6) = 1988.0d0

*** From ITRF94 to ITRF93
*** Update scale in v3.3.0 to match current value published by IGN 
*** (2022) "Transformation parameters from ITRF2020 to past ITRFs" at
*** https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt
*** This update also affects the scale for transformations from ITRF94 
*** to all earlier ITRFs in HTDP.
*** Note that this transformation is not published in IERS TN20 (1996).
      tx(7) =     0.006d0
      ty(7) =    -0.005d0
      tz(7) =    -0.015d0
      dtx(7) =   -0.0029d0
      dty(7) =    0.0004d0
      dtz(7) =    0.0008d0
      rx(7) =     0.00039d0 / rhosec
      ry(7) =    -0.00080d0 / rhosec
      rz(7) =     0.00096d0 / rhosec
      drx(7) =    0.00011d0 / rhosec
      dry(7) =    0.00019d0 / rhosec
      drz(7) =   -0.00005d0 / rhosec
      scale(7) =  0.049d-8    !Previous 0.040d-8 (prior to v3.3.0)
      dscale(7) = 0.0d0
      refepc(7) = 1988.0d0

*** From ITRF94 to ITRF96
      tx(8) =     0.0d0
      ty(8) =     0.0d0
      tz(8) =     0.0d0
      dtx(8) =    0.0d0
      dty(8) =    0.0d0
      dtz(8) =    0.0d0
      rx(8) =     0.0d0 / rhosec
      ry(8) =     0.0d0 / rhosec
      rz(8) =     0.0d0 / rhosec
      drx(8) =    0.0d0 / rhosec
      dry(8) =    0.0d0 / rhosec
      drz(8) =    0.0d0 / rhosec
      scale(8) =  0.0d0
      dscale(8) = 0.0d0
      refepc(8) = 1997.0d0    !Previous 1996.0d0 (prior to v3.4.0)

*** From ITRF94 to ITRF97, based on IGS adopted values
*** According to IGS:  ITRF97 <> ITRF96 = ITRF94
      tx(9) =     0.00207d0
      ty(9) =     0.00021d0
      tz(9) =    -0.00995d0
      dtx(9) =   -0.00069d0
      dty(9) =    0.00010d0
      dtz(9) =   -0.00186d0
      rx(9) =    -0.00012467d0 / rhosec
      ry(9) =     0.00022355d0 / rhosec
      rz(9) =     0.00006065d0 / rhosec
      drx(9) =   -0.00001347d0 / rhosec
      dry(9) =    0.00001514d0 / rhosec
      drz(9) =   -0.00000027d0 / rhosec
      scale(9) =  0.93496d-9
      dscale(9) = 0.19201d-9
      refepc(9) = 1997.0d0

*** From ITRF94 to WGS72 (composition of ITRF94 -> NAD_83 -> WGS72)
*** Remove WGS72 and replace with new WGS84 original (Transit) in v3.4.0
c      tx(10) = 0.9910d0
c      ty(10) = -1.9072d0
c      tz(10) = -0.5129d0 - 4.5d0
c      dtx(10) = 0.d0
c      dty(10) = 0.d0
c      dtz(10) = 0.d0
c      rx(10) = 1.25033d-7
c      ry(10) = 0.46785d-7
c      rz(10) = 0.56529d-7 + 26.85868d-7
c      drx(10) = 0.00258d-7
c      dry(10) = -0.03599d-7
c      drz(10) = -0.00153d-7
c      scale(10) = 0.d0 - 0.2263d-6
c      dscale(10) = 0.0d0
c      refepc(10) = 1997.0d0
      
*** From ITRF94 to WGS84 original (Transit), added in v3.4.0.
*** Based on IERS-published transformation from ITRF90 to WGS84 (Transit)
*** in Table 3.1 of McCarthy (1992) "IERS standards," IERS Tech. Note 13, 
*** Observatoire de Paris, France: Central Bureau of IERS. 
      tx(10) =      0.078d0
      ty(10) =     -0.505d0
      tz(10) =     -0.253d0
      dtx(10) =     0.0d0
      dty(10) =     0.0d0
      dtz(10) =     0.0d0
      rx(10) =     -0.0183d0 / rhosec
      ry(10) =      0.0003d0 / rhosec
      rz(10) =     -0.0070d0 / rhosec
      drx(10) =     0.0d0 / rhosec
      dry(10) =     0.0d0 / rhosec
      drz(10) =     0.0d0 / rhosec
      scale(10) = -10.010d-9
      dscale(10) =  0.0d0
      refepc(10) =  1997.0d0

*** From ITRF94 to ITRF2000 (also IGS00 and IGb00)
*** assumes that ITRF94 = ITRF96 and
*** uses IGS values for ITRF96 -> ITRF97
*** and IERS values for ITRF97 -> ITRF2000
      tx(11)     = -0.00463d0
      ty(11)     = -0.00589d0
      tz(11)     =  0.00855d0
      dtx(11)    = -0.00069d0
      dty(11)    =  0.00070d0
      dtz(11)    = -0.00046d0
      rx(11)     = -0.00012467d0 / rhosec
      ry(11)     =  0.00022355d0 / rhosec
      rz(11)     =  0.00006065d0 / rhosec
      drx(11)    = -0.00001347d0 / rhosec
      dry(11)    =  0.00001514d0 / rhosec
      drz(11)    =  0.00001973d0 / rhosec
      scale(11)  = -0.61504d-9
      dscale(11) =  0.18201d-9
      refepc(11) =  1997.0d0

*** From ITRF94 to PACP00
*** use PA/ITRF2000 rotation rates from Beavan et al. (2002)
*** In v3.4.0, edited parameters to make consistent with parameters used
*** for ITRF94 to ITRF97 transformation which differed due to rounding.
      tx(12)     =  0.90557d0               !Previous  0.9056d0 (prior to v3.4.0)
      ty(12)     = -2.01999d0               !Previous -2.0200d0 (prior to v3.4.0)
      tz(12)     = -0.55165d0               !Previous -0.5516d0 (prior to v3.4.0)
      dtx(12)    = -0.00069d0
      dty(12)    =  0.00070d0
      dtz(12)    = -0.00046d0
      rx(12)     =  0.02761633d0 / rhosec   !Previous  0.027616d0 (prior to v3.4.0)
      ry(12)     =  0.01369255d0 / rhosec   !Previous  0.013692d0 (prior to v3.4.0)
      rz(12)     =  0.00277265d0 / rhosec   !Previous  0.002773d0 (prior to v3.4.0)
      drx(12)    = -0.00039747d0 / rhosec   !Previous -0.000397d0 (prior to v3.4.0)
      dry(12)    =  0.00102214d0 / rhosec   !Previous  0.001022d0 (prior to v3.4.0)
      drz(12)    = -0.00216627d0 / rhosec   !Previous -0.002166d0 (prior to v3.4.0)
      scale(12)  = -0.61504d-9
      dscale(12) =  0.18201d-9
      refepc(12) =  1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
*** In v3.4.0, edited parameters to make consistent with parameters used
*** for ITRF94 to ITRF97 transformation which differed due to rounding.
      tx(13)     =  0.90557d0               !Previous  0.9056d0 (prior to v3.4.0)
      ty(13)     = -2.01999d0               !Previous -2.0200d0 (prior to v3.4.0)
      tz(13)     = -0.55165d0               !Previous -0.5516d0 (prior to v3.4.0)
      dtx(13)    = -0.00069d0
      dty(13)    =  0.00070d0
      dtz(13)    = -0.00046d0
      rx(13)     =  0.02884633d0 / rhosec   !Previous  0.028847d0 (prior to v3.4.0)
      ry(13)     =  0.01064355d0 / rhosec   !Previous  0.010644d0 (prior to v3.4.0)
      rz(13)     =  0.00898865d0 / rhosec   !Previous  0.008989d0 (prior to v3.4.0)
      drx(13)    = -0.00003347d0 / rhosec   !Previous -0.000033d0 (prior to v3.4.0)
      dry(13)    =  0.00012014d0 / rhosec   !Previous  0.000120d0 (prior to v3.4.0)
      drz(13)    = -0.00032727d0 / rhosec   !Previous -0.000327d0 (prior to v3.4.0)
      scale(13)  = -0.61504d-9
      dscale(13) =  0.18201d-9
      refepc(13) =  1997.00d0

*** From ITRF94 to ITRF2005 (also IGS05)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
      tx(14)     = -0.00533d0
      ty(14)     = -0.00479d0
      tz(14)     =  0.00895d0
      dtx(14)    = -0.00049d0
      dty(14)    =  0.00060d0
      dtz(14)    =  0.00134d0
      rx(14)     = -0.00012467d0 / rhosec
      ry(14)     =  0.00022355d0 / rhosec
      rz(14)     =  0.00006065d0 / rhosec
      drx(14)    = -0.00001347d0 / rhosec
      dry(14)    =  0.00001514d0 / rhosec
      drz(14)    =  0.00001973d0 / rhosec
      scale(14)  = -0.77504d-9
      dscale(14) =  0.10201d-9
      refepc(14) =  1997.0d0

*** From ITRF94 to ITRF2008 (also IGS08 and IGb08)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
      tx(15)     = -0.00243d0
      ty(15)     = -0.00389d0
      tz(15)     =  0.01365d0
      dtx(15)    = -0.00079d0
      dty(15)    =  0.00060d0
      dtz(15)    =  0.00134d0
      rx(15)     = -0.00012467d0 / rhosec
      ry(15)     =  0.00022355d0 / rhosec
      rz(15)     =  0.00006065d0 / rhosec
      drx(15)    = -0.00001347d0 / rhosec
      dry(15)    =  0.00001514d0 / rhosec
      drz(15)    =  0.00001973d0 / rhosec
      scale(15)  = -1.71504d-9
      dscale(15) =  0.10201d-9
      refepc(15) =  1997.0d0

*** From ITRF94 to ITRF2014 (also IGS14 and IGb14)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
*** uses IERS values for ITRF2008 -> ITRF2014
      tx(16)     = -0.01430d0
      ty(16)     =  0.00201d0
      tz(16)     =  0.02867d0
      dtx(16)    = -0.00079d0
      dty(16)    =  0.00060d0
      dtz(16)    =  0.00144d0
      rx(16)     = -0.00029978d0 / rhosec
      ry(16)     =  0.00042037d0 / rhosec
      rz(16)     =  0.00031714d0 / rhosec
      drx(16)    = -0.00001347d0 / rhosec
      dry(16)    =  0.00001514d0 / rhosec
      drz(16)    =  0.00001973d0 / rhosec
      scale(16)  = -0.36891d-9
      dscale(16) =  0.07201d-9
      refepc(16) =  2010.0d0

*** From ITRF94 to ITRF2020 (also IGS20)
*** assumes that ITRF94 = ITRF96
*** uses IGS values for ITRF96 -> ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000-> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
*** uses IERS values for ITRF2008 -> ITRF2014
*** uses IERS values for ITRF2014 -> ITRF2020
      tx(17)     = -0.01290d0
      ty(17)     =  0.00241d0
      tz(17)     =  0.02827d0
      dtx(17)    = -0.00079d0
      dty(17)    =  0.00070d0
      dtz(17)    =  0.00124d0
      rx(17)     = -0.00029978d0 / rhosec
      ry(17)     =  0.00042037d0 / rhosec
      rz(17)     =  0.00031714d0 / rhosec
      drx(17)    = -0.00001347d0 / rhosec
      dry(17)    =  0.00001514d0 / rhosec
      drz(17)    =  0.00001973d0 / rhosec
      scale(17)  =  0.05109d-9
      dscale(17) =  0.07201d-9
      refepc(17) =  2010.0d0

*** From ITRF94 to WGS84(G1150), added in v3.5.0.
*** Uses IGS values for ITRF96 -> ITRF97.
*** Treated as biased with respect to ITRF2000.
*** Based on NGA-published transformation between WGS84(G1150)
*** and WGS84(G1762) = ITRF2008 in Table 2.5 of NGA.STND.0036_1.0.0_WGS84,
*** "Department of Defense World Geodetic System 1984: its definition and 
*** relationships with local geodetic systems", v1.0.0, 2014.  See also
*** Kelly and Dennis (2022) "Transforming Between WGS84 Realizations,"
*** J Surv Eng, doi: 10.1061/(ASCE)SU.1943-5428.0000389.
      tx(18)     = -0.00580d0
      ty(18)     = -0.00019d0
      tz(18)     = -0.00513d0
      dtx(18)    = -0.00069d0
      dty(18)    =  0.00070d0
      dtz(18)    = -0.00046d0
      rx(18)     = -0.00029978d0 / rhosec
      ry(18)     =  0.00042037d0 / rhosec
      rz(18)     =  0.00031714d0 / rhosec
      drx(18)    = -0.00001347d0 / rhosec
      dry(18)    =  0.00001514d0 / rhosec
      drz(18)    =  0.00001973d0 / rhosec
      scale(18)  =  4.83109d-9
      dscale(18) =  0.18201d-9
      refepc(18) =  2010.0d0

*** From ITRF94 to WGS84(G1674), added in v3.5.0.
*** Uses IGS values for ITRF96 -> ITRF97.
*** Treated as biased with respect to ITRF2008.
*** Based on NGA-published bias transformation between WGS84(G1674)
*** and WGS84(G1762) = ITRF2008 in Table 2.5 of NGA.STND.0036_1.0.0_WGS84,
*** "Department of Defense World Geodetic System 1984: its definition and 
*** relationships with local geodetic systems", v1.0.0, 2014.
*** See also Kelly and Dennis (2022) "Transforming Between WGS84 Realizations,"
*** J Surv Eng, doi: 10.1061/(ASCE)SU.1943-5428.0000389.
      tx(19)     = -0.00870d0
      ty(19)     =  0.00091d0
      tz(19)     =  0.02707d0
      dtx(19)    = -0.00079d0
      dty(19)    =  0.00060d0
      dtz(19)    =  0.00134d0
      rx(19)     = -0.00056978d0 / rhosec
      ry(19)     =  0.00069037d0 / rhosec
      rz(19)     = -0.00006286d0 / rhosec
      drx(19)    = -0.00001347d0 / rhosec
      dry(19)    =  0.00001514d0 / rhosec
      drz(19)    =  0.00001973d0 / rhosec
      scale(19)  =  6.51109d-9
      dscale(19) =  0.10201d-9
      refepc(19) =  2010.0d0

C*************************************************************************************************************************
C  Parameters computed with the IERS convention ITRF96 = ITRF97

*** From ITRF94 to NAD 83
*** In v3.4.0, changed rotations from radians to mas to make consistent
*** with published parameters.
      tx1(1) =     0.9910d0        
      ty1(1) =    -1.9072d0      
      tz1(1) =    -0.5129d0      
      dtx1(1) =    0.d0        
      dty1(1) =    0.d0       
      dtz1(1) =    0.d0      
      rx1(1) =     0.02579d0 / rhosec     !Previous  1.25033d-7 rad (prior to v3.4.0)
      ry1(1) =     0.00965d0 / rhosec     !Previous  0.46785d-7 rad (prior to v3.4.0)
      rz1(1) =     0.01166d0 / rhosec     !Previous  0.56529d-7 rad (prior to v3.4.0)
      drx1(1) =    0.0000532d0 / rhosec   !Previous  0.00258d-7 rad (prior to v3.4.0)
      dry1(1) =   -0.0007423d0 / rhosec   !Previous -0.03599d-7 rad (prior to v3.4.0)
      drz1(1) =   -0.0000316d0 / rhosec   !Previous -0.00153d-7 rad (prior to v3.4.0)
      scale1(1) =  0.d0
      dscale1(1) = 0.0d0
      refepc1(1) = 1997.0d0

*** From ITRF94 to ITRF88
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx1(2) =     0.018d0
      ty1(2) =     0.000d0
      tz1(2) =    -0.092d0
      dtx1(2) =    0.0d0    
      dty1(2) =    0.0d0
      dtz1(2) =    0.0d0
      rx1(2) =    -0.0001d0 / rhosec
      ry1(2) =     0.0d0 / rhosec
      rz1(2) =     0.0d0 / rhosec
      drx1(2) =    0.0d0 / rhosec
      dry1(2) =    0.0d0 / rhosec
      drz1(2) =    0.0d0 / rhosec
      scale1(2) =  0.749d-8    !Previous 0.74d-8 (prior to v3.3.0)
      dscale1(2) = 0.0d0
      refepc1(2) = 1988.0d0

*** From ITRF94 to ITRF89
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx1(3) =     0.023d0
      ty1(3) =     0.036d0
      tz1(3) =    -0.068d0
      dtx1(3) =    0.0d0    
      dty1(3) =    0.0d0
      dtz1(3) =    0.0d0
      rx1(3) =     0.0d0 / rhosec
      ry1(3) =     0.0d0 / rhosec
      rz1(3) =     0.0d0 / rhosec
      drx1(3) =    0.0d0 / rhosec
      dry1(3) =    0.0d0 / rhosec
      drz1(3) =    0.0d0 / rhosec
      scale1(3) =  0.439d-8    !Previous 0.43d-8 (prior to v3.3.0)
      dscale1(3) = 0.0d0
      refepc1(3) = 1988.0d0

*** From ITRF94 to ITRF90
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx1(4) =     0.018d0
      ty1(4) =     0.012d0
      tz1(4) =    -0.030d0
      dtx1(4) =    0.0d0    
      dty1(4) =    0.0d0
      dtz1(4) =    0.0d0
      rx1(4) =     0.0d0 / rhosec
      ry1(4) =     0.0d0 / rhosec
      rz1(4) =     0.0d0 / rhosec
      drx1(4) =    0.0d0 / rhosec
      dry1(4) =    0.0d0 / rhosec
      drz1(4) =    0.0d0 / rhosec
      scale1(4) =  0.099d-8    !Previous 0.09d-8 (prior to v3.3.0)
      dscale1(4) = 0.0d0
      refepc1(4) = 1988.0d0

*** From ITRF94 to ITRF91
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx1(5) =     0.020d0
      ty1(5) =     0.016d0
      tz1(5) =    -0.014d0
      dtx1(5) =    0.0d0    
      dty1(5) =    0.0d0
      dtz1(5) =    0.0d0
      rx1(5) =     0.0d0 / rhosec
      ry1(5) =     0.0d0 / rhosec
      rz1(5) =     0.0d0 / rhosec
      drx1(5) =    0.0d0 / rhosec
      dry1(5) =    0.0d0 / rhosec
      drz1(5) =    0.0d0 / rhosec
      scale1(5) =  0.069d-8    !Previous 0.06d-8 (prior to v3.3.0)
      dscale1(5) = 0.0d0
      refepc1(5) = 1988.0d0

*** From ITRF94 to ITRF92
*** Update scale in v3.3.0 due to ITRF94 to ITRF93 transformation update.
      tx1(6) =     0.008d0
      ty1(6) =     0.002d0
      tz1(6) =    -0.008d0
      dtx1(6) =    0.0d0    
      dty1(6) =    0.0d0
      dtz1(6) =    0.0d0
      rx1(6) =     0.0d0 / rhosec
      ry1(6) =     0.0d0 / rhosec
      rz1(6) =     0.0d0 / rhosec
      drx1(6) =    0.0d0 / rhosec
      dry1(6) =    0.0d0 / rhosec
      drz1(6) =    0.0d0 / rhosec
      scale1(6) = -0.071d-8    !Previous -0.08d-8 (prior to v3.3.0)
      dscale1(6) = 0.0d0
      refepc1(6) = 1988.0d0

*** From ITRF94 to ITRF93
*** Update scale in v3.3.0 to match current value published by IGN 
*** (2022) "Transformation parameters from ITRF2020 to past ITRFs" at
*** https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt
*** This update also affects the scale for transformations from ITRF94 
*** to all earlier ITRFs in HTDP.
*** Note that this transformation is not published in IERS TN20 (1996).
      tx1(7)     =  0.006d0
      ty1(7)     = -0.005d0
      tz1(7)     = -0.015d0
      dtx1(7)    = -0.0029d0
      dty1(7)    =  0.0004d0
      dtz1(7)    =  0.0008d0
      rx1(7)     =  0.00039d0 / rhosec
      ry1(7)     = -0.00080d0 / rhosec
      rz1(7)     =  0.00096d0 / rhosec
      drx1(7)    =  0.00011d0 / rhosec
      dry1(7)    =  0.00019d0 / rhosec
      drz1(7)    = -0.00005d0 / rhosec
      scale1(7)  =  0.049d-8    !Previous 0.04d-8 (prior to v3.3.0)
      dscale1(7) =  0.0d0
      refepc1(7) =  1988.0d0

*** From ITRF94 to ITRF96
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx1(8)     = 0.d0
      ty1(8)     = 0.d0
      tz1(8)     = 0.d0
      dtx1(8)    = 0.d0
      dty1(8)    = 0.d0
      dtz1(8)    = 0.d0
      rx1(8)     = 0.d0 / rhosec
      ry1(8)     = 0.d0 / rhosec
      rz1(8)     = 0.d0 / rhosec
      drx1(8)    = 0.d0 / rhosec
      dry1(8)    = 0.d0 / rhosec
      drz1(8)    = 0.d0 / rhosec
      scale1(8)  = 0.d0
      dscale1(8) = 0.0d0
      refepc1(8) = 1997.0d0    !Previous 1996.0d0 (prior to v3.4.0)

*** From ITRF94 to ITRF97 (also IGS97), based on IERS adopted values
*** According to IERS:  ITRF97 = ITRF96 = ITRF94
      tx1(9)     = 0.0d0
      ty1(9)     = 0.0d0
      tz1(9)     = 0.0d0
      dtx1(9)    = 0.0d0
      dty1(9)    = 0.0d0
      dtz1(9)    = 0.0d0
      rx1(9)     = 0.0d0 / rhosec
      ry1(9)     = 0.0d0 / rhosec
      rz1(9)     = 0.0d0 / rhosec
      drx1(9)    = 0.0d0 / rhosec
      dry1(9)    = 0.0d0 / rhosec
      drz1(9)    = 0.0d0 / rhosec
      scale1(9)  = 0.0d0
      dscale1(9) = 0.0d0
      refepc1(9) = 1997.0d0    !Previous 2000.0d0 (prior to v3.4.0)

*** From ITRF94 to WGS72 (composition of ITRF94 -> NAD_83 -> WGS72)
*** Remove WGS72 and replace with new WGS84 original (Transit) in v3.4.0
c      tx1(10) = 0.9910d0
c      ty1(10) = -1.9072d0
c      tz1(10) = -0.5129d0 - 4.5d0
c      dtx1(10) = 0.d0
c      dty1(10) = 0.d0
c      dtz1(10) = 0.d0
c      rx1(10) = 1.25033d-7
c      ry1(10) = 0.46785d-7
c      rz1(10) = 0.56529d-7 + 26.85868d-7
c      drx1(10) = 0.00258d-7
c      dry1(10) = -0.03599d-7
c      drz1(10) = -0.00153d-7
c      scale1(10) =  0.d0 - 0.2263d-6
c      dscale1(10) = 0.0d0
c      refepc1(10) = 1997.0d0
      
*** From ITRF94 to WGS84 original (Transit), added in v3.4.0.
*** Based on IERS-published transformation from ITRF90 to WGS84 (Transit)
*** in Table 3.1 of McCarthy (1992) "IERS standards," IERS Tech. Note 13, 
*** Observatoire de Paris, France: Central Bureau of IERS. 
      tx1(10) =      0.078d0
      ty1(10) =     -0.505d0
      tz1(10) =     -0.253d0
      dtx1(10) =     0.0d0
      dty1(10) =     0.0d0
      dtz1(10) =     0.0d0
      rx1(10) =     -0.0183d0 / rhosec
      ry1(10) =      0.0003d0 / rhosec
      rz1(10) =     -0.0070d0 / rhosec
      drx1(10) =     0.0d0    / rhosec
      dry1(10) =     0.0d0    / rhosec
      drz1(10) =     0.0d0    / rhosec
      scale1(10) = -10.010d-9
      dscale1(10) =  0.0d0
      refepc1(10) =  1997.0d0

*** From ITRF94 to ITRF2000 (also IGS00 and IGb00)
*** assumes that         ITRF94 = ITRF96 and
*** uses IERS convention ITRF96 = ITRF97
*** and  IERS values for ITRF97 -> ITRF2000
*** Changed epoch to 1997.0 to make consistent with trapa epoch (v3.4.0)
      tx1(11)     = -0.0067d0
      ty1(11)     = -0.0061d0             !Previous -0.00430d0 (prior to v3.4.0)
      tz1(11)     =  0.0185d0             !Previous  0.02270d0 (prior to v3.4.0)
      dtx1(11)    =  0.0000d0
      dty1(11)    =  0.0006d0
      dtz1(11)    =  0.0014d0
      rx1(11)     =  0.0d0     / rhosec
      ry1(11)     =  0.0d0     / rhosec
      rz1(11)     =  0.0d0     / rhosec   !Previous  0.00006000d0 (prior to v3.4.0)
      drx1(11)    =  0.0d0     / rhosec
      dry1(11)    =  0.0d0     / rhosec
      drz1(11)    =  0.00002d0 / rhosec
      scale1(11)  = -1.55d-9              !Previous -1.58000d-9 (prior to v3.4.0)
      dscale1(11) = -0.01d-9
      refepc1(11) =  1997.0d0             !Previous 2000.0d0 (prior to v3.4.0)

*** From ITRF94 to PACP00
*** use PA/ITRF2000 rotation rates from Beavan et al. (2002)
      tx1(12)     =  0.9035d0
      ty1(12)     = -2.0202d0
      tz1(12)     = -0.5417d0
      dtx1(12)    =  0.0d0
      dty1(12)    =  0.00060d0
      dtz1(12)    =  0.00140d0
      rx1(12)     =  0.027741d0 / rhosec
      ry1(12)     =  0.013469d0 / rhosec
      rz1(12)     =  0.002712d0 / rhosec
      drx1(12)    = -0.000384d0 / rhosec
      dry1(12)    =  0.001007d0 / rhosec
      drz1(12)    = -0.002166d0 / rhosec
      scale1(12)  = -1.55d-9
      dscale1(12) = -0.01d-9
      refepc1(12) =  1997.0d0

*** From ITRF94 to MARP00
*** Use velocity of GUAM
      tx1(13)     =  0.9035d0
      ty1(13)     = -2.0202d0
      tz1(13)     = -0.5417d0
      dtx1(13)    =  0.0d0
      dty1(13)    =  0.00060d0
      dtz1(13)    =  0.00140d0
      rx1(13)     =  0.028971d0 / rhosec
      ry1(13)     =  0.010420d0 / rhosec
      rz1(13)     =  0.008928d0 / rhosec
      drx1(13)    = -0.000020d0 / rhosec
      dry1(13)    =  0.000105d0 / rhosec
      drz1(13)    = -0.000327d0 / rhosec
      scale1(13)  = -1.55d-9
      dscale1(13) = -0.01d-9
      refepc1(13) =  1997.00d0

*** From ITRF94 to ITRF2005 (also IGS05)
*** assumes that         ITRF94 = ITRF96
*** uses IERS convention ITRF96 = ITRF97
*** uses IERS values for ITRF97   -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
*** Changed epoch to 1997.0 to make consistent with trapa epoch (v3.4.0)
      tx1(14)     = -0.0074d0             !Previous -0.00680d0 (prior to v3.4.0)
      ty1(14)     = -0.0050d0             !Previous -0.00350d0 (prior to v3.4.0)
      tz1(14)     =  0.0189d0             !Previous  0.02850d0 (prior to v3.4.0)
      dtx1(14)    =  0.0002d0
      dty1(14)    =  0.0005d0
      dtz1(14)    =  0.0032d0
      rx1(14)     =  0.0d0     / rhosec
      ry1(14)     =  0.0d0     / rhosec
      rz1(14)     =  0.0d0     / rhosec   !Previous  0.00006000d0 (prior to v3.4.0)
      drx1(14)    =  0.0d0     / rhosec
      dry1(14)    =  0.0d0     / rhosec
      drz1(14)    =  0.00002d0 / rhosec
      scale1(14)  = -1.71d-9              !Previous -1.98000d-9 (prior to v3.4.0)
      dscale1(14) = -0.09d-9
      refepc1(14) =  1997.0d0             !Previous 2000.0d0 (prior to v3.4.0)

*** From ITRF94 to ITRF2008 (also IGS08 and IGb08)
*** assumes that ITRF94 = ITRF96
*** uses IERS convention ITRF96 = ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
*** Changed epoch to 1997.0 to make consistent with trapa epoch (v3.4.0)
      tx1(15)     = -0.0045d0             !Previous -0.00480d0 (prior to v3.4.0)
      ty1(15)     = -0.0041d0             !Previous -0.00260d0 (prior to v3.4.0)
      tz1(15)     =  0.0236d0             !Previous  0.03320d0 (prior to v3.4.0)
      dtx1(15)    = -0.0001d0
      dty1(15)    =  0.0005d0
      dtz1(15)    =  0.0032d0
      rx1(15)     =  0.0d0     / rhosec
      ry1(15)     =  0.0d0     / rhosec
      rz1(15)     =  0.0d0     / rhosec   !Previous  0.00006000d0 (prior to v3.4.0)
      drx1(15)    =  0.0d0     / rhosec
      dry1(15)    =  0.0d0     / rhosec
      drz1(15)    =  0.00002d0 / rhosec
      scale1(15)  = -2.65d-9              !Previous -2.92d-9 (prior to v3.4.0)
      dscale1(15) = -0.09d-9
      refepc1(15) =  1997.0d0             !Previous 2000.0d0 (prior to v3.4.0)

*** From ITRF94 to ITRF2014 (also IGS14 and IGb14)
*** assumes that ITRF94 = ITRF96
*** uses IERS convention ITRF96 = ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
*** uses IERS values for ITRF2008 -> ITRF2014
      tx1(16)     = -0.0074d0
      ty1(16)     =  0.0005d0
      tz1(16)     =  0.0628d0
      dtx1(16)    = -0.0001d0
      dty1(16)    =  0.0005d0
      dtz1(16)    =  0.0033d0
      rx1(16)     =  0.0d0     / rhosec
      ry1(16)     =  0.0d0     / rhosec
      rz1(16)     =  0.00026d0 / rhosec
      drx1(16)    =  0.0d0     / rhosec
      dry1(16)    =  0.0d0     / rhosec
      drz1(16)    =  0.00002d0 / rhosec
      scale1(16)  = -3.80d-9
      dscale1(16) = -0.12d-9
      refepc1(16) =  2010.0d0

*** From ITRF94 to ITRF2020 (also IGS20)
*** assumes that ITRF94 = ITRF96
*** uses IERS convention ITRF96 = ITRF97
*** uses IERS values for ITRF97 -> ITRF2000
*** uses IERS values for ITRF2000 -> ITRF2005
*** uses IERS values for ITRF2005 -> ITRF2008
*** uses IERS values for ITRF2008 -> ITRF2014
*** uses IERS values for ITRF2014 -> ITRF2020
      tx1(17)     = -0.0060d0
      ty1(17)     =  0.0009d0
      tz1(17)     =  0.0624d0
      dtx1(17)    = -0.0001d0
      dty1(17)    =  0.0006d0
      dtz1(17)    =  0.0031d0
      rx1(17)     =  0.0d0     / rhosec
      ry1(17)     =  0.0d0     / rhosec
      rz1(17)     =  0.00026d0 / rhosec
      drx1(17)    =  0.0d0     / rhosec
      dry1(17)    =  0.0d0     / rhosec
      drz1(17)    =  0.00002d0 / rhosec
      scale1(17)  = -3.38d-9
      dscale1(17) = -0.12d-9
      refepc1(17) =  2010.0d0

*** From ITRF94 to WGS84(G1150), added in v3.5.0.
*** Uses IERS convention ITRF96 = ITRF97.
*** Treated as biased with respect to ITRF2000.
*** Based on NGA-published transformation between WGS84(G1150)
*** and WGS84(G1762) = ITRF2008 in Table 2.5 of NGA.STND.0036_1.0.0_WGS84,
*** "Department of Defense World Geodetic System 1984: its definition and 
*** relationships with local geodetic systems", v1.0.0, 2014.  See also
*** Kelly and Dennis (2022) "Transforming Between WGS84 Realizations,"
*** J Surv Eng, doi: 10.1061/(ASCE)SU.1943-5428.0000389.
      tx1(18)     =  0.0011d0
      ty1(18)     = -0.0017d0
      tz1(18)     =  0.0290d0
      dtx1(18)    =  0.0d0
      dty1(18)    =  0.0006d0
      dtz1(18)    =  0.0014d0
      rx1(18)     =  0.0d0     / rhosec
      ry1(18)     =  0.0d0     / rhosec
      rz1(18)     =  0.00026d0 / rhosec
      drx1(18)    =  0.0d0     / rhosec
      dry1(18)    =  0.0d0     / rhosec
      drz1(18)    =  0.00002d0 / rhosec
      scale1(18)  =  1.40d-9
      dscale1(18) = -0.01d-9
      refepc1(18) =  2010.0d0

*** From ITRF94 to WGS84(G1674), added in v3.5.0.
*** Uses IERS convention ITRF96 = ITRF97.
*** Treated as biased with respect to ITRF2008.
*** Based on NGA-published bias transformation between WGS84(G1674)
*** and WGS84(G1762) = ITRF2008 in Table 2.5 of NGA.STND.0036_1.0.0_WGS84,
*** "Department of Defense World Geodetic System 1984: its definition and 
*** relationships with local geodetic systems", v1.0.0, 2014.
*** See also Kelly and Dennis (2022) "Transforming Between WGS84 Realizations,"
*** J Surv Eng, doi: 10.1061/(ASCE)SU.1943-5428.0000389.
      tx1(19)     = -0.0018d0
      ty1(19)     = -0.0006d0
      tz1(19)     =  0.0612d0
      dtx1(19)    = -0.0001d0
      dty1(19)    =  0.0005d0
      dtz1(19)    =  0.0032d0
      rx1(19)     = -0.00027d0 / rhosec
      ry1(19)     =  0.00027d0 / rhosec
      rz1(19)     = -0.00012d0 / rhosec
      drx1(19)    =  0.0d0     / rhosec
      dry1(19)    =  0.0d0     / rhosec
      drz1(19)    =  0.00002d0 / rhosec
      scale1(19)  =  3.08d-9
      dscale1(19) = -0.09d-9
      refepc1(19) =  2010.0d0

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

*** (x1, y1, z1) --> input ITRF94 coordinates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

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

*** (x1, y1, z1) --> input ITRF94 coordinates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

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

*** (x1, y1, z1) --> input coordinates (meters)
*** (x2, y2, z2) --> output  ITRF94 coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of input reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

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
      parameter (numref = 19)

      common /tranpa1/ tx1(numref), ty1(numref), tz1(numref), 
     &                dtx1(numref), dty1(numref), dtz1(numref),
     &                rx1(numref), ry1(numref), rz1(numref), 
     &                drx1(numref), dry1(numref), drz1(numref),
     &                scale1(numref), dscale1(numref), refepc1(numref)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

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
      nframe(4) = 'WGS84 original (Transit)'
      
      iframe(5) = 5
      nframe(5) = 'WGS84(G730)             '
      
      iframe(6) = 8
      nframe(6) = 'WGS84(G873)             '
      
      iframe(7) = 18
      nframe(7) = 'WGS84(G1150)            '
      
      iframe(8) = 19
      nframe(8) = 'WGS84(G1674)            '
      
      iframe(9) = 15
      nframe(9) = 'WGS84(G1762)            '
      
      iframe(10)= 16
      nframe(10)= 'WGS84(G2139)            '
      
      
c      iframe(11)= 5                          !Included with ITRF91 in v3.4.0
c      nframe(11)= 'SIO/MIT_92              ' !Included with ITRF91 in v3.4.0

      iframe(11)= 2
      nframe(11)= 'ITRF88                  '
      
      iframe(12)= 3
      nframe(12)= 'ITRF89                  '
      
      iframe(13)= 4
      nframe(13)= 'ITRF90                  '
      
      iframe(14)= 5
      nframe(14)= 'ITRF91                  '
      
      iframe(15)= 6
      nframe(15)= 'ITRF92                  '
      
      iframe(16)= 7
      nframe(16)= 'ITRF93                  '
      
      iframe(17)= 8
      nframe(17)= 'ITRF94                  '
      
      iframe(18)= 8
      nframe(18)= 'ITRF96                  '
      
      iframe(19)= 9
      nframe(19)= 'ITRF97                  '
      
      iframe(20)= 11
      nframe(20)= 'ITRF2000 or IGS00/IGb00 '
      
      iframe(21)= 14
      nframe(21)= 'ITRF2005 or IGS05       '
      
      iframe(22)= 15
      nframe(22)= 'ITRF2008 or IGS08/IGb08 '
      
      iframe(23)= 16
      nframe(23)= 'ITRF2014 or IGS14/IGb14 '
      
      iframe(24)= 17
      nframe(24)= 'ITRF2020 or IGS20       '

      write(luout, 100)  
  100 format(
     1'  1...NAD_83(2011/CORS96/2007)  North America plate fixed     '/
     1'  2...NAD_83(PA11/PACP00)       Pacific plate fixed           '/
     1'  3...NAD_83(MA11/MARP00)       Mariana plate fixed           '/
     1'                                                              '/
     1'  4...WGS84 original (Transit)                                '/
     1'  5...WGS84(G730)   ITRF91 used                               '/
     1'  6...WGS84(G873)   ITRF94=ITRF96=ITRF97 used                 '/
     1'  7...WGS84(G1150)  Biased with respect to ITRF2000           '/
     1'  8...WGS84(G1674)  Biased with respect to ITRF2008           '/
     1'  9...WGS84(G1762)  ITRF2008=IGS08=IGb08 used                 '/
     1' 10...WGS84(G2139)  ITRF2014=IGS14=IGb14 used                 '/
     1'                                                              '/
     1' 11...ITRF88                      18...ITRF96 (=ITRF94=ITRF97)'/
     1' 12...ITRF89                      19...ITRF97 (=ITRF94=ITRF96)'/
     1' 13...ITRF90 (or PNEOS90/NEOS90)  20...ITRF2000 or IGS00/IGb00'/
     1' 14...ITRF91 (or SIO/MIT_92)      21...ITRF2005 or IGS05      '/
     1' 15...ITRF92                      22...ITRF2008 or IGS08/IGb08'/
     1' 16...ITRF93                      23...ITRF2014 or IGS14/IGb14'/
     1' 17...ITRF94 (=ITRF96=ITRF97)     24...ITRF2020 or IGS20      '/)

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
      parameter (numref = 19)
      parameter (nbbdim = 10000)
      CHARACTER    OLDBB*80, NEWBB*80, NAMEIF*80
      CHARACTER    NAME24*24
      CHARACTER    OPT*1, ANSWER*1, BBTYPE*1, VOPT*1
      CHARACTER    LATDIR*1, LONDIR*1, LATDR*1, LONDR*1
      character    frame1*24, frame2*24
      character    HTDP_version*8
      LOGICAL      TEST
      LOGICAL      Is_iopt_NAD83
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      WRITE(LUOUT,20)
   20 FORMAT(' ********************************************'/
     1   ' Please enter the time to which the updated'/     
     1   ' positions and/or observations are to correspond:')
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
     1   ' Specify the reference frame of the input positions:'/)
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
     5   '     2...Update positions for Bluebook stations.'/
     6   '     3...Update values for Bluebook observations.'/
     7   '     4...Update both the positions for Bluebook stations'/
     9   '         and the values for Bluebook observations.'/
     9   '     5...Update positions for multiple points contained'/
     9   '         in a file in the format LAT,LON,EHT,TEXT:'/
     9   '         LAT = latitude in degrees (positive north)'/
     9   '         LON = longitude in degrees (positive west)'/
     9   '         EHT = ellipsoid height in meters'/
     9   '         TEXT = descriptive text (maximum 24 characters)'/
     9   '         Example: 40.731671553,112.212671753,34.241,SALT'/)
      READ(LUIN,'(A1)',err=501,iostat=ios) OPT
      if (ios /= 0) goto 501
      IF(OPT .eq. '0') THEN
         RETURN
      ELSEIF(OPT .eq. '1') THEN
         WRITE(LUOUT,1020)
 1020    FORMAT(' ************************************'/
     1       ' Enter the time',
     2       ' to which the input positions will correspond:  ')
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
         READ(LUIN,'(A80)',err=502,iostat=ios) NEWBB
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

*** Updating a Bluebook 

      ELSEIF(OPT .eq. '2' .or. OPT .eq. '3' .or. OPT .eq. '4') THEN

          if (nbbdim .eq. 10000) then
   90        WRITE(LUOUT,91)
   91        FORMAT(
     1       ' Identify type of Bluebook:'/
     2       '    1...Standard (4-digit SSN)'/
     3       '    2...Non-standard (5-digit SSN)  ')
             READ(LUIN,'(A1)',err=503,iostat=ios) BBTYPE
             if (ios /= 0) goto 503
             IF(BBTYPE .NE. '1' .AND. BBTYPE .NE. '2') GO TO 90
          else
             BBTYPE = '1'
          endif

          WRITE(LUOUT,100)
  100     FORMAT(' Enter name of the Bluebook file to be updated.'/)
          READ(LUIN,110,err=502,iostat=ios) OLDBB
          if (ios /= 0) goto 502
  110     FORMAT(A80)
          WRITE(LUOUT,120)
  120     FORMAT(' Enter name for the new Bluebook file that is to '/
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
     1           ' to which the input positions correspond:'/)
  123        CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
             IF(TEST) then
              write(luout,*) ' Do you wish to re-enter the time? (y/n)'
              read(luin, '(A1)',err=500,iostat=ios) ANSWER
              if (ios /= 0) goto 500
              IF (ANSWER .eq. 'y' .or. ANSWER .eq. 'Y') GO TO 123
              RETURN
             ENDIF
          ENDIF

*** Retrieve geodetic positions from old Bluebook file

         IF(BBTYPE .EQ. '1') THEN
              CALL GETPO4(IOPT, DATE1)
         ELSEIF(BBTYPE .EQ. '2') THEN
              CALL GETPO5(IOPT, DATE1)
         ENDIF

*** Create new Bluebook file

         OPEN(I2,FILE = NEWBB, STATUS = 'UNKNOWN')

         write (i2,127) HTDP_version
  127    format (' ***CAUTION: This file was processed using HTDP',
     &           ' version ', a8, '***')
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
           CALL UPBB4(MIN1,MIN2,OPT)
         ELSEIF(BBTYPE .EQ. '2') THEN
           CALL UPBB5(MIN1,MIN2,OPT)
         ENDIF
         CLOSE(I1, STATUS = 'KEEP')
         CLOSE(I2, STATUS = 'KEEP')

*** Update G-file

         IF(OPT .eq. '3' .or. OPT .eq. '4') THEN
  605      WRITE(LUOUT,610)
  610      FORMAT(/' Is there a G-file to be updated? (y/n)  ')
           READ(LUIN,'(A1)',err=500,iostat=ios) ANSWER
           if (ios /= 0) goto 500
           IF(ANSWER .eq. 'N' .or. ANSWER .eq. 'n') THEN
             CONTINUE
           ELSEIF(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y')THEN
             WRITE(LUOUT,620)
  620        FORMAT(' Enter name of input G-file to be updated.  ')
             READ(LUIN,'(A80)',err=502,iostat=ios) OLDBB
             if (ios /= 0) goto 502
             WRITE(LUOUT,630)
  630        FORMAT(' Enter name for the updated output G-file.  ')
             READ(LUIN,'(A80)',err=502,iostat=ios) NEWBB
             if (ios /= 0) goto 502
             OPEN(I1,FILE=OLDBB,STATUS='OLD')
             OPEN(I2,FILE=NEWBB,STATUS='UNKNOWN')
             WRITE(I2,130) MONTH2, IDAY2, IYEAR2, DATE2
  130        FORMAT(' ***CAUTION: Observations in this file have been ',
     *       'updated to ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
        
  634        WRITE(LUOUT, 632)
  632        FORMAT(/' ***************************'/
     *              ' Specify the reference frame for the updated'
     *              ' G-file vectors:'//
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
     1                ,' HTDP version ', a8, ' ***')
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
         read(luin, '(A80)',err=502,iostat=ios)NEWBB
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

*** Retrieve geodetic coordinates from the Bluebook

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

*** Retrieve geodetic position from Bluebook with 5-digit SSN

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
      SUBROUTINE UPBB4(MIN1,MIN2,OPT)

*** Update Bluebook

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
     1            ' will need to divided into two or more Bluebooks'/
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
  240        FORMAT(' Bluebook file has incorrect structure.'/
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

*** Unrecognized Bluebook record
*** Associated ELSE statement with label 500 deleted in v3.3.0.

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
      SUBROUTINE UPBB5(MIN1,MIN2,OPT)

*** Update 5-digit Bluebook

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
             ISEC1 = IDINT(SEC)           !Convert SEC to integer (v3.3.0)
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
  240           FORMAT(' Bluebook file has incorrect structure.'/
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
           ISEC1 = IDINT(SEC)             !Convert SEC to integer (v3.3.0)
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

*** Unrecognized Bluebook record

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
C*******************************************************************
      SUBROUTINE UPGFI4(DATE2, MIN2, IOPT, KOPT, MONTH, IDAY, IYEAR)

*** Update G-file of Bluebook

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER     CARD*120
      CHARACTER     TYPE*1
      CHARACTER     NRF*2
      CHARACTER     ZT*2
      CHARACTER     CHAR14*14
      LOGICAL       TEST
      LOGICAL Is_inp_NAD83, Is_out_NAD83, Is_out1_NAD83
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      ZT = 'ZT'
      
C  Initialize NAD 83 input/poutput logical variables as false:
      Is_inp_NAD83 = .FALSE.
      Is_out_NAD83 = .FALSE.
      Is_out1_NAD83 = .FALSE.

*** Obtain Bluebook reference frame identifier
*** corresponding to KOPT  (output frame for vectors)
      IF (KOPT .NE. -1) THEN
         CALL RFCON1(KOPT, NBBREF)
         WRITE(NRF,10) NBBREF
   10    FORMAT(I2.2)
      ENDIF

   90 READ(I1,100,END=200,err=300,iostat=ios) CARD
      if (ios /= 0) goto 300
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
         CALL RFCON(IBBREF, JREF)      !JREF is the current frame or frame of input
         IF (KOPT .NE. -1) THEN
            CARD(52:53) = NRF
         ENDIF
         Is_inp_NAD83 = (JREF == 1)
         Is_out_NAD83 = (IOPT == 1)    !IOPT is the frame of the positions
         Is_out1_NAD83 = (KOPT == 1)   !KOPT is the frame of the vectors (comment corrected v3.4.0)
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

*** Update G-file of 5-digit Bluebook

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

*** Obtain Bluebook reference frame identifier
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
C        CALL TOTIME(IYEAR1, 1, 1, MIN00)
C        DECYR1 = DBLE(IYEAR1) + DBLE(MINO1 - MIN00)/525600.D0
C        CALL TOTIME(IYEAR2,MONTH2,IDAY2,MINO2)
C        CALL TOTIME(IYEAR2, 1, 1, MIN00)
         CALL IYMDMJ(IYEAR1,MONTH1,IDAY1,MJD1)
         MINO1 = MJD1 * 24 * 60
         CALL IYMDMJ(IYEAR1, 1, 1, MJD0)
         DECYR1 = DBLE(IYEAR1) + DBLE(MJD1 - MJD0)/365.D0
         CALL IYMDMJ(IYEAR2, MONTH2,IDAY2, MJD2)
         MINO2 = MJD2 * 24 * 60
         CALL IYMDMJ(IYEAR2, 1, 1, MJD0)
         DECYR2 = DBLE(IYEAR2) + DBLE(MJD2 - MJD0)/365.D0
C        DECYR2 = DBLE(IYEAR2) + DBLE(MINO2 - MIN00)/525600.D0
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
*** a reference frame identifier in the Bluebook
*** to a reference frame identifier in HTDP and back

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 19 )
      COMMON /REFCON/ IRFCON(40), JRFCON(numref)

*** From Bluebook identifier to HTDP indentifier
*** WGS72 Precise
c     IRFCON(1) = 10
      IRFCON(1) = 1
C HTDP no longer supports WGS72. Hence, if a BlueBook
C file contains WGS72 coordinates, HTDP treats these
C coordinates as if they were NAD 83(2011) coordinates.

*** WGS84 (orig) Precise
      IRFCON(2) = 10   !Changed in v3.4.0 (no longer equal to NAD 83)

*** WGS72 Broadcast
c     IRFCON(3) = 10
      IRFCON(3) = 1
C HTDP no longer supports WGS72. Hence, if a BlueBook
C file contains WGS72 coordinates, HTDP treats these
C coordinates as if they were NAD 83(2011)coordinates.

*** WGS84 (orig) Broadcast
      IRFCON(4) = 10   !Changed in v3.4.0 (no longer equal to NAD 83)

*** ITRF89
      IRFCON(5) = 3

*** PNEOS 90 or NEOS 91.25 (set equal to ITRF90)
      IRFCON(6) = 4

*** NEOS 90 (set equal to ITRF90)
      IRFCON(7) = 4

*** ITRF91
      IRFCON(8) = 5

*** SIO/MIT 92.57 (set equal to ITRF91)   !This ooption no longer available as of v3.4.0
      IRFCON(9) = 5

*** ITRF91
      IRFCON(10) = 5

*** ITRF92
      IRFCON(11) = 6

*** ITRF93
      IRFCON(12) = 7

*** WGS84 (G730) Precise (set equal to ITRF91)
      IRFCON(13) = 5

*** WGS84 (G730) Broadcast (set equal to ITRF91)
      IRFCON(14) = 5

*** ITRF94
      IRFCON(15) = 8

*** WGS84 (G873) Precise  (set equal to ITRF94)
      IRFCON(16) = 8

*** WGS84 (G873) Broadcast (set equal to ITRF94)
      IRFCON(17) = 8

*** ITRF96
      IRFCON(18) = 8

*** ITRF97
      IRFCON(19) = 9

*** IGS97
      IRFCON(20) = 9

*** ITRF2000
      IRFCON(21) = 11

*** IGS00
      IRFCON(22) = 11

*** WGS84 (G1150)
      IRFCON(23) = 18

*** IGb00
      IRFCON(24) = 11

*** ITRF2005
      IRFCON(25) = 14

*** IGS05
      IRFCON(26) = 14

*** IGS08
      IRFCON(27) = 15

*** IGb08
      IRFCON(28) = 15

*** ITRF2008
      IRFCON(29) = 15

*** WGS84 (G1674)
      IRFCON(30) = 19

*** WGS84 (G1762)
      IRFCON(31) = 15

*** ITRF2014
      IRFCON(32) = 16

*** IGS14
      IRFCON(33) = 16

*** NAD83 (2011/2007/CORS96/FBN/HARN)
      IRFCON(34) = 1

*** NAD83 (PA11)
      IRFCON(35) = 12

*** NAD83 (MA11)
      IRFCON(36) = 13

*** IGb14
      IRFCON(37) = 16

*** WGS84 (G2139)
      IRFCON(38) = 16

*** ITRF2020
      IRFCON(39) = 17

*** IGS20
      IRFCON(40) = 17

*** From HTDP identifier to Bluebook identifier.
*** NAD 83 (2011/2007/CORS96/...) referenced to North America plate.
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

*** WGS72
      JRFCON(10) = 1

*** ITRF2000
      JRFCON(11) = 21

*** NAD 83 (PACP00) or NAD 83 (PA11)
      JRFCON(12) = 35

*** NAD 83 (MARP00) or NAD 83 (MA11)
      JRFCON(13) = 36

*** ITRF2005 or IGS05
      JRFCON(14) = 26

*** ITRF2008 or IGS08
      JRFCON(15) = 27

*** IGb08
      JRFCON(15) = 28

*** ITRF2014 or IGS14
      JRFCON(16) = 33

*** IGb14
      JRFCON(16) = 37

*** ITRF2020 or IGS20
      JRFCON(17) = 40

      RETURN
      END
***************************************************
      SUBROUTINE RFCON(IBBREF, JREF)

*** Convert reference frame identifier from
*** system used in the Bluebook to the
*** system used in HTDP

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 19 )
      COMMON /REFCON/ IRFCON(40), JRFCON(numref)

      IF (1 .LE. IBBREF .AND. IBBREF .LE. 40) THEN
          JREF = IRFCON(IBBREF)
      ELSE
          WRITE(6, 10) IBBREF
   10     FORMAT(' Improper reference frame identifier (=',
     1      I4, ')' /
     1      ' appearing in B-record of the G-file')
          STOP
       ENDIF

       RETURN
       END
*******************************************************
      SUBROUTINE RFCON1(JREF, IBBREF)

*** Convert reference frame identifier from
*** system used in HTDP to the system
*** used in the Bluebook

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 19 )
      COMMON /REFCON/ IRFCON(40), JRFCON(numref)

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
      parameter (numref = 19)

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
      parameter (numref = 19)

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
  100 FORMAT(' Enter name for point (maximum 24 characters).')
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
     2    ' with north positive, for example 35 17 28.3 or 35,17,28.3 '/
     2    ' (for the southern hemisphere, enter a minus sign before'/
     3    ' each value, for example -35 -17 -28.3 or -35,-17,-28.3)')
     
        READ(LUIN,*,err=202,iostat=ios) LATD, LATM, SLAT
        if (ios /= 0) goto 202

*** START insertion of new code and comments for v3.3.0 to handle 
*** incorrect input latitude and longitude.
***
*** Test whether latitude magnitude is > 90 degrees. If it is, prompt 
*** for re-entry.
        TEST = DBLE(ABS(LATD)) + DBLE(ABS(LATM))/60.D0 
     &         + DABS(SLAT)/3600.D0
        DO WHILE (DABS(TEST) .GT. 90.D0)
          WRITE(*,*)'  Latitude magnitude (degrees) =', TEST
          WRITE(*,*)
     &      '  Magnitude cannot exceed 90. Please re-enter all values.'
          READ(LUIN,*,err=202,iostat=ios) LATD, LATM, SLAT
          if (ios /= 0) goto 202
          TEST = DBLE(ABS(LATD)) + DBLE(ABS(LATM))/60.D0 
     &           + DABS(SLAT)/3600.D0
        END DO

*** Non-zero latitude minutes and seconds must be the same sign, and
*** if latitude degrees is non-zero, minutes and seconds must be the 
*** same sign as latitude degrees. Otherwise prompt for re-entry.
        IF (LATD .EQ. 0) THEN
          IF (DBLE(LATM)*SLAT .LT. 0.D0) THEN
            WRITE(*,*)'  Minutes and seconds must be the same sign.',
     &        ' Please re-enter minutes and seconds.'
            READ(LUIN,*,err=202,iostat=ios) LATM, SLAT
            if (ios /= 0) goto 202
          ENDIF
        ENDIF
        IF (LATD .GT. 0) THEN
          IF ((LATM .LT. 0).OR.(SLAT .LT. 0.D0))
     &       WRITE(*,*)'  Latitude degrees are positive.'
          DO WHILE (LATM .LT. 0)
              WRITE(*,*)
     &     '  Minutes must also be positive. Please re-enter minutes.'
            READ(LUIN,*,err=202,iostat=ios) LATM
            if (ios /= 0) goto 202
          END DO
          DO WHILE (SLAT .LT. 0.D0)
            WRITE(*,*)
     &     '  Seconds must also be positive. Please re-enter seconds.'
            READ(LUIN,*,err=202,iostat=ios) SLAT
            if (ios /= 0) goto 202
          END DO
        ENDIF

        IF (LATD .LT. 0) THEN
          IF ((LATM .GT. 0).OR.(SLAT .GT. 0.D0))
     &       WRITE(*,*)'  Latitude degrees are negative.'
          DO WHILE (LATM .GT. 0)
            WRITE(*,*)
     &     '  Minutes must also be negative. Please re-enter minutes.'
            READ(LUIN,*,err=202,iostat=ios) LATM
            if (ios /= 0) goto 202
          END DO
          DO WHILE (SLAT .GT. 0.D0)
            WRITE(*,*)
     &     '  Seconds must also be negative. Please re-enter seconds.'
            READ(LUIN,*,err=202,iostat=ios) SLAT
            if (ios /= 0) goto 202
          END DO
        ENDIF

*** Test whether latitude minute or second magnitudes are >= 60. 
*** If either is, prompt for re-entry.
        DO WHILE (ABS(LATM) .GE. 60)
          WRITE(*,*)'  Minutes magnitude must be < 60.',
     &              ' Please re-enter latitude minutes.'
          READ(LUIN,*,err=202,iostat=ios) LATM
          if (ios /= 0) goto 202
        END DO

        DO WHILE (DABS(SLAT) .GE. 60.D0)
          WRITE(*,*)'  Seconds magnitude must be < 60.',
     &              ' Please re-enter latitude seconds.'
          READ(LUIN,*,err=202,iostat=ios) SLAT
          if (ios /= 0) goto 202
        END DO
        
        WRITE(LUOUT,120)
  120   FORMAT(
     1    ' Enter longitude degrees-minutes-seconds in free format'/,
     2    ' with west being positive.  To express a longitude measured'/
     3    ' eastward, enter a minus sign before each value.')
        READ(LUIN,*,err=203,iostat=ios) LOND, LONM, SLON
          if (ios /= 0) goto 203
          
*** Test whether longitude magnitude is > 360 degrees. If it is, prompt 
*** for re-entry.
        TEST = DBLE(ABS(LOND)) + DBLE(ABS(LONM))/60.D0 
     &         + DABS(SLON)/3600.D0
        DO WHILE (DABS(TEST) .GT. 360.D0)
          WRITE(*,*)'  Longitude magnitude (degrees) =', TEST
          WRITE(*,*)
     &      '  Magnitude cannot exceed 360. Please re-enter all values.'
          READ(LUIN,*,err=203,iostat=ios) LOND, LONM, SLON
          if (ios /= 0) goto 203
          TEST = DBLE(ABS(LOND)) + DBLE(ABS(LONM))/60.D0 
     &           + DABS(SLON)/3600.D0
        END DO

*** Non-zero longitude minutes and seconds must be the same sign, and
*** if longitude degrees is non-zero, minutes and seconds must be the 
*** same sign as longitude degrees. Otherwise prompt for re-entry.
        IF (LOND .EQ. 0) THEN
          IF (DBLE(LONM)*SLON .LT. 0.D0) THEN
            WRITE(*,*)'  Minutes and seconds must be the same sign.',
     &        ' Please re-enter minutes and seconds.'
            READ(LUIN,*,err=203,iostat=ios) LONM, SLON
            if (ios /= 0) goto 203
          ENDIF
        ENDIF

        IF (LOND .GT. 0) THEN
          IF ((LONM .LT. 0).OR.(SLON .LT. 0.D0))
     &      WRITE(*,*)'  Longitude degrees are negative.'
          DO WHILE (LONM .LT. 0)
            WRITE(*,*)
     &     '  Minutes must also be positive. Please re-enter minutes.'
            READ(LUIN,*,err=203,iostat=ios) LONM
            if (ios /= 0) goto 203
          END DO
          DO WHILE (SLON .LT. 0.D0)
            WRITE(*,*)
     &     '  Seconds must also be positive. Please re-enter seconds.'
            READ(LUIN,*,err=203,iostat=ios) SLON
            if (ios /= 0) goto 203
          END DO
        ENDIF

        IF (LOND .LT. 0) THEN
          IF ((LONM .GT. 0).OR.(SLON .GT. 0.D0))
     &       WRITE(*,*)'  Longitude degrees are negative.'
          DO WHILE (LONM .GT. 0)
            WRITE(*,*)
     &     '  Minutes must also be negative. Please re-enter minutes.'
            READ(LUIN,*,err=203,iostat=ios) LONM
            if (ios /= 0) goto 203
          END DO
          DO WHILE (SLON .GT. 0.D0)
            WRITE(*,*)
     &     '  Seconds must also be negative. Please re-enter seconds.'
            READ(LUIN,*,err=203,iostat=ios) SLON
            if (ios /= 0) goto 203
          END DO
        ENDIF
          
*** Test whether longitude minute or second magnitudes are >= 60. 
*** If either is, prompt for re-entry.
        DO WHILE (ABS(LONM) .GE. 60)
          WRITE(*,*)'  Minutes magnitude must be < 60.',
     &              ' Please re-enter longitude minutes.'
          READ(LUIN,*,err=203,iostat=ios) LONM
          if (ios /= 0) goto 203
        END DO

        DO WHILE (DABS(SLON) .GE. 60.D0)
          WRITE(*,*)'  Seconds magnitude must be < 60.',
     &              ' Please re-enter longitude seconds.'
          READ(LUIN,*,err=203,iostat=ios) SLON
          if (ios /= 0) goto 203
        END DO

*** END insertion of new code and comments for v3.3.0 to handle 
*** incorrect input latitude and longitude.

        WRITE(LUOUT, 125)
  125   FORMAT(
     1  ' Enter ellipsoid height in meters. (Note that'/,
     1  ' predicted motions are independent of this height.)  ')
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
     1       '     2...user will specify the global X-Y-Z'/
     1       '         components of site velocity.' )
      else
        write(luout, 211)
  211   FORMAT(' How do you wish to specify the velocity: '/
     1       '     1...north-east-up components.'/
     1       '     2...global X-Y-Z components.' )
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
*** Input VN, VE, and VU are in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NDLOC = 2740)

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
       
** Compute the contribution due to earthquakes.
** It is assumed that the components of displacement,
** (DNORTH, DWEST, DUP) do not vary from one reference
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

      CALL RADII(YLAT,RADMER,RADPAR)

      DN =  RADMER * (YLAT2 - YLAT1)
      DE = -RADPAR * (YLON2 - YLON1)
      DU =  HT2 - HT1

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
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

** Get reference latitude (RLAT) and reference longitude (RLON)
C  The following 2 lines were added on 07/22/2015 after Rich found this bug
         elon = -ylon
         call TOXYZ(ylat, elon, eht, x, y, z)

         IF(IOPT .EQ. 15) THEN                   !Velocity grids are in ITRF2008
            RLAT = YLAT
***     Added conditional statement to get rid of out-of-region error when
***     YLON is negative for reference frame (v3.3.0):
            IF(YLON .GE. 0.D0) THEN
              RLON = YLON
            ELSE
              RLON = YLON + TWOPI
            ENDIF
         ELSE
             CALL XTO08 (X,Y,Z,RLAT,RLON,EHTNAD,DATE,IOPT)   !Positions are in ITRF2008
         ENDIF

** Get deformation region

         CALL GETREG(RLAT,RLON,JREGN)
         IF (JREGN .EQ. 0) THEN
           VN = 0.D0
           VE = 0.D0
           VU = 0.D0
           RETURN
         ENDIF
         CALL COMVEL( RLAT, RLON, JREGN, VN, VE, VU)       !Those velocities are in ITRF2008

** Convert  velocity to reference of iopt, if frame not ITRF2008

         Is_iopt_NAD83 = (iopt == 1)
         IF (IOPT .NE. 15) THEN
           CALL TOVXYZ( YLAT, ELON, VN, VE, VU, VX, VY, VZ)
           if (Is_iopt_NAD83) then
             CALL VTRANF( X, Y, Z, VX, VY, VZ, 15, IOPT)
           else
             CALL VTRANF_IERS( X, Y, Z, VX, VY, VZ, 15, IOPT)
           endif
           CALL TOVNEU( YLAT, ELON, VX, VY, VZ, VN, VE, VU)
         ENDIF

         RETURN
         END

****************************************************************
      subroutine TRFVEL

*** Transform velocities from one reference frame to another

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 19)
      character    nameif*80,name24*80
      character    NAMEF*80
      character    frame1*24, frame2*24
      character    option*1
      character    vopt*1, LATDIR*1,LONDIR*1
      character    record*120                 
      LOGICAL      Is_inp_NAD83, Is_out_NAD83
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

C  Initialize NAD 83 input/poutput logical variables as false:
      Is_inp_NAD83 = .FALSE.
      Is_out_NAD83 = .FALSE.

      write( luout, 100)
  100 format(
     1  ' Please enter the name of the file to contain '/
     1  ' the transformed velocities. ')
      read( luin, '(A80)',err=600,iostat=ios) namef
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

      IF(iopt1 .EQ. 1) Is_inp_NAD83 = .TRUE.
c      Is_inp_NAD83 = (iopt1 ==  1)    !Comment out v3.4.0

  115 write( luout, 120)
  120 format( /' Enter the reference frame for the output velocities.')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
      write( luout, *) 'Improper selection -- try again.'
      go to 115
      endif
      
      IF(iopt2 .EQ. 1) Is_out_NAD83 = .TRUE.
c      Is_out_NAD83 = (iopt2 ==  1)    !Comment out v3.4.0

      write( i2, 125) frame1, frame2
  125 format( ' TRANSFORMING VELOCITIES FROM ', A24, ' TO ', a24//
     1   16X, ' INPUT VELOCITIES      OUTPUT VELOCITIES'/)

  130 write( luout, 140)
  140 format(' ************************************************'/
     1  ' Velocities will be transformed at each specified point.'/
     1  ' Please indicate how you wish to input points.'/
     1  '    0...No more points.  Return to main menu.'/
     1  '    1...Individual points entered interactively.'/
     1  '    2...Transform velocities contained in file of delimited'/
     1  '        records of the form LAT,LON,VN,VE,VU,TEXT:' /
     1  '        LAT = latitude in degrees (positive north)'/
     1  '        LON = longitude in degrees (positive west)'/
     1  '        VN = northward velocity in mm/yr '/
     1  '        VE = eastwars velocity in mm/yr '/
     1  '        VU = upward velocity in mm/yr '/
     1  '        TEXT = descriptive text (maximum 24 characters) '/
     1  '        Example: '/
     1  '          40.731671553,112.212671753,3.7,3.8,-2.4,SALT AIR '/)
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
         LINE = 0   !Inititalize line number of input file.
  210    read(i1,'(a)',end = 220,err=603,iostat=ios) record          
         call interprate_velocity_record (record,xlat,xlon,vn,ve,vu,
     &                                    name24)
         if (ios /= 0) goto 603
         LINE = LINE + 1   !Increment line number of input file.
        
*** If latitude magnitudes > 90 degrees or longitude magnitude > 360 degrees,
*** write error message and terminate program (added for v3.3.0).
         IF ((DABS(xlat) .GT. 90).OR.(DABS(xlon) .GT. 360)) THEN
           WRITE(*,*)'***********************************************'
           WRITE(*,*)'Invalid latitude or longitude in input file'
           WRITE(*,2101)'on line', LINE, ':'
           WRITE(*,2102) xlat, xlon, vn, ve, vu, name24 
           WRITE(*,*)'Please check your input file and try again.'
           WRITE(*,*)'***********************************************'
 2101      FORMAT (1X, A, I8, A/)
 2102      FORMAT (1X, 2F17.10, 3F8.2, 2X, A)
           STOP
         ENDIF

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
      parameter (numref = 19)
      common /tranpa/ tx(numref), ty(numref), tz(numref),
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref),
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      IF(IOPT1 .le. numref .and. IOPT2. le. numref
     &   .and. IOPT1 .gt. 0 .and. IOPT2 .gt. 0 ) THEN

*** Convert from mm/yr to m/yr
         VX = VX / 1000.d0
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
      parameter (numref = 19)
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
  100 format(' ************************************************'/
     1  ' New northward velocity = ', f8.2, ' mm/yr' /
     1  ' New eastward velocity  = ', f8.2, ' mm/yr'/
     1  ' New upward velocity    = ', f8.2, ' mm/yr'/
     1  ' New x velocity         = ', f8.2, ' mm/yr'/
     1  ' New y velocity         = ', f8.2, ' mm/yr'/
     1  ' New z velocity         = ', f8.2, ' mm/yr')
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
       character  HTDP_version*8
       COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
       COMMON /VERSION/ HTDP_version                        

       WRITE(I2, 10) HTDP_version
   10  FORMAT(' HTDP OUTPUT, VERSION ',a / )
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
     1        '     1. Month-day-year in free format.'/
     1        '        For example 5 12 1979 or 5,12,1979'/
     1        '        represents May 12, 1979'/
     1        '     2. Decimal year.'/
     1        '        For example 2010.000 represents UTC midnight'/
     1        '        at the beginning of January 1, 2010.'/)

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
          write(luout,*) 'Enter decimal year '
          READ (LUIN, *,err=102,iostat=ios) DATE
          if (ios /= 0)  goto 102
          
          IF (DATE .lt. 1906.0d0) then
             write (luout, 10)
             TEST = .TRUE.
             RETURN
          ENDIF

          IYEAR = IDINT(DATE)   !Convert DATE to integer (v3.3.0)
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
              IDAY = IDINT(REMDAY) - IBEGIN + 1   !Convert REMDAY to integer (v3.3.0)
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
                  IDAY = IDINT(REMDAY) - IBEGIN + 1   !Convert REMDAY to integer (v3.3.0)
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
C       COMMENTS:              REQUIRES 4-DIGIT YEAR (E.G. 1992 NOT 92).  
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
*** Following 4 lines edited in v3.3.0 to address compiler issue with
*** real to integer conversion.
C      A=  IYRP*0.01D0
C      B=  2 - A + DINT( A*0.25D0 )
C      C=  365.25D0*IYRP
C      D=  30.6001D0*(IMOP + 1)
C      MJD =  (B + C + D + IDAY - 679006) 
C      
C      WRITE(*,*)'Old values:'
C      WRITE(*,*)'A, B, C, D =', A, B, C, D
C      WRITE(*,*)'MJD =', MJD
      
      A = IDINT(DBLE(IYRP)*0.01D0)
      B = 2 - A + IDINT(DBLE(A)*0.25D0)
      C = IDINT(365.25D0*DBLE(IYRP))
      D = IDINT(30.6001D0*(DBLE(IMOP + 1)))
      MJD = B + C + D + IDAY - 679006

C      WRITE(*,*)'New values:'
C      WRITE(*,*)'A, B, C, D =', A, B, C, D
C      WRITE(*,*)'MJD =', MJD
C      WRITE(*,*)
 
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
C PURPOSE:     RETURNS THE INDICES OF THE LOWER-LEFT
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
C PURPOSE:     RETRIEVES THE AMPLITUDES OF THE FOUR
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
C  But "trim" were not available in f77.

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

         integer*4   i,length
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

         integer*4   i,length
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

         integer*4   i,length
         real*8      x,y
         character   name*24,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80

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

         integer*4   i,length
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
C-----------------------------------------------------------------------