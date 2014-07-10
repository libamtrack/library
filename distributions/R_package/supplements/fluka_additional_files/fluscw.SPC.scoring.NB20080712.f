*$ CREATE FLUSCW.FOR
*COPY FLUSCW
*                                                                      *
*=== fluscw ===========================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION FLUSCW ( IJ    , PLA   , TXX   , TYY   ,
     &                                   TZZ   , WEE   , XX    , YY    ,
     &                                   ZZ    , NREG  , IOLREG, LLO   ,
     &                                   NSURF )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1989-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*     New version of Fluscw for FLUKA9x-FLUKA200x:                     *
*                                                                      *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*     !!! This is a completely dummy routine for Fluka9x/200x. !!!     *
*     !!! The  name has been kept the same as for older  Fluka !!!     *
*     !!! versions for back-compatibility, even though  Fluscw !!!     *
*     !!! is applied only to estimators which didn't exist be- !!!     *
*     !!! fore Fluka89.                                        !!!     *
*     !!! User  developed versions  can be used for  weighting !!!     *
*     !!! flux-like quantities at runtime                      !!!     *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*           Ij = (generalized) particle code (Paprop numbering)        *
*          Pla = particle laboratory momentum (GeV/c) (if > 0),        *
*                or kinetic energy (GeV) (if <0 )                      *
*    Txx,yy,zz = particle direction cosines                            *
*          Wee = particle weight                                       *
*     Xx,Yy,Zz = position                                              *
*         Nreg = (new) region number                                   *
*       Iolreg = (old) region number                                   *
*          Llo = particle generation                                   *
*        Nsurf = transport flag (ignore!)                              *
*                                                                      *
*     Output variables:                                                *
*                                                                      *
*       Fluscw = factor the scored amount will be multiplied by        *
*       Lsczer = logical flag, if true no amount will be scored        *
*                regardless of Fluscw                                  *
*                                                                      *
*     Useful variables (common SCOHLP):                                *
*                                                                      *
*     Flux like binnings/estimators (Fluscw):                          *
*          ISCRNG = 1 --> Boundary crossing estimator                  *
*          ISCRNG = 2 --> Track  length     binning                    *
*          ISCRNG = 3 --> Track  length     estimator                  *
*          ISCRNG = 4 --> Collision density estimator                  *
*          ISCRNG = 5 --> Yield             estimator                  *
*          JSCRNG = # of the binning/estimator                         *
*                                                                      *
*----------------------------------------------------------------------*
*   
      INCLUDE '(SCOHLP)'
      INCLUDE '(USRTRC)'
      INCLUDE '(TRACKR)'
      INCLUDE '(PAPROP)'
      INCLUDE '(FHEAVY)'
*
       LOGICAL LFIRST
       SAVE LFIRST
       DATA LFIRST /.TRUE./

* ../../flutil/ldpm3qmd fluscw.f
* return message from first call
       IF (LFIRST) THEN
          WRITE(LUNOUT,*) 'FLUSCW TEST_NB'
          LFIRST = .FALSE.
       END IF
       
       FLUSCW=ZERZER
       
**     WRITE(LUNOUT,*) '_NB ',TITUTC(JSCRNG),IJ,JSCRNG
       
*     All ions:
       IF (TITUTC(JSCRNG) .EQ. 'ALL     ') THEN
          FLUSCW=ONEONE
       END IF
       
*     all pions
       IF ((IJ .EQ. 13) .OR. (IJ .EQ. 14) .OR. (IJ .EQ. 23)) THEN
          IF ((TITUTC(JSCRNG) .EQ. 'PI*    ') .OR. (TITUTC(JSCRNG) .EQ. 
     $         'H*    ')) THEN
             FLUSCW = ONEONE
             PLA = PLA / 0.149834553
          END IF
       END IF

*     all kaons
       IF ((IJ .EQ. 15) .OR. (IJ .EQ. 16) .OR. (IJ .EQ. 24)) THEN
          IF ((TITUTC(JSCRNG) .EQ. 'K*     ') .OR. (TITUTC(JSCRNG) .EQ. 
     $         'H*    ')) THEN
             FLUSCW = ONEONE
             PLA = PLA / 0.529894834
          END IF
       END IF

       
* all p,pbar
       IF ((IJ .EQ. 2) .OR. (IJ .EQ. 1)) THEN
          IF (TITUTC(JSCRNG) .EQ. 'H*    ') THEN
             FLUSCW = ONEONE
             PLA = PLA / 1.0
          END IF
       END IF
*     d
       IF ((IJ .EQ. -3) .AND. (TITUTC(JSCRNG) .EQ. 'H*    ')) THEN
          FLUSCW = ONEONE
          PLA = PLA / 2.0
       END IF
*     t
       IF ((IJ .EQ. -4) .AND. (TITUTC(JSCRNG) .EQ. 'H*    ')) THEN
          FLUSCW = ONEONE
          PLA = PLA / 3.0
       END IF
       
* d and t seperate
      IF ((IJ .EQ. -3) .AND. (TITUTC(JSCRNG) .EQ. 'D       ')) THEN
         FLUSCW=ONEONE
         PLA = PLA / 2.0
      END IF
      IF ((IJ .EQ. -4) .AND. (TITUTC(JSCRNG) .EQ. 'T       ')) THEN
         FLUSCW=ONEONE
         PLA = PLA / 3.0
      END IF


* All Z = 2
      IF (IJ .EQ. -5) THEN 
         IF ((TITUTC(JSCRNG) .EQ. 'HE3     ') .OR. (TITUTC(JSCRNG) .EQ. 
     $        'HE*     ')) THEN
            FLUSCW=ONEONE
            PLA = PLA/3.0
         END IF
      END IF
      IF (IJ .EQ. -6) THEN 
         IF ((TITUTC(JSCRNG) .EQ. 'HE4     ') .OR. (TITUTC(JSCRNG) .EQ. 
     $       'HE*     ')) THEN
            FLUSCW=ONEONE
            PLA = PLA/4.0
         END IF
      END IF

* testing, recording NOT per nucleon
      IF (IJ .EQ. -6) THEN
         IF ((TITUTC(JSCRNG) .EQ. 'HE4NO   '))  THEN
            FLUSCW=ONEONE
*     PLA = PLA/4.0
         END IF
      END IF
      
* Heavy ions:


      IF (IJ .EQ. -2) THEN
         IA = IBARCH(IJ)
         PLA = PLA/DBLE(IA)
         CALL USRDCI(IJ,IONA,IONZ,IONM)

c$$$      IF (IJ .GT. -2) THEN
c$$$         CALL USRDCI(IJ,IONA,IONZ,IONM)
c$$$         IF (-6 .LE. IJ ) THEN
c$$$            IA = IBARCH(IJ)
c$$$         ELSE IF (IJ .LT. -6) THEN
c$$$            IA = IBHEAV(-IJ)
c$$$         END IF
c$$$         PLA = PLA/DBLE(IA)
* All Li
         IF ((IONZ .EQ. 3) .AND. (TITUTC(JSCRNG) .EQ. 'LI*     ')) THEN
            FLUSCW=ONEONE
         END IF   
         IF ((IONZ .EQ. 4) .AND. (TITUTC(JSCRNG) .EQ. 'BE*     ')) THEN
            FLUSCW=ONEONE
         END IF   
         IF ((IONZ .EQ. 5) .AND. (TITUTC(JSCRNG) .EQ. 'B*      ')) THEN
            FLUSCW=ONEONE
         END IF   
         IF ((IONZ .EQ. 6) .AND. (TITUTC(JSCRNG) .EQ. 'C*      ')) THEN
            FLUSCW=ONEONE
         END IF
         IF ((IONZ .EQ. 7) .AND. (TITUTC(JSCRNG) .EQ. 'N*      ')) THEN
            FLUSCW=ONEONE
         END IF
         IF ((IONZ .EQ. 8) .AND. (TITUTC(JSCRNG) .EQ. 'O*      ')) THEN
            FLUSCW=ONEONE
         END IF   
      END IF

      LSCZER = .FALSE.
      RETURN
*=== End of function Fluscw ===========================================*
      END
