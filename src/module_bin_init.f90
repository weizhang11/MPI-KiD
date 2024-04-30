!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module contains routine used to initialise and set-up
! the Tel-Aviv university (TAU) warm,
! size-resolved cloud microphysics scheme described in Tzivion et al,
! JAS (1987,JAS), Feingold et al (1988, JAS), and Tzivion et al (1989,
! JAS).
!
! In this version of the TAU microphysics the cloud drop size
! distribution is divided into 34 bins with a radii range of 1.56 to
! 3200 microns and mass doubling from one bin to the next. The
! method of moments (Tzivion et al. 1987, JAS) is used to solve for
! mass and number concentration in each size bin that result from
! diffusional growth (Tzivion et al 1989, JAS), collision-coalescence
! and collisional breakup (Tzivion et al, 1987 and Feingold et al,
! 1989, JAS). Sedimentation is performed with a first-order upwind scheme.
! Aerosol are represented by a single prognostic variable that is
! assumed to be ammonium sulfate with a log-normal distribution (Stevens
! et al 1996, JAS).
!
! The numerical methods and code in this module have been used in
! a variety of 2-D and 3-D dynamical frameworks to investigate a number
! of cloud microphysical problems. For example, drizzle production in marine
! Sc (Feingold et al, 1996), the dynamic and microphysical
! details of non-precipitating and precipitating marine Sc
! (Stevens et al,JAS, 1996 & 1998), the effect of drizzle on cloud optical
! depth and susceptibility (Feingold et al, JGR, 1997),
! the role of giant CCN in marine Sc,
! (Feingold et al, JAS, 1999), the role of giant CCN in cumulus (Yin et al,
! Atmospheric Research, 2000), turbulence, condensation and liquid water
! transport in non-precipitating marine Sc (Wang et al, JAS, 2003) and
! aerosol-cloud interactions (Feingold et al, GRL, 2005; Jiang et al, JGR,
! 2006; Xue and Feingold,JAS,2006; Xue et al, JAS, 2008; Hill et al, JAS,
! 2009)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module module_bin_init

  Use column_variables
  Use physconst, only : p0, this_r_on_cp=>r_on_cp, pi
  Use mphys_tau_bin_declare
  Use switches, only: l_sediment, mphys_var, l_fix_aerosols &
       , l_act, l_cond_evap, l_coll_coal, l_break, l_fix_supersat
  Use namelists, only: aero_sig_init, aero_rd_init
!  Use lem2_4_tau

contains

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
!DECK CLOUDSINIT
!***********************************************************
      SUBROUTINE bin_init
!***********************************************************
! this subroutine is called if bin resolved moicrophysics is 
! required. CLOUNDSINIT, declares the required variables, 
! defines and calculates the mass categories for usage in the bin
! model that is called in DYNVIS. The original code was IMPLICIT
! double, I've IMPLICIT noned it, and I compile double precision


      IMPLICIT NONE 

      REAL,PARAMETER :: RRV=461 !gas constant of moist air J/kg/K
      REAL,PARAMETER :: AL=597  !latent heat of evaporation and cond
      INTEGER K,L,J, ih, imom
      REAL,PARAMETER :: T00=273.15 !freezing point of water in K,from
                                   !the LES
      REAL,PARAMETER :: CPBIN=0.24 !specific heat of water vapor
      REAL,PARAMETER :: AAR=0.24E-3
      REAL,PARAMETER :: UM = 1.72E-04 !dynamic viscosity of air
      REAL ALCP!latent heat of evap divided by specific heat of water vapor
      REAL ZRU
           
!*******************************************************************
!     USING MKS UNITS: Kg,Mts,Sec,KCal
!*******************************************************************
!******************************************************************
!      DEFINING MASS CATEGORIES (XK) & AVERAGE CATEGORY (XAVE)
!      AND XI PARAMETER
!******************************************************************
      ZI=1./3.
      XI=0.5*(1.+(1.+0.5*ZI)**(1.-ZI)*(1.-0.5*(1.-ZI))**ZI)
      ZI=1./2.
      XI2=0.5*(1.+(1.+0.5*ZI)**(1.-ZI)*(1.-0.5*(1.-ZI))**ZI)
      ZI=2./3.
      XI23=0.5*(1.+(1.+0.5*ZI)**(1.-ZI)*(1.-0.5*(1.-ZI))**ZI)
      DIONE(1)=3.125E-6
      XK(1)=3.14159*DIONE(1)*DIONE(1)*DIONE(1)*1000./6.0
      XK_gr(1) = XK(1)*1.e3
      DO 18 L=2,LK
        XK(L)=2.0*XK(L-1)
        XK_gr(l) = XK(l)*1.e3
        DIONE(L)=2.0*DIONE(L-1)
 18   CONTINUE
      DO 109 L=1,LK
        XKK1(L)=((XK(L)*6./3141.59)**(1./3.))/2.
        ZRU=XKK1(L)*2.*1E6
        PRINT 108,L,XKK1(L),XK(L),ZRU
 108    FORMAT (' CATEGORY:',I4,3X,'RADIUS(m): ',E10.4,3X,'MASS: ',     &
     &  E10.4,3X,'DIAMETER(microns):',F11.5)
 109  CONTINUE
!calculate mean size of each bin
      DO L=1, LK-1
        XKmean(L) = 0.5*(XK(L)+XK(L+1))
      ENDDO
        XKmean(LK) = 0.5*(XK(LK)+(2*XK(LK)))

!******************************************************************
!    CALCULATING THE MASS CATEGORIES TO BE USED IN SUBROUTINES
!******************************************************************
      DIAM(1)=1.5625*2.E-04
      X_bin(1)=PI/6.*DIAM(1)**3
      DO 44 L=2,LKDD
        X_bin(L)=X_bin(1)*2.**(L-1)
        DIAM(L)=(6./PI*X_bin(L))**(1./3.)
44    CONTINUE
      do l=1,LKDD
         X2(l)=X_bin(l)*X_bin(l)
         X13(l)=X_bin(l)**(1./3.)
         X23(l)=X_bin(l)**(2./3.)
         X12(l)=sqrt(X_bin(l))
      enddo
!********************************************************
!        CALCULATING COEFFICIENTS FOR DROPS VELOCITY
!********************************************************
      DO K=1,LK
        IF(K >= 1.AND.K <  15)THEN
          BET(K)=2./3.
          ALP(K)=0.45795E+06
        ELSEIF(K >= 15.AND.K <  25)THEN
          BET(K)=1./3.
          ALP(K)=4.962E+03
        ELSEIF(K >= 25.AND.K <  32)THEN
          BET(K)=1./6.
          ALP(K)=1.732E+03
        ELSEIF(K >= 32)THEN
          BET(K)=0.
          ALP(K)=917.
        ENDIF
      ENDDO
!*********************************************************
! Reynolds parameters for evaporation in the bin model   
!********************************************************
       do k = 1,kkp
         SC23 = (RHON(K)/1.64)**(2./3.)
         DO L=1,LK
            RE(L,K)=(DIAM(L)*(RHON(K)*1.e-03)/1.72e-04) &
     &               *ALP(L)*X_bin(L)**BET(L)
            IF(RE(L,K).LT.2.)THEN
               VNTF(L,K)=1.00+0.108*SC23*RE(L,K)
            ELSE
               VNTF(L,K)=0.78+0.308*SQRT(SC23*RE(L,K))
            ENDIF 
          ENDDO
          DV(K)=0.211*(TREF(K)/273.15)**1.94*(101250./PREFN(1,K))
          SC(K)=UM/(RHON(K)*1.E-03*DV(K))
        ENDDO

!**********************************************************
! KiD - set the aerosol distribution
!*********************************************************
        ! assumes radius in entered in namelist in m
        ! converts to diameter in cm
        DG1 = (aero_rd_init(1)*2.)*100.
        SG1 = aero_sig_init(1)
        

      END subroutine bin_init
!
!DECK DATA
!************************************************************
       SUBROUTINE DATA
!************************************************************
      IMPLICIT NONE

      INTEGER I,J,K
      REAL XP
    OPEN(30,FILE='./src/tau_data/KIJ',STATUS='old')
    OPEN(31,FILE='./src/tau_data/KBARF',STATUS='old')
    OPEN(32,FILE='./src/tau_data/PLL',STATUS='old') 
    OPEN(33,FILE='./src/tau_data/ACON1',STATUS='old')

!      PI=4.*ATAN(1.)
    
!************** SCON COEFFICIENTS   ACON1 *******
      DO 456 I=1,18
        DO 456 J=1,I
!       READ (3,801) A(I,J)
          READ(33,707) ACON(I,J)
456   CONTINUE
707   FORMAT(2X,E12.5)
!************** SXY COEFFICIENTS  KBARF ********
      DO 849 I=1,LK-1
        DO 849 J=1,I
!       IF(J >  15)GO TO 849
          READ(31,808)KBAR(I,J)
           IF (KBAR(I,J) <  0.0)KBAR(I,J)=0.0
849   CONTINUE
808   FORMAT(2X,E12.5)
!************************
      AD(1)=0.01
      AX(1)=PI*AD(1)**3/6.
      XP=2.
      DO 1 I=2,20
 4      AX(I)=XP*AX(I-1)
        AD(I)=(6.*AX(I)/PI)**(1./3.)
1     CONTINUE
!**************  BREAKUP COEFFICIENTS  KIJ *******
      DO 556 I=9,18
        DO 556 J=6,I
          READ (30,801) KIJ(I,J)
801     FORMAT(2X,E12.5)
556   CONTINUE
!************* BREAKUP COEFFICIENTS PLL ********
      DO 555 I=9,18
        DO 555 J=6,I
          DO 555 K=1,I+1
            READ (32,800)PLL(I,J,K)
              IF(I == 13.AND.J == 13.AND.K >= 8)                        &
     &          PLL(I,J,K)=PLL(I,J,K)*PI**2*AD(I)**2*AD(J)**2/4.
555   CONTINUE
800   FORMAT(2X,E12.5)
      DO 855 I=9,18
         DO 855 J=6,I
           PLM(I,J)=0.
             DO 855 K=1,I+1
               PLM(I,J)=PLM(I,J)+PLL(I,J,K)*AX(K)**2
855   CONTINUE
      CLOSE(30)
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
      
      END subroutine DATA
!


end module
