!
! Common blocks used to run the tau bin micro in the 
! 1-D framework

module mphys_tau_bin_declare

  Use typeKind, only : wp
  Use parameters, only : nz, nx, dt, num_h_moments, num_h_bins, &
       nspecies, mom_names, h_names, mom_units, max_char_len, &
       num_aero_moments,num_aero_bins, aero_mom_init

#ifdef USE_MPI
  Use mpi_mod
#endif
  Implicit none

! Control parameters for the model:-
! Those which require a user value:
INTEGER, PARAMETER::IIP = 1  
                                      ! Number of grid points in x-direc
                                      ! Number of grid points in y-direc
#ifdef USE_MPI
!!!!!
! This 'pnx' is temporarily defined. This code is needed to be parallelized.
!!!!!
INTEGER, PARAMETER::JJP = pnx
#else
INTEGER, PARAMETER::JJP = nx  
#endif
                                      ! Number of grid points in z-direc
INTEGER, PARAMETER::KKP = nz  
                                      ! Number of processors used (PEs)
INTEGER, PARAMETER::NPES = 1  
INTEGER, PARAMETER::IBSCATP = 0  
                                      ! Switch for Backscatter 1=ON 0=OF
INTEGER, PARAMETER::NQSCTP = 0  
                                      ! Switch for Backscatter on Q 1=ON
INTEGER, PARAMETER::NBEGSCATP = 0  
                                      ! Timestep to begin Backscatter
INTEGER, PARAMETER::ITVDUVWP = 1  
                     ! Switch for TVD advection on momentum 1=TVD 0=Pias
INTEGER, PARAMETER::ITVDSCALP = 1  
                     ! Switch for TVD advection on scalars 1=TVD 0=Piasc
INTEGER, PARAMETER::IFORSCALP = 1  
                ! Switch for forward stepping of scalars 1=Forward 0=Cen
INTEGER, PARAMETER::IFORUVWP = 1  
                ! Switch for forward stepping of momentum 1=Forward 0=Ce
INTEGER, PARAMETER::IBAROCLP = 0  
                ! Switch for baroclinicity: 1=Geostrophic shear 0=None
INTEGER, PARAMETER::IUSETHP = 1  
                     ! Switch to use TH variable 1=Use 0=Do not use
INTEGER, PARAMETER::IPASTHP = 0  
          ! Switch for passive TH 1=Passive 0=Active - must be 1 if IUSE
INTEGER, PARAMETER::IUSEQP = 1  
                     ! Switch to use Q variables 1=Use 0=Do not use
INTEGER, PARAMETER::IPASQP = 0  
          ! Switch for passive Q 1=Passive 0=Active - must be 1 if IUSEQ
INTEGER, PARAMETER :: Ln2=1
                     ! Ln2 is the number of bins for aerosol. Used the same
                     ! notation as the cloud model for now.
INTEGER, PARAMETER :: LK=34
                     ! Number of bins for cloud. In the model the number of
                     ! cloud bins is LK*2, first 34 is mass and the second
                                ! number concentration
integer, parameter :: lk_cloud = 15 
                     ! the maximum cloud bin seperation
INTEGER, PARAMETER :: LK_big = LK-15 
                   ! Number of large bins used in SCONC and BREAK      
INTEGER, PARAMETER :: NQP_cond = 3 !the number of NQP fields 
                                         !conditional diags are calced on
INTEGER, PARAMETER :: AERO_BIN=1
                     ! position after which ccn occur in Q array
                     !i.e after variable 5 in the Q array CCN bins begin
INTEGER, PARAMETER :: LKCO=19
INTEGER, PARAMETER :: LK4=8
      !cut off for small particles. After Mordy (1959), any particles 
      !lower than 0.12microns, if nucleated are moved to first cloud
      !bin
INTEGER, PARAMETER :: LNUC=Ln2+1
INTEGER, PARAMETER :: LKDD=LK+1
REAL, PARAMETER :: bin_res = 2.0 !bin resolution, i.e. 2.0 will result
                                !in a bin distribution that is mass doubling
REAL, PARAMETER :: DTMIC=2.0 
                     !timestep for bin microphysics 
INTEGER, PARAMETER::IANELP = 0  
          ! Switch to use anelastic equations 1=Anelastic 0=Boussinesq
INTEGER, PARAMETER::IDAMPP = 0  
                     ! Switch to use damping layer 1=ON 0=OFF
INTEGER, PARAMETER::INOVISP = 0  
          ! Switch for viscosity 1=Inviscid 0=Viscosity and Diffusion ON
INTEGER, PARAMETER::IGALOFFP = 1  
                     ! Switch for Galilean transformation 1=OFF 0=ON
INTEGER, PARAMETER::INOSURFP = 1  
          ! Switch for surface fluxes 1=Fluxes OFF 0=Fluxes ON
INTEGER, PARAMETER::ISCBCP = 1 
          ! Flag for type of surface boundary condition for scalars
          ! 1=fixed surface fluxes 2=fixed surface values
INTEGER, PARAMETER::ISATSURFP = 1 
                     ! If =1 then saturated surface for ISCBCP=2
INTEGER, PARAMETER::IADJANELP = 1  
          ! Flags for various options for setting up anelastic profiles
          ! 1 = No adjustment
          ! 2 = ensure P0=PSF by adjusting THREF profile by constant fac
          ! 3 = ensure P0=PSF by adjusting PSF (not advised)
          ! 4 = ensure P0=PSF by adjusting PTOP
          ! Values only relevant if IANELP=1
                                      ! Number of bins for each time ser
INTEGER, PARAMETER::NTIMP = 0  
                                      ! Switch for timeseries 1=ON 0=OFF
INTEGER, PARAMETER::ITSERP = 0  
                                      ! Switch for spectra 1=ON 0=OFF
INTEGER, PARAMETER::ISPECP = 0  
                                      ! Maximum number of spectra
INTEGER, PARAMETER::NSPECP = 0  
                              ! Parameters to control microphysics
     ! I'm making these varaiables...
!INTEGER, PARAMETER::   &
!     IRAINP = 10 &
!      ,ICECLP = 0 &
!      , ISNOWP =0 &
!      , IGRAUP =0 &
!      , INICEP = 0 &
!      , INSNOWP = 0 &
!      , INRAINP = 1 &
!      , INGRAUP = 0 &
!      , IVGRAUP = 0 &
!      , INLIQP = 0
                      ! Flag for microphysics, 0=no rain, 1=Kessler
                      ! 2=Lee 3=Swann 10=three-phase microphysics
                      ! Following must be 0 unless irainp=10
                      ! 1 = use ice cloud mixing ratios
                      ! 1 = use snow mixing ratios
                      ! 1 = use graupel mixing ratios
                      ! 1 = use ice cloud number concentrations
                      ! 1 = use snow number concentrations
                      ! 1 = use rain no. conc.
                      ! 1 = use graupel number concentrations
                      ! 1 = use graupel volume
                      ! 1 = use liquid water no. conc. (not available)
INTEGER, PARAMETER::IREPBARP = 0
          ! Switch for BAR calculations 1=reproducible BARS with arb NPE
          !                             0=non-reproducible with arb NPES
INTEGER, PARAMETER::ITWOFILEP = 1
          ! Switch to dump to one of two files 0=filea only 1=filea and
                             ! Parameters to control radiation
INTEGER, PARAMETER::IRADNP = 0, IRADSWP = 0, IRADLWP = 0, &
 IRAD1DP = 0
                         ! 0 = no radiation, 1 = radiation
                         ! 1 = sw radiation on
                         ! 1 = lw radiation on
                         ! 1 = horizontally homogeneous radiation
! Parameters with default values that can be overridden by the user
                              ! Parameters to control large-scale forcin
INTEGER, PARAMETER::ISUBP = 0, IFORCEP = 0, ITMFORP = 0  
                      ! Switch for large-scale subsidence
                      ! Switch for prescribed forcing
                      ! Switch for time-varying forcing
                              ! Parameters controlling the size of the
INTEGER, PARAMETER::NTMFORP = (120 - 1) * ITMFORP + 1, KTMFORP = &
 (KKP - 1) * ITMFORP + 1
                              ! time-varying forcing arrays
                                     ! Maximum number of time values
                                     ! Maximum number of levels
                             ! Parameters to control diagnostics
INTEGER, PARAMETER::IDGDEEP = 0, ISURFDG = 1, ICLASSP = 0, &
 ILATENTDGP = 0, NUPDOWNMAXP = 1, IVARDGP = 0
                        ! 1 = switch on conditional diagnostics
                        ! 1 = switch on surface diagnostics
                        ! 1 = use column classified diagnostics
                        ! 1 = output microphysics heating rates
                        ! 2 = and process conversion rates
                        ! Max number of conditional partitions
                        ! >0 Output conditional variance diagnostics
                        !    for various moist variables
INTEGER, PARAMETER::NQFBUDGPMAX = 5  
                     ! Max number of Q variables for flux budgets
INTEGER, PARAMETER::MOISTRIP = 2  
           ! Flag for choice of moist Richardson Number calculation
           ! 1=not recommended 2=recommended
INTEGER, PARAMETER::IFBCHGP = 0  
           ! Switch for time-varying surface values 1=ON 0=OFF
INTEGER, PARAMETER::NDGSP = 600  
           ! Maximum number of time-averaged diagnostics
INTEGER, PARAMETER::IEACHDUMPP = 0  
           ! Switch for dump 0=all PEs dump to same file
           !                 1=each PE to dump separately
INTEGER, PARAMETER::ISMTHP = 1  
           ! Flag for choice of time filtering
           ! 0=no filtering 1=Robert filter 2=original (V2.0 and before)
                                           ! Length of APARMS array
INTEGER, PARAMETER::NAPARMSP = 1024  
      ! MAXQSP removed from code at V2.3
                                        ! Length of Monin-Obukov lookup
INTEGER, PARAMETER::LOOKP = 80  
                                           ! A small number
REAL, PARAMETER::SMALLP = 1.E-14  
                                           ! A large number
REAL, PARAMETER::RLARGEP = 1.E37  
INTEGER, PARAMETER::NTIMPDGP = 100  
           ! Dimension of TIMPDG: times for unaveraged diagnostic files
INTEGER, PARAMETER::NTIMRDGP = 100  
           ! Dimension of TIMRDG: times for averaged diagnostic files
INTEGER, PARAMETER::NTIMDUMP = 100  
           ! Dimension of TIMDUM: specified times for model dumps
INTEGER, PARAMETER::NTMHALTP = 10  
           ! Dimension of TIMHALT: times to halt model run
INTEGER, PARAMETER::NTIMHFP = 20  
           ! Maximum number of input times for time-varying surface flux
INTEGER, PARAMETER::NINITP = 20  
           ! Maximum number of input levels for initial profiles
INTEGER, PARAMETER::QLINITP = 0  
           ! =1 THINIT and QINIT(1) are input as theta_l and q_t
           !    and iteration calculates QINIT(1), QINIT(2)
           !    and THINIT as initial qv,ql, and theta values
           ! =0 THINIT and QINIT(1) are input as theta and q_v
INTEGER, PARAMETER::JMINP = 1, JMAXP = JJP  
           ! Lower and upper bounds for certain J loops
INTEGER, PARAMETER::IRADSNW = 3  
         ! Radiation treatment of hydrometeors:
         ! >=1 treat snow, >=2 treat rain, >=3 treat graupel
INTEGER, PARAMETER::IEXPLWP = 0  
!        ! LW radiation calculated using simple exponential scheme
!        ! using Q(IQL) if IEXPLWP = 1
!        ! using Q(1)   if IEXPLWP = 2
INTEGER, PARAMETER::IEXPSWP = 0  
!        ! 1 = SW radiation using EUROCS "sunray" code
                               ! Parameters to control rad diagnostics
INTEGER, PARAMETER::RAD_DIAGS = IRADNP, RAD_DIAGS2 = IRADNP  
                               ! >=2 timeseries of toa & surface fluxes
                               ! >=3 averaged profiles of fluxes
                               ! >=1 2-D fields of heating rates
                               ! >=2 2-D fields of toa & surf fluxes
                               ! >=3 2-D fields of fluxes thro'out
                            ! Parameters to control the animation diagno
INTEGER, PARAMETER::IMOVIEP = 0, JMOVMINP = 1, JMOVMAXP = JJP, &
 KMOVMINP = 2, KMOVMAXP = KKP, IJINCP = 2
                            ! Switch to activate diagnostics =1 ON
                            ! Minimum of range in J direction
                            ! Maximum of range in J direction
                            ! Minimum of range in K direction
                            ! Maximum of range in K direction
                            ! Horizontal grid increment
                                !Parameters to control the NC animation
INTEGER, PARAMETER::IMOVIE_NC = 0, KMINP_MNC = 2, KMAXP_MNC = KKP  
                                !Switch to activate NC movie
                                !Minimum of range in K direction
                                !Maximum of range in K direction
! seed for forte compiler needs to be at least 4 integer
! long - marc 20/11/02
!      INTEGER, PARAMETER :: ISDSZP = 2 ! Maximum pos. size of seed arra
                                         ! Maximum pos. size of seed arr
INTEGER, PARAMETER::ISDSZP = 256  
                                       ! which is found by a call to
                                       ! RANDOM_SEED - intrinsic routine
                                    ! 0=Ice aggregation based on Ferrier
INTEGER, PARAMETER::IAGGP = 0  
                                    !   JAS 51, 249--280.
                                    ! 1=R. Cotton's ice aggregation sche
      ! Aerosol flags: 1=activated, 0=disabled.
                                         ! CCN aerosol
INTEGER, PARAMETER::INCCNP = 0  
                                         ! Deposition ice nuclei
INTEGER, PARAMETER::INDEPP = 0  
                                         ! CCN that are capable of evapo
INTEGER, PARAMETER::INEVFP = 0  
                                         !     freezing
                                         ! Contact-freezing aerosol
INTEGER, PARAMETER::INCONP = 0  
                                         ! Aerosol
INTEGER, PARAMETER::INAERP = 0  
! Parameters calculated from above parameters, not to be overridden by t
! user unless you REALLY know what you're doing !
INTEGER, PARAMETER::ISAVP = 1 - (ISMTHP - 1) * (ISMTHP - 1)  
           ! =1 for Robert, 0 otherwise: collapses 3D arrays to 1D
                                                ! Flag for 3D 0=2D 1=3D
INTEGER, PARAMETER::I3DP = 2 - (IIP + 1) / IIP  
INTEGER, PARAMETER::IIPEP = I3DP * (IIP / NPES) + (1 - I3DP)  
           ! Number of slices per PE for 3-D; =1 for 2-D
INTEGER, PARAMETER::JJTOTP = JJP * (1 + (NPES - 1) * (1 - I3DP) )  
           ! Total number of columns; = JJP for 3-D; =JJP*NPES for 2-D
INTEGER, PARAMETER::JJPEP = JJP / NPES, JJPEMINP = 1 - I3DP, &
 JJPEMAXP = JJPEP * I3DP + 1
           ! parameters used for transposed tridiagonal algorithm
           ! JJPEMINP and JJPEMAXP cover possibility for
           ! 2-D run on single processor
INTEGER, PARAMETER::IDIMMINP = ( - 1 - IBSCATP - 1) * I3DP + 1, &
 IDIMMAXP = (IIPEP + 2 + IBSCATP - 1) * I3DP + 1
           ! Maximum and minimum I dimensions for main field arrays
           ! Effectively the number of halo slices required
INTEGER, PARAMETER::IMINP = 0, IMAXP = IIP * I3DP + JJTOTP * &
 (1 - I3DP) + 1
          ! IMAXP and IMINP are dimensions for PT (transformed pressure)
          ! in POISSON and account for 2-D run (when P not transformed)
INTEGER, PARAMETER::IPMINP = 1 - I3DP, IPMAXP = IIPEP + I3DP  
          ! x-direction dimension of the pressure array
INTEGER, PARAMETER::NSP = 2 * I3DP + 1  
          ! Number of slices held for source (S) fields
INTEGER, PARAMETER::NVISP = 2 * I3DP + 1 + IBSCATP  
          ! Number of slices held for VIS and DIFF fields
INTEGER, PARAMETER::NDISP = 4 * IBSCATP + 1  
          ! Number of slices held for Backscatter dissipation fields
INTEGER, PARAMETER::IFIRSTP = - 1 - 2 * IBSCATP  
          ! Index of first slice for I loop in NNSTEPS
                                               ! Total size of vertical
INTEGER, PARAMETER::LENP = (JJP + 2) * KKP  
INTEGER, PARAMETER::IBSCATTP = IBSCATP * IUSETHP  
          ! Switch for Backscatter on TH (check if using TH)
INTEGER, PARAMETER::IBSCATQP = IBSCATP * IUSEQP  
          ! Switch for Backscatter on Q (check if using Q)
INTEGER, PARAMETER::NFLDSCTP = (IBSCATP + IBSCATTP + NQSCTP * &
 IBSCATQP - 1) * IBSCATP + 1
          ! Total number of fields for Backscatter; momentum counts as 1
!INTEGER, PARAMETER::IMICROP = MIN (1, IRAINP)  
!          ! Switch for microphysics 1=ON 0=OFF
INTEGER, PARAMETER::NMETEORP = 10 !IMICROP + ICECLP + ISNOWP + IGRAUP  
          ! Number of hydrometeors (Q fields with fall velocity)
INTEGER, PARAMETER::NMETP = 10 ! + (NMETEORP - 1) * IMICROP  
          ! Number of hydrometeors =1 if IMICROP=0
!INTEGER, PARAMETER::NAEROSOL = INCCNP + INDEPP + INEVFP + INCONP + &
! INAERP
          ! Number of aerosol fields
INTEGER, PARAMETER :: NQP=AERO_BIN+Ln2+LK+LK+1
!INTEGER, PARAMETER::NQP = 10 !MAX (1 + IUSEQP, IMICROP * (NMETEORP + &
! 2) + INICEP + INRAINP + INGRAUP + INSNOWP + IVGRAUP + NAEROSOL)
          ! Number of Q variables (must be 1 if IUSEQP=0)
INTEGER, PARAMETER::NPRC = 70  
          ! Number of cloud microphysical process terms
INTEGER, PARAMETER::NSERP = 50 ! + IMICROP + 4 * ICECLP + ISNOWP + &
! IGRAUP + 4 * ICLASSP + 3 * IRADLWP + 4 * IRADSWP + 2 * NQP + &
! NMETEORP
          ! Maximum number of timeseries
INTEGER, PARAMETER::NTMP = (NTIMP - 1) * ITSERP + 1, NSRP = &
 (NSERP - 1) * ITSERP + 1
          ! Number of timeseries bins, minimised if ITSERP=0
INTEGER, PARAMETER::INDTVDSP = 2 - (ITVDSCALP + 2) / (ITVDSCALP + &
 1), INDTVDMP = 2 - (ITVDUVWP + 2) / (ITVDUVWP + 1)
          ! To minimise the size of ADV_ arrays if TVD not used

real :: eps ! machine dependent small number, epsilon(1d.0)

COMMON / MOIST / QSAT, QSATFAC, TSTARPR, DQSATDT, TREF, QSATPW  
! Arrays for moist processes
! Note that at V2.3 arrays only used in MOISTRI(2) and Taylor series
! expanded over horizontal mean only. Pointwise values used in LATHEAT
REAL, DIMENSION (KKP) ::QSAT, QSATFAC, TSTARPR, DQSATDT, TREF  
                     ! Saturation mixing ratio (kg/kg)
                     ! 1/(1 + RLVAP_ON_CP*DQSATDT)
                     ! Temperature about which Taylor Expansion
                     ! conducted for qsat calculation
                     ! Gradient of saturation mixing ratio with temperat
                     ! Reference temperature (K)
REAL, DIMENSION (JMINP:JMAXP, KKP) ::QSATPW  
                     ! Pointwise QSAT calculated in LATHEAT and then
                     ! used in MICROPHYS for sat. def. calculation
!
COMMON / MEANS / OLUBAR, OLZUBAR, OLVBAR, OLZVBAR, OLTHBAR, &
 OLZTHBAR, OLQBAR, OLZQBAR
! Contains the true instantaneous horizontal mean fields for:
                            ! Current U mean
REAL :: OLUBAR (KKP)  
                            ! Previous timestep U mean
REAL :: OLZUBAR (KKP)  
                            ! Current V mean
REAL :: OLVBAR (KKP)  
                            ! Previous timestep V mean
REAL :: OLZVBAR (KKP)  
                            ! Current TH mean
REAL :: OLTHBAR (KKP)  
                            ! Previous timestep TH mean
REAL :: OLZTHBAR (KKP)  
                            ! Current Q means
REAL :: OLQBAR (KKP, NQP)  
                            ! Previous timestep Q means


REAL :: OLZQBAR (KKP, NQP)  
!
! Parameters for the three-phase microphysics routines
!
!   Various physical constants:
!   ***************************
REAL, PARAMETER::rhoW = 1000., Rw = 461.5, THcond = 0.0243, &
 diffWV = 2.26E-5, visair = 1.44E-5, Cwater = 4187., Cice = 2093.
                                 ! density of water
                                ! gas conbstant for water vapour
                                ! thermal conductivity of air
                                ! diffusivity of water vapour in air
                                ! kinematic viscosity of air
                                ! specific heat capacity of liquid water
                                ! specific heat capacity of ice
!
!   Intercept and Shape parameters:
!   *******************************
!       RnaX and RnbX are NOT used for double-moment representation
!            instead the number concentration is predicted.
!
REAL, PARAMETER::Rna_R = 1.1E15, Rnb_R = 0., alph_R = 2.5  
! Tuning for tropical convection small snowflakes with higher concentrat
REAL, PARAMETER::Rna_S = 2.E27, Rnb_S = - 3.5, alph_S = 2.5  
! For midlatitude, bigger snow flakes - see Swann'98 (Atmos. Res.:
!      REAL, PARAMETER :: Rna_S=6.E20,Rnb_S=-2,0,alph_S=2.5
! Rutledge and Hobbs'83 :much used crap snow scheme -sensitivity study
!      REAL, PARAMETER :: Rna_S=4.E6,Rnb_S=0.0,alph_S=0.0
! Graupel as in Swann'98 (Atmos. Res.):
REAL, PARAMETER::Rna_G = 5.E25, Rnb_G = - 4.0, alph_G = 2.5  
REAL, PARAMETER::alph_I = 0.0  
!
!   Mass--Diameter relationships:
!   *****************************
!                             [Mass of X particle]=cX*[DiameterX]**dX
REAL, PARAMETER::c_S = 52.36, d_S = 3.  
REAL, PARAMETER::c_R = 523.6, d_R = 3.  
! Tuning for tropical convection, denser graupel (due to high QL:
!      REAL, PARAMETER :: c_G=366.5,d_G=3.
                         !c_G not used if IVGRAUP=1(triple-moment)
! For midlatitude soft graupel:
REAL, PARAMETER::c_G = 261.8, d_G = 3.  
                         !c_G not used if IVGRAUP=1(triple-moment)
REAL, PARAMETER::c_I = 104.0, d_I = 3.  
!
! Fallspeed--Diameter relationships
! *********************************
!          [Vel of X particle]=aX*[DiameterX]**bX*exp(-[DiameterX]*fX)
REAL, PARAMETER::a_R = 362., b_R = 0.65, f_R = 0.  
                                                     !f_S must be 0.0
REAL, PARAMETER::a_S = 4.84, b_S = 0.25, f_S = 0.  
                                                     !f_I must be 0.0
REAL, PARAMETER::a_I = 71.34, b_I = 0.6635, f_I = 0.  
                                                             !f_G must b
REAL, PARAMETER::a_G = 253., b_G = 0.734, f_G = 0., g_G = 0.422  
!
!   Various constants for ice parameterization:
!   *******************************************
! mass of pristine ice crystal (equivalent to radius of order 1 micron).
! New value introduced with R.Cotton's ice aggregation
!      REAL, PARAMETER :: RMI0=1.E-15
REAL::RMI0 = 1.E-18  
! constant in aerosol concentration eq.
REAL, PARAMETER::RNa0 = 2.E3  
! maximum mean Diameter of ice crystals before conversion to snow
REAL, PARAMETER::DImax = 0.0003  
! Diameter of ice crystals that convert to snow
REAL, PARAMETER::DI2S = 0.00033  
! minimum lambI before conversion to snow
REAL, PARAMETER::RlambImin = (1. + alph_I + d_I) / DImax  
! max mass-mean diameter of snow before breakup
REAL, PARAMETER::DSmax = 0.004  
! minimum lambS before breakup
REAL, PARAMETER::RlambSmin = (1. + alph_S + d_S) / DSmax  
!
!   Thresholds for autoconversion and rate coefficient for PRAUT:
!   *************************************************************
! liquid cloud droplet concentration
REAL, PARAMETER::RNc = 50E6 !240E6  
! droplet diameter required for autoconversion
REAL, PARAMETER::DliqtoR = 0.00002  
! Not very sensitive to this for single-moment rain:
REAL, PARAMETER::autoR = 0.001  
! Radius of newly formed raindrops
REAL, PARAMETER::Ro = 0.000025  
      ! Coefficients used in most recent microphysics conversions
REAL, PARAMETER::ACFAC11 = 2.20, QREXP11 = 0.875, ACFAC12 = 7.58, &
 QLEXP12 = 1.029, QREXP12 = 1.042, ACFAC13 = 3.9176, QREXP13 = &
 0.9702, AUTFAC12 = 1.44, AUTEXP12 = 2.36, AUTFAC14 = 1350, &
 AUTEXP14 = 2.47, RNCEXP14 = - 1.79
!
!   Collection efficiencies:
!   ************************
REAL, PARAMETER::ERW = 1., ERG = 1., ERS = 1., ERI = 1., ESW = 1., &
 EGW = 1., EIW = 1., EGS_WET = 1., EGI_WET = 1.
REAL, PARAMETER::ERR = 1.0  
REAL, PARAMETER::EGS1_DRY = 0.2, EGI1_DRY = 0.2, ESI1 = 0.2, EII1 &
 = 0.2
REAL, PARAMETER::EGS2_DRY = 0.08, EGI2_DRY = 0.08, ESI2 = 0.08, &
 EII2 = 0.08
REAL, PARAMETER::ESS1 = 0.2, ESS2 = 0.08  
!
!   Constants used in Bigg freezing:
!   ********************************
REAL, PARAMETER::ABIGG = 0.66, BBIGG = 100.  
!
!   These are multipliers for standard rates of:
!   ********************************************
! HALMOS - Hallett-Mossop process
!          (default rate equivalent to 350 splinters per milligram of ri
! CONTACT- contact nucleation rate, default value now 0.0001,due to
!          Young(73) vastly overestimating this process.
! AUTOCON- threshold QL for onset of autoconversion. It is better to adj
!          this by means of the threshold droplet size DliqtoR.
! FLETCHER- fletcher primary nucleation, set to 0.0 to use Meyers'
! RMEYER-   Meyers'  primary nucleation, set to 0.0 to use Fletcher's
! RHOMOG-   Switch for homogeneous nucleation of liquid water at -38C,
!            set to zero to de-activate this process as a source of ice
!
! "Homogeneous" nucleation only - set RMEYER=FLETCHER=0, RHOMOG=1.0
!
REAL, PARAMETER::HALMOS = 1.0, CONTACT = 0.0001, FLETCHER = 0.0, &
 RMEYER = 1.0, RHOMOG = 1.0, AUTOCON = 1.0
!
!   Not applicable until D-M rain (INRAINP=1) option is coded:
!   ***********************************************
                                       ! mean diameter of rain from melt
REAL, PARAMETER::DRMELT = 0.001  
!
!   Others.. only applicable if IVGRAUP=1:
!   ******************************************
                                       !    Each of these
REAL, PARAMETER::RHOGAUT = 200.  
                                       !    is the density of
REAL, PARAMETER::RHOGFR = 900.  
                                       !    each source of graupel.
REAL, PARAMETER::RHORACI = 900.  
REAL, PARAMETER::RHOGACS = 200.  
                                       !    Each sink to the graupel mas
REAL, PARAMETER::RHOGACI = 200.  
                                       !    is assumed not to change the
REAL, PARAMETER::RHOSACR = 500.  
                                       !    graupel density.
REAL, PARAMETER::RHOGDEP = 200.  
REAL, PARAMETER::RHOGACW = 500.  


REAL, PARAMETER::RHOGACR = 800.  
!   Factors used in Grapuel fall/linit calculations
!   ***********************************************

REAL, PARAMETER::gl1fac = 1.e-12, gl2fac = 1.0 / 100., gl3fac = &
 1.0 / 910., gl4fac = 6.0, gl5fac = 1.e-12, gl6fac = - 0.625, &
 gl7fac = 0.578
!   Various factors for triple-moment graupel
!   *****************************************

REAL ::gfacmass, gfacfall, gfacwet  
                      ! (Actual rhoG)/(rhoG_from_c_G)
                      ! gfacmass**0.578
                      ! (1.E-12 + gfacmass)**(-0.625)


COMMON / TMGRAUP / gfacmass, gfacfall, gfacwet  
!   Thresholds for Activating Microphysics
!   **************************************


REAL ::qsmll = 1.e-8, qsmll2 = 1.e-4, qvsmll = 1.e-13
!REAL :: qsmlli = 0.99 * RMI0   ! changed to a variable
REAL :: qsmlli                  ! set in set_constants
                            ! critical hydrometeor content for calc
                            ! should be combined with similar in diags
                            ! small but a bit bigger REVIEW!
!     & , qsmlli = 1.e-11   ! Small critical value for ice REVIEW!
                            ! Introduced with R.Cotton's ice aggregation
                            ! smaller critical value of olqbar REVIEW
!   Useful Thermodynamic Quantities
!   *******************************


REAL, PARAMETER::TK0C = 273.15, pa2mb = 1. / 100.0  
                          ! Temperature of freezing in K
                          ! Factor in converting pascals to mb
!   Default PSD Parameters
!   **********************




REAL, PARAMETER::Lambda_default = 1.0e7, N0_default = 1.0e7  
                                ! Default value of Lambda_x
                                ! Default value of N0_x
!**********************************************************************
! ----------------------------------------------------------------------
! Microphysical Process Rates
! ----------------------------------------------------------------------
! These  process rates are described in Swann (1998), PhD Thesis,
! Meteorology Dept, Reading Univ.
! See also the Scientific Documentation of the LEM, APR Turbulence and
! Diffusion Note 276.
! Also on-line within the Met Office at:
!           http://www-nwp/~lemdev/CONV/LEM/Scientific/sci.shtml
! ----------------------------------------------------------------------
!
!
! Autoconversion Process Rates
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PRaut, RRaut, PSaut, RSaut, PGaut, &
 RGaut
                           ! Conversion of cloud liquid water to rain.
                           ! Conversion of cloud ice to snow.
                           ! Conversion of rimed snow to graupel.
! Aggregation and Breakup
! ----------------------------------------------------------------------


REAL, DIMENSION (JJP) ::RRacR, RRbrk, RIacI, RSacS, RSbrk  
                           ! Aggregation of rain drops.
                           ! Breakup of rain drops.
                           ! Cloud ice aggregation.
                           ! Snow aggregation.
                           ! Breakup of snow flakes.
! Collection of Liquid Cloud
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PRacW, RRacW, PIacW, RIacW, PSacW, &
 RSacW, PGacW, RGacW
                           ! Collection by rain.
                           ! Collection by cloud ice.
                           ! Collection by snow.
                           ! Collection by graupel.
! Collisions between Rain and Cloud Ice
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PRacI_S, PRacI_G, PIacR_S, PIacR_G, &
 RRacI_S, RRacI_G, RIacR_S, RIacR_G
                           ! Rate at which cloud ice is converted to
                           ! snow by collisions with rain.
                           ! Rate at which cloud ice is converted to
                           ! graupel by collisions with rain.
                           ! Rate at which rain is converted to snow
                           ! by collisions with cloud ice.
                           ! Rate at which rain is converted to graupel
                           ! by collisions with cloud ice.
                           ! Since all collisions are assumed to be
                           ! binary, this is just the same as RIacR_S.
                           ! See comment for RRacI_S.
                           ! Rate at which binary collisions between
                           ! rain drops and cloud ice particles creates
                           ! new snow flakes.
                           ! Rate at which binary collisions between
                           ! rain drops and cloud ice particles creates
                           ! new graupel particles.
! Collisions between Rain and Ice Aggregates (Snow and Graupel)
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PSacR, PRacS, RSacR, PGacR, RGacR  
                           ! Rate at which rain is collected by snow, to
                           ! form graupel (or snow if graupel not used).
                           ! Rate at which snow is collected by rain, to
                           ! form graupel (or snow if graupel not used).
                           ! Rate at which binary collisions between
                           ! rain and snow form graupel particles.
                           ! NB. Both above are sources of graupel.
                           ! Rate at which graupel collects rain.
! Collisions Among the Ice Categories
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PSacI, RSacI, PGacI, RGacI, PGacS, &
 RGacS
                           ! Rate at which cloud ice is collected by
                           ! snow.
                           ! Rate at which cloud ice is collected by
                           ! graupel.
                           ! Rate at which graupel collects snow.
! Graupel Growth Rates and Shedding
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PGwet, PGdry, PGshd, RGshd  
                           ! Maximum rate of growth of graupel, assuming
                           ! graupel surface is wet, ie. at 0 degrees C.
                           ! Graupel growth rate assuming graupel
                           ! surface is dry, ie. < 0 degrees C.
                           ! Excess liquid water shed by graupel when
                           ! the wet growth limit is reached.
! Phase Changes
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PRevp, RRevp, PGfrR, RGfrR, PIfrW, &
 RIfrW, PImlt, RImlt, PSmlt, RSmlt, PGmlt, RGmlt, PIsub, RIsub, &
 PSsub, RSsub, PGsub, RGsub, PIdep, RIdep, PSdep, RSdep, PGdep, &
 RGdep
                           ! Evaporation of rain.
                           ! Freezing of rain to form graupel.
                           ! Freezing of cloud liquid water.
                           ! Melting of cloud ice.
                           ! Melting of snow.
                           ! Melting of graupel.
                           ! Sublimation of cloud ice.
                           ! Sublimation of snow.
                           ! Sublimation of graupel.
                           ! Deposition onto cloud ice.
                           ! Deposition onto snow.
                           ! Depsoition onto graupel.
! Cloud Ice Nucleation
! ----------------------------------------------------------------------

REAL, DIMENSION (JJP) ::PIprm, RIprm, PIcnt, RIcnt, PIhal, &
 RIhal
                           ! Primary ice nucleation.
                           ! Freezing of liquid cloud in contact with
                           ! aerosols.
                           ! Hallet-Mossop process (rime splintering).
! R. Cotton's Ice Aggregation Scheme
! ----------------------------------------------------------------------


REAL, DIMENSION (JJP) ::PIdcf, PIEvf, PIFrH, RCCNI  
                           !
                           !
                           !
                           !
! Microphysics increments to Q variables
! ----------------------------------------------------------------------


REAL, Dimension (JJP) ::dqV, dqL, dqR, dqI, dqS, dqG, dqNL, &
 dqNR, dqNI, dqNS, dqNG, dqVG, dqNCCN, dqNDep, dqNEvF, dqNCon, &
 dqNAer
! Common blocks for passing microphysics tendencies & process rates
! ----------------------------------------------------------------------

!COMMON / PROCRATES / PRaut, RRaut, PSaut, RSaut, PGaut, RGaut, &
! RRacR, RRbrk, RIacI, RSacS, RSbrk, PRacW, RRacW, PIacW, RIacW, &
! PSacW, RSacW, PGacW, RGacW, PRacI_S, PRacI_G, PIacR_S, PIacR_G, &
! RRacI_S, RRacI_G, RIacR_S, RIacR_G, PSacR, PRacS, RSacR, PGacR, &
! RGacR, PSacI, RSacI, PGacI, RGacI, PGacS, RGacS, PGwet, PGdry, &
! PGshd, RGshd, PRevp, RRevp, PGfrR, RGfrR, PIfrW, RIfrW, PImlt, &
! RImlt, PSmlt, RSmlt, PGmlt, RGmlt, PIsub, RIsub, PSsub, RSsub, &
! PGsub, RGsub, PIdep, RIdep, PSdep, RSdep, PGdep, RGdep, PIprm, &
! RIprm, PIcnt, RIcnt, PIhal, RIhal, PIdcf, PIEvf, PIfrH, RCCNI

COMMON / QTEND / dqV, dqL, dqR, dqI, dqS, dqG, dqNL, dqNR, dqNI, &
 dqNS, dqNG, dqVG, dqNCCN, dqNDep, dqNEvF, dqNCon, dqNAer
COMMON / RAIN / QTHRESH, EXPZ_ON_SCL, ABLTEMP, PUDDLE, bin_pcpt, &
     bin_ppt_level
! Additional arrays used for the warm rain microphysics calculations
                            ! Autoconversion threshold = 0.001/RHON (kg/
REAL :: QTHRESH (KKP)  
                            ! Factor for accretion calculation
REAL :: EXPZ_ON_SCL (KKP)  
                            ! Factor for evaporation calculation
REAL :: ABLTEMP (KKP)  
! And diagnostic array for all microphysics:

REAL :: PUDDLE (JMINP:JMAXP, NMETP, 0: (IIPEP + 1) )
         ! Diagnostic of total surface hydrometeor fallout for that
         ! model grid point over the entire integration (kg/m2)
REAL, DIMENSION(JMINP:JMAXP) :: bin_pcpt
REAL, DIMENSION(KKP,JMINP:JMAXP) :: bin_ppt_level

!
! Contains values of all physical constants as parameters
!      CP, R, RATIO_MOL_WTS, G, RLVAP, RLSUB,
! and functions of the above
!      RLVAP_ON_CP, RLSUB_ON_CP, R_ON_CP, G_ON_2
! this is to prevent any of these values being accidently changed
REAL, PARAMETER::CP = 1005.  
                   ! Specific heat of gas at constant pressure
REAL, PARAMETER::R = 287.05  
                   ! Universal gas constant
REAL, PARAMETER::RATIO_MOL_WTS = 1.608  
                   ! Ratio of molecular weights of air and water vapour
REAL, PARAMETER::G = 9.81  
                   ! Acceleration due to gravity
REAL, PARAMETER::RLVAP = 2.501e6  
                   ! Latent heat of vapourisation
REAL, PARAMETER::RLSUB = 2.834e6  
                   ! Latent heat of sublimation
REAL, PARAMETER::RLFUS = RLSUB - RLVAP  
                   ! Latent heat of fusion
REAL, PARAMETER::RLVAP_ON_CP = RLVAP / CP  
REAL, PARAMETER::RLSUB_ON_CP = RLSUB / CP  
REAL, PARAMETER::RLFUS_ON_CP = RLSUB_ON_CP - RLVAP_ON_CP  
REAL, PARAMETER::R_ON_CP = R / CP  
REAL, PARAMETER::G_ON_2 = 0.5 * G  
! Time data which is useful to keep constant (formally in RAD_MOD)

REAL, PARAMETER::secinhr = 3600., hrinday = 24., secsinday = &
 86400., rsecsinday = 1. / secsinday
                                ! second in hour
                                ! hours in day
                                ! seconds in day
!
COMMON / MICROCOEF / physcal, phys1d  
! Coefficients used in 3-phase microphysics
INTEGER, PARAMETER::nphyscalp = 11, nphys1dp = 78  
                           ! Number of values of physcal
                           ! Number of phys1d
REAL, DIMENSION (nphyscalp) ::physcal  
                           ! Exchange coefficients used in microx

REAL, DIMENSION (kkp, nphys1dp) ::phys1d  
                           ! Coefficient also fn of height
      !
COMMON / PARAMS / DXX, DYY, VK, FCORIOL, VG0, UG0, TSMTH, BETAM, &
 GAMMAM, BETAH, GAMMAH, DFBMAX, CQ, Z0, PSF, PSFR, KGD, HGD, ZZTOP, &
 THREF0, ZNREFGD, THREFGD, ZNREF_READ, THREF_READ, ALPHAH, Z0TH, &
 DUGDZ, DVGDZ, DUGDZFAC, DVGDZFAC, F_ON_G, RNHPTS, TSMTHC1, &
 TSMTHC2, ZLOGM, ZLOGTH, CZN, COCVIS, VK_ON_ZLOGM, ZLOGTH_ON_VK, &
 !EXPTFC, TIMEAV, DTM_X2, DTM_X4, TIME, TIMEDG, DTM, DTMOLD, DTMNEW, &
 EXPTFC, TIMEAV, DTM_X2, DTM_X4, TIMEDG, DTM, DTMOLD, DTMNEW, &
 RINCMAX, DTMMAX, DTMMIN, UGAL, VGAL, FTH, FQ, FBUOY, FBUOYNEW, &
 SHFLX_SEN, SHFLX_LAT, SCT, SCTT, SCTQ, RHOBOUS, THCONA, THCONB, &
 TOL, RDTM, R_2DTM, THSURF, QSURF, THVSURF, EXP_LW, LWTOP_IN, &
 LWBASE_IN, BHBC, CMBC, RCMBC, RHMBC, DDBC, R2DDBC, DDBC_X4, EECON, &
 TSTRCONA, TSTRCONB, ELLMOCON, X4CON, XX0CON, Y2CON, YY0CON
! A wide selection of model variables defining the numerics and
! physics of the model, mostly input via the namelists
! from Namelist NUMERICS:-
                   ! Grid spacing in x-direction (m)
REAL :: DXX  
                   ! Grid spacing in y_direction (m)
REAL :: DYY  
                   ! Time-smoothing factor
REAL :: TSMTH  
                   ! Incremental factor for increasing timestep
REAL :: RINCMAX  
                   ! Maximum value for timestep (s)
REAL :: DTMMAX  
                   ! Minimum value for timestep (s)
REAL :: DTMMIN  
! from Namelist DYNAMICS:-
                   ! Coriolis parameter (s-1)
REAL :: FCORIOL  
                    ! Components of surface geostrophic wind (ms-1)
REAL :: UG0, VG0  
                        ! Rate of change of components of
REAL :: DUGDZ, DVGDZ  
                        ! geostrophic wind with height (s-1)
! from Namelist PHYSICS:-
                   ! von Karman constant (default value=0.4)
REAL :: VK  
                    ! Monin-Obukhov coefficents; default values:
REAL :: BETAM, GAMMAM, BETAH, GAMMAH, ALPHAH  
                    ! =4.8
                    ! =19.3
                    ! =7.8
                    ! =12.0
                    ! =1.0
                   ! Maximum change in buoyancy flux between iterations
REAL :: DFBMAX  
                   ! in CHGBUOY.
                   ! Coefficients for "virtual" contributions to
REAL :: CQ (NQP)  
                   ! buoyancy from Q fields.
                 ! ratio of the circumference of a circle to its diamete
!REAL :: PI  
! from Namelist INPUT:-
                          ! surface roughness length for momentum (m)
REAL :: Z0  
                          ! surface roughness length for scalars (m)
REAL :: Z0TH  
                          ! surface pressure (Pa)
REAL :: PSF  
                          ! surface reference pressure (Pa)
REAL :: PSFR  
                          ! density for Boussinesq run (kg/m3)
REAL :: RHOBOUS  
                          ! specified surface sensible heat flux (W/m2)
REAL :: SHFLX_SEN  
                          ! specified surface latent heat fluxes (W/m2)
REAL :: SHFLX_LAT (NQP)  
                          ! specified value of surface theta (K)
REAL :: THSURF  
                          ! specified value fo surface vapour mixing rat
REAL :: QSURF  
! from Namelist GRID:-
                          ! Array of heights for setting vertical grid (
REAL :: HGD (NINITP)  
                          ! Gridpoints for which HGD values correspond
INTEGER :: KGD (NINITP)  
                       ! Height of domain top (m)
REAL :: ZZTOP  
! from Namelist THPROF:-
                          ! Reference potential temnperature for
REAL :: THREF0  
                          ! Boussinesq run (IANELP=0)
                               ! Array of heights for setting reference
REAL :: ZNREF_READ (NINITP)  
                               ! potential temperature profile (THREF)
                               ! Corresponding array of values for setti
REAL :: THREF_READ (NINITP)  
                               ! THREF
! from Namelist SUBMODEL:-
                         ! Backscatter coefficient for momentum
REAL :: SCT  
                         ! Backscatter coefficient for TH
REAL :: SCTT  
REAL :: SCTQ (NQP)  
! other variables
                              ! Array of heights for THREF calculation (
REAL :: ZNREFGD (NINITP + 1)  
                              ! THREF values corresponding to ZNREFGD (K
REAL :: THREFGD (NINITP + 1)  
                        ! =DUGDZ*THREF0*F_ON_G ) for calculation of
REAL :: DUGDZFAC  
                        ! =DVGDZ*THREF0*F_ON_G ) baroclinicity
REAL :: DVGDZFAC  
                        ! FCORIOL/G
REAL :: F_ON_G  
                        ! Reciprocal of number of horizontal
REAL :: RNHPTS  
                        ! grid points = 1/REAL(IIP*JJTOTP)
                        ! 1 - TSMTH    ) used for timestep
REAL :: TSMTHC1  
                        ! 2*TSMTH - 1  ) smoothing
REAL :: TSMTHC2  
                        ! Galilean transformation for u (m/s)
REAL :: UGAL  
                        ! Galilean transformation for v (m/s)
REAL :: VGAL  
! following used in LOWERBC
                        ! LOG(1 + ZN(2)/Z0)
REAL :: ZLOGM  
                        ! LOG(1 + ZN(2)/Z0TH)
REAL :: ZLOGTH  
                        ! DZN(2)*0.5
REAL :: CZN  
                        ! Coefficient for viscous Courant number in LOWE
REAL :: COCVIS  
                        ! VK/ZLOGM
REAL :: VK_ON_ZLOGM  
                        ! ZLOGTH/VK
REAL :: ZLOGTH_ON_VK  
                        ! Convective limit = -0.5
REAL :: EXPTFC  
! Time variables
                        ! Current model time (s)
!REAL :: TIME  
                        ! Current model timestep (s)
REAL :: DTM  
                        ! Current time over which diagnostics averaged
REAL :: TIMEAV  
                        ! DTM * 2
REAL :: DTM_X2  
                        ! DTM * 4
REAL :: DTM_X4  
                        ! Time advanced since AVDG last called
REAL :: TIMEDG  
                        ! Timestep on last model step (s)
REAL :: DTMOLD  
                        ! Timestep for next model step (s)
REAL :: DTMNEW  
                        ! Reciprocal of timstep
REAL :: RDTM  
                        ! Reciprocal of twice the timestep
REAL :: R_2DTM  
! Flux and buoyancy variables
                        ! Surface temperature flux
REAL :: FTH  
                        ! Surface moisture flux
REAL :: FQ (NQP)  
                        ! Surface buoyancy flux
REAL :: FBUOY  
                        ! NOT REALLY SURE
REAL :: FBUOYNEW  
! variables for fixed temperature surface boundary condition
                        ! RATIO_MOL_WTS*THREF0 (or THREF(K))
REAL :: THCONA  
                        ! THCONA-THREF0 (or equiv for IANEL=1)
REAL :: THCONB  
                        ! Tolerence for surface flux iteration
REAL :: TOL  
                        ! Surface virtual potential temperature
REAL :: THVSURF  
                        ! ALPHAH*ZLOGTH
REAL :: BHBC  
                        ! BETAM*ZN(2)*G*VK/THVSURF
REAL :: CMBC  
                        ! 1/CMBC
REAL :: RCMBC  
                        ! BETAH*(ZN(2)+Z0-Z0TH)/(BETAM*ZN(2))
REAL :: RHMBC  
                        ! ZLOGM*(BHBC-RHMBC*ZLOGM)
REAL :: DDBC  
                        ! 0.5/DDB
REAL :: R2DDBC  
                        ! 4.*DDBC
REAL :: DDBC_X4  
                        ! 2.*RHMBC*ZLOGM-BHBC
REAL :: EECON  
                        ! VK/BHBC
REAL :: TSTRCONA  
                        ! VK/ALPHAH
REAL :: TSTRCONB  
                        ! THVSURF/(G*VK)
REAL :: ELLMOCON  
                        ! GAMMAM*(ZN(2)+Z0)
REAL :: X4CON  
                        ! GAMMAM*Z0
REAL :: XX0CON  
                        ! GAMMAH*(ZN(2)+Z0)
REAL :: Y2CON  
                        ! GAMMAH*Z0TH
REAL :: YY0CON  
! variables for exponential radiation
                        ! LW exponential decay factor
REAL :: EXP_LW  
                        ! cloud top LW flux (Wm-2)
REAL :: LWTOP_IN  
                        ! cloud base LW flux (Wm-2)


REAL :: LWBASE_IN  
!
REAL ::a_IX, b_IX, c_IX, d_IX  

REAL ::a_SX, b_SX, c_SX, d_SX  

REAL, DIMENSION (7) ::GamFN_I, GamFN_S  

REAL ::dndt_I, dndt_S  

COMMON / AGGPRM / a_IX, b_IX, c_IX, d_IX, a_SX, b_SX, c_SX, d_SX, &
 dndt_I, dndt_S, GamFN_I, GamFN_S
      ! Additional variables used in the deep diagnostics
COMMON / dgstore / l_dodgs, l_updwn4d, dth_dt, dq_dt, dthl_dt_p, dqt_dt_p, &
 dmse_dt_p, olqlbar!, fallq, procrate, nupdown, cupdown
                                                 ! redundant icl,
      !
                         ! number of domain decompositions
!INTEGER ::nupdown  
      !
REAL, DIMENSION (kkp) ::olqlbar  
REAL, DIMENSION (jjp, kkp) ::dth_dt, dthl_dt_p, dmse_dt_p  
REAL, DIMENSION (jjp, kkp, nqp) ::dq_dt!, fallq  
REAL, DIMENSION (jjp, kkp) ::dqt_dt_p  
!REAL, DIMENSION (jjp, kkp, nprc) ::procrate  
      !!
LOGICAL, DIMENSION (nupdownmaxp, jjp, kkp, iipep) ::l_updwn4d  
                       ! flag for each grid point and each catagory

LOGICAL ::l_dodgs  
                       ! Flag for diagnostic run or not

!CHARACTER(LEN = 3), DIMENSION (nupdownmaxp) ::cupdown  
                       ! Name of each category
common / dgchars / chav, chts, chsp, ctemp1, ctemp2, ctemp3, &
 procchar
      !
                                                       !

CHARACTER(LEN = 7), DIMENSION (nprc) ::procchar  
CHARACTER(LEN = 11), DIMENSION (ndgsp) ::chav  
                   !
CHARACTER(LEN = 9), DIMENSION (nserp) ::chts  
                   !
                                                                      !
CHARACTER(LEN = 10), DIMENSION ( (nspecp - 1) * ispecp + 1) :: &
 chsp
      ! temporary characters
CHARACTER(LEN = 3) ::ctemp3  
CHARACTER(LEN = 2) ::ctemp2  

CHARACTER(LEN = 1) ::ctemp1  
!
! History fields for homogeneous freezing

INTEGER, PARAMETER::IIP_CCN = MAX (1, IIP * INCCNP), JJP_CCN = &
 MAX (1, JJP * INCCNP), KKP_CCN = MAX (1, KKP * INCCNP)
REAL ::TC_2 (1:JJP_CCN, 1:KKP_CCN, 1:IIP_CCN), TC_1 (1:JJP_CCN, &
 1:KKP_CCN, 1:IIP_CCN), sW_1 (1:JJP_CCN, 1:KKP_CCN, 1:IIP_CCN)

REAL ::DT_2, DT_1  

!!!Declarations for TAU bin model and interface
  real :: rdt
  real :: prefn(jjp,KKP)
  real :: dzn(KKP)
  real :: prefrcp(jjp,KKP)
  real :: rprefrcp(jjp,KKP)
  real :: rhon(KKP)
  real :: rdz_on_rhon(KKP)
  integer::iqv = 0
!!$  integer::iql = 0  
!!$  integer::iqr = 0  
!!$  integer::iqs = 0  
!!$  integer::iqg = 0  
!!$  integer::iqi = 0  
!!$  integer::iqn = 0  
!!$  integer::iqng = 0
!!$  integer::iqns = 0  
!!$  integer::iqnr = 0  
!!$  integer::iqvg = 0  
!!$  integer::iqnl = 0 
!!$  integer::iqndep = 0  
!!$  integer::iqnccn = 0  
!!$  integer::iqnevf = 0 
!!$  integer::iqncon = 0
!!$  integer::inbinql = 0 
  integer::iqss = 0 

!!$  integer i
!!$
!!$  !Logical switches 
!!$  logical :: l_ice=.False.
!!$  logical :: micro_unset=.True.

  !Parameters now treated as variables
  integer ::         &
         IRAINP  = 0 &
       , ICECLP  = 0 &
       , ISNOWP  = 0 &
       , IGRAUP  = 0 &
       , INICEP  = 0 &
       , INSNOWP = 0 &
       , INRAINP = 0 &
       , INGRAUP = 0 &
       , IVGRAUP = 0 &
       , INLIQP  = 0

  !Parameters for calling bin micro treated as vars
  integer ::         &
         IMICROBIN  = 0 &
       , IRAINBIN   = 0 &
       , IG_FORCE   = 1 &
       , IINHOM_mix = 0 &
       , IEVAP_MIX_JR06 = 0 
  INTEGER :: arr_loop
  integer,DIMENSION(Ln2):: &
          IAERO_BIN=(/(arr_loop, &
                     arr_loop = AERO_BIN+1,AERO_BIN+LN2)/)
  integer,DIMENSION(LK):: &
          ICDKG_BIN=(/(arr_loop, &
                     arr_loop = AERO_BIN+LN2+1,AERO_BIN+LN2+LK)/)
  integer,DIMENSION(LK):: &
          ICDNC_BIN=(/(arr_loop, &
                     arr_loop = AERO_BIN+LN2+LK+1,AERO_BIN+LN2+LK+LK)/)

!!$  type qindex
!!$     integer :: ispecies ! Species index
!!$     integer :: imoment  ! moment index
!!$  end type qindex
!!$     
!!$  type(qindex), allocatable :: qindices(:)
!!$  integer :: nqs ! number of qindices (possibly different from nqp)

  
  integer:: iq, ih, imom
  character(max_char_len) :: name, units
  real :: value


!*COMDECK AEROSOL
!*/AEROSOL COMMON BLOCK - contains the physical properties of the 
!*/aerosol for kohler calculations in readnuc
common /AEOROSOL/rhos,nu,amss,rhos2,nu2,amss2
 REAL RHOS!density of fine aerosol (assumed to be ammonium sulphate)
 REAL NU !number of ions resulting from dissociation of ammonium
 !sulphate
 REAL AMSS!molecular weight of ammonium sulphate 
 REAL RHOS2!density of coarse aerosol (assumed to be NaCl)
 REAL NU2!number of ions resulting from dissociation of NaCl
 REAL amss2!molecular wieght of NaCl
!
!*COMDECK AEROPROF
!*/AEROPROF COMMON BLOCK - contains info on the initial aerosol
!*/profile
COMMON/AEROPROF/RN,CCNINIT,BXsul, BXnacl,xn,xa
 REAL,DIMENSION(lnuc)::RN !aerosol bin boundaries
 !set in BLKDATA
 REAL,DIMENSION(Ln2) ::CCNINIT!initial aerosol distribution
  !from readnuc, read from external file
 REAL,DIMENSION(Ln2) :: xn !aerosol number
 REAL,DIMENSION(lnuc) :: xa !aerosol mass
 
! KiD : moved mean and stdev of aersol dist from module_mp_tau to here
!       so that value can be initialised in namelist. 
 REAL DG1 ! mean aerosol diameter in cms, converted in bin_init
 REAL SG1 ! standard deviation of the aerosol distribution
!
 REAL BXsul!solute term for fine particles (assumed to 
 !be ammonium sulphate)
 REAL BXnacl!solute term for coarse particles (assumed to be
  !NaCl)
! 
!*COMDECK BA
COMMON/BA/CCN
 REAL,DIMENSION(JMINP:JMAXP,KKP,Ln2)::CCN
!CCN is new number of aerosol following activation
! 
!*COMDECK CC
COMMON/CC/XKK1,DUS,DUS1
 REAL,DIMENSION(LK) :: XKK1
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: DUS,DUS1
! 
!*COMDECK X123
common/X123/X2,X12,X13,X23
 REAL,DIMENSION(LKDD)::X2,X12,X13,X23
! 
!*COMDECK DYN
COMMON/DYN/DS
REAL,DIMENSION(JMINP:JMAXP,KKP):: DS
   !absolute value of supersat liquid in Yan's model
! 
!*COMDECK SD
COMMON/SD/ANK,AMK,AMN,AMM
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: ANK!number concentration
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: AMK!mass concentration
 REAL,DIMENSION(LK) :: AMN,AMM
!
!*COMDECK SDOLD
COMMON/SDOLD/ANKOLD,AMKOLD,AMKORIG,ANKORIG,CCNOLD 
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: ANKOLD !no. conc previous dt
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: AMKOLD !mass previous dt
!ANKOLD and AMKOLD are updated after activation
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: AMKORIG !no. conc previous dt
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK):: ANKORIG !mass previous dt
 !AMKORIG and ANKORIG are not updated during microphys, these
 !values are used to work out the source term
 REAL,DIMENSION(JMINP:JMAXP,KKP,LN2) :: CCNOLD

 !Aerosol number concentration before activation
! 
!*COMDECK CON 
COMMON/CON/AP,XI2,XI,YYY,ZZZ,DDS,XI23,ZI
 REAL AP
 REAL XI2
 REAL XI
 REAL YYY
 REAL ZZZ
 REAL DDS
 REAL XI23
 REAL ZI
! 
!*COMDECK WM
COMMON/WM/AM1,AN1,AN2
REAL,DIMENSION(JMINP:JMAXP,KKP)::AM1,AN1
REAL(wp),DIMENSION(JMINP:JMAXP,KKP)::AN2
! 
!*COMDECK QZ 
 COMMON/QZ/QST
 REAL,DIMENSION(JMINP:JMAXP,KKP):: QST
 !saturation specific humidity for bin model 
! 
!*COMDECK TI 
 COMMON/TI/DIEMD,DIEND,DSEDMDT,DSEDNDT
 REAL,DIMENSION(LK)::DIEMD,DIEND!change in mass and number conc
   !due to collision_coalescence and pcpt (bin)
 REAL,DIMENSION(JMINP:JMAXP,KKP,LK)::DSEDMDT,DSEDNDT
  !change in mass and number conc due to sedimentation
! 
!*COMDECK XD 
 COMMON/XD/X_bin,XX,DIAM,BET,ALP
 REAL,DIMENSION(LKDD)::X_bin,XX,DIAM
 REAL,DIMENSION(LK) ::BET,ALP
! 
!*COMDECK RI
  COMMON/RI/XK,DIEMC,DIENC,DIONE,XKmean, XKlocal, &
  XKlayer,XK_gr
  REAL,DIMENSION(LK)::XK,DIEMC,DIENC,DIONE,XKmean,XK_gr
  REAL, DIMENSION(JMINP:JMAXP,KKP,LK) :: XKlocal
  REAL, DIMENSION(KKP,LK) :: XKlayer
! 
!*COMDECK AH
 COMMON/AH/KBAR,ACON
 REAL,DIMENSION(LKCO,LKCO)::ACON
 REAL KBAR(LK,LK)
! 
!*COMDECK LB
 COMMON/LB/PLL,KIJ
 REAL PLL(LKCO,LKCO,LKCO)
 REAL KIJ(LKCO,LKCO)
! 
!*COMDECK SK
 COMMON/SK/SMCK
 REAL SMCK(LKCO)
! 
!*COMDECK LM
 COMMON/LM/PLM
 REAL PLM(LKCO,LKCO)
! 
!*COMDECK S
 COMMON/S/AX,AD
 REAL,DIMENSION(20)::AX,AD
! 
!*COMDECK SSS
 COMMON/SSS/AZ0,AM3,AM4,CA,CB,CC,CD
 REAL,DIMENSION(LK)::AZ0,AM3,AM4
 REAL,DIMENSION(LKCO)::CA,CB,CC,CD
! 
!*COMDECK XY
 COMMON/XY/CM,CN,SCM
 REAL,DIMENSION(LK)::CM,CN,SCM
! 
!*COMDECK ZM
 COMMON/ZM/PSI2,F2 !rm vars F and PSI, see comments in SXY
 REAL,DIMENSION(LK)::PSI2,F2
! 
!*COMDECK A
 COMMON/A/SN0,SN,SM0,SM,SMC
 REAL,DIMENSION(LKCO)::SN0,SN,SM0,SM,SMC
! 
!*COMDECK CONEV
! variables for condensation and evaporation subroutines 
 COMMON/CONEV/AM0,AN0,B,C,DD,EE,ZK,YK
 REAL,DIMENSION(LK)::AM0,AN0,B,C,DD,EE,ZK,YK
!
!*COMDECK RW
! COMMON/RW/RE,DV,SC,VNTF,SC23
 REAL,DIMENSION(KKP)::DV,SC
 REAL,DIMENSION(LK,KKP)::RE
 REAL,DIMENSION(LK,KKP):: VNTF
 REAL :: SC23
!
!*COMDECK COL_COAL
!variables used in SXY,SCONC - I'm not sure what they do
 COMMON/COL_COAL/SM1,SM2,SM3,SM4,SM5,SN1,SN2,SN3,SN4
 REAL :: SM1,SM2,SM3,SM4,SM5,SN1,SN2,SN3,SN4
!*COMDECK SEDIM
 !variables for sedimentation.
 COMMON/SEDIM/AMS,VT,VBAR, QL_SED, QLN_SED
 REAL,DIMENSION(JMINP:JMAXP,KKP):: AMS,VT,VBAR 
 REAL, DIMENSION(JMINP:JMAXP,KKP,NQP) :: QL_SED
 REAL, DIMENSION(JMINP:JMAXP,KKP,NQP) :: QLN_SED!
!*COMDECK NUCFACT
 COMMON/NUCFACT/K_factor, K_rad, K_mass, RNstar, R_CRIT, RN_CRIT &
,k_radorig, radeq
 REAL, DIMENSION(LN2):: K_factor !factor to multiply aerosol by to
 !calculate the size of cloud drop resulting from aerosol at 100%RH
 REAL, DIMENSION(LN2):: K_rad, k_radorig, radeq 
 ! radius of cloud droplet resulting 
 ! from K_factor*rad_aersol
 REAL, DIMENSION(LN2):: K_mass! mass of cloud droplet
 REAL, DIMENSION(LN2):: RNstar
                 !critical radius for nuc (not the same as rstar in kohler)
                 !based on Ivanova et al (1977)   
  REAL, DIMENSION(LN2):: R_CRIT, RN_CRIT
!
!*COMDECK NUCMAIN
 COMMON/NUCMAIN/scrit,ankcc,amkcc,CCINTER
 REAL, DIMENSION(LNUC) :: SCRIT !critical supersaturation for nuc
 REAL, DIMENSION(LN2) :: CCINTER !change aerosol number due to nucleation
 REAL,DIMENSION(LK) :: ANKCC !CDNC resulting from nuc
 REAL,DIMENSION(LK) :: AMKCC !CDmass resulting from nuc 
!
!*COMDECK BIN_RAIN
 COMMON/BIN_RAIN/RATE,ZTOTAL,REFF,ZW_bin
 REAL, DIMENSION(JMINP:JMAXP,KKP) :: RATE,ZTOTAL,REFF 
 REAL, DIMENSION(JMINP:JMAXP,KKP,LK) :: ZW_bin

!*COMDECK REGEN
 COMMON/REGEN/CCNLOSS, CDNCEVAP, CCNREGEN, BINEVAP,ccnorig, &
 TOTCCNORIG, CCNfrac, CCNDIFF, NUCREG_COUNT,LREG, CCNNEWTOT,& 
 CCNORIGTOT, CCNORIGAVG, CCNNEWAVG, CCNxcess, CCNNEWAFT, AEROTOT,&
 qvtot,AN1OLD,AM1OLD
 REAL, DIMENSION(LN2) ::CCNNEWTOT,CCNORIGTOT, CCNORIGAVG, CCNNEWAVG
 REAL, DIMENSION(LN2) :: CCNLOSS, CCNNEWAFT 
 REAL, DIMENSION(JMINP:JMAXP,KKP,Ln2):: CCNREGEN, ccnorig
 REAL, DIMENSION(JMINP:JMAXP,KKP,Ln2):: CCNfrac, CCNDIFF 
 REAL, DIMENSION(JMINP:JMAXP,KKP,Ln2):: NUCREG_COUNT, CCNxcess
 REAL, DIMENSION(JMINP:JMAXP,KKP,LK) :: BINEVAP
 REAL, DIMENSION(JMINP:JMAXP,KKP) ::CDNCEVAP,AN1OLD,TOTCCNORIG,AM1OLD
 REAL :: AEROTOT, qvtot
 INTEGER :: LREG
!*COMDECK AERO_DIAG
 COMMON /AERO_DIAG/ CCN_NUC, CCN_REG, dqn_act, dqn_reg,ssat,dcrit 
 REAL, DIMENSION(JMINP:JMAXP,KKP,LN2) :: CCN_NUC,CCN_REG
 REAL, DIMENSION(JMINP:JMAXP,KKP) :: dqn_act, dqn_reg
 REAL, DIMENSION(JMINP:JMAXP,KKP) :: SSAT,dcrit 
!*COMDECK CCN_TEST
 COMMON /CCN_TEST/adv_nd,adv_w_nd,adv_v_nd,ccn_init
 REAL, DIMENSION(JMINP:JMAXP, KKP) ::adv_nd,adv_w_nd,adv_v_nd, ccn_init
!*COMDECK rand_test
 common /rand_test/randarr
! real randarr(jjp,kkp,iip)
 real randarr(JMINP:JMAXP*IUSETHP,(KKP-1)*IUSETHP+1,& 
(IDIMMINP-1)*IUSETHP+1:(IDIMMAXP-1)*IUSETHP+1)
!
!*COMDECK RAD_TIME
 COMMON /RAD_TIME/zradtime,zdradtime,ziradswitch
 REAL :: zradtime,zdradtime
 INTEGER :: ziradswitch
!
!*COMDECK JR2006
 COMMON /JR2006/t_mix,t_mix_approx,tevap_bb,tevap_squires,da_no &
,R_bar, R_int,da_no_rat,evap_a
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: t_mix  
                             !mixing timescale calced using eddy dis
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: t_mix_approx 
                             !mixing timescale calced using vert vel only
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: tevap_bb
                             !evaporation timescale calced using method
                             !based on brenguier+burnet
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: tevap_squires
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: da_no   
                             !dahmkohler number,ratio of t_mix to t_evap
 REAL,DIMENSION(JMINP:JMAXP,KKP) ::  R_bar   
                             !mean cloud drop radius
 REAL,DIMENSION(JMINP:JMAXP,KKP) ::  R_int   
                             !integral radius
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: da_no_rat, evap_a
!
!*COMDECK MICRO_MIX
!MICRO_MIX contains the subgrid and resolved tendency in lwc
!and nd due to cond/evap, i.e. this does not include activation
!so sum will not equal dq_dt(IQL) or dq_dt(ND) in the present
!diags
 COMMON /MICRO_MIX/ce_dqdt_sub,ce_dqdt_res,am1_diag,an1_diag &
 ,am2_diag,an2_diag 
 REAL,DIMENSION(JMINP:JMAXP,KKP,NQP) :: ce_dqdt_sub,ce_dqdt_res
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: am1_diag,an1_diag
 REAL,DIMENSION(JMINP:JMAXP,KKP) :: am2_diag,an2_diag 
 
  common /auto_con_diag/rmass_cw, rnum_cw, auto, acw, revp, rcond
 REAL,DIMENSION(KKP) :: rmass_cw, rnum_cw 
 REAL, DIMENSION(5,KKP) :: auto, acw, revp, rcond


end module mphys_tau_bin_declare
