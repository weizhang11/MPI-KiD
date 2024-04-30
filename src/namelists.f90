! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to deal with namelists for input
!
module namelists

  Use parameters
  Use header_data, only : mphys_id
  Use switches

#if SHIPWAY_MICRO == 1
  ! Temporary for adding in 4a switches
  Use mphys_switches, only: iopt_act, option, aerosol_option        &
     , l_aaut, l_aacc, l_aevp, l_ased, l_warm                       &
     , l_inuc, iopt_rcrit, iopt_inuc                                &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , l_sed_3mdiff, l_cons, l_abelshipway, l_condensation          &
     , l_active_inarg2000, l_separate_rain,                         &
     l_g, l_sg , l_pcond, l_praut, l_pracw  , l_pracr  , l_prevp  , & 
     l_psedl  , l_psedr  , l_ptidy  , l_ptidy2 , l_pinuc , l_pidep, & 
     l_piacw  , l_psaut  , l_psdep  , l_psacw  , l_pgdep , l_pseds, & 
     l_psedi  , l_psedg  , l_psaci  , l_praci  , l_psacr , l_pgacr, & 
     l_pgacw  , l_pgaci  , l_pgacs  , l_piagg  , l_psagg , l_pgagg, & 
     l_psbrk  , l_pgshd  , l_pihal  , l_psmlt  , l_pgmlt , l_phomr, & 
     l_phomc  , l_pssub  , l_pgsub  , l_pisub  , l_pimlt            &
     , max_step_length, max_sed_length, process_level,   &
     l_subseds_maxv, cfl_vt_max
  Use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3
#endif

  implicit none
  
  namelist/mphys/num_h_moments, num_h_bins, mom_init, &
       h_names, mom_names, mom_units,num_aero_moments,num_aero_bins, &
       aero_mom_init, aero_N_init, aero_sig_init, aero_rd_init, aero_names
  
  namelist/control/dt, dg_dt, mphys_scheme, mphys_var &
       , wctrl, zctrl, tctrl, pctrl_z, pctrl_v, pctrl_T, ipctrl &
       , xctrl, lhf_ctrl, shf_ctrl, diaglevel, dgstart &
       , init_hydrometeors

  namelist/case/input_file, l_input_file, ifiletype, icase

  namelist/switch/l_mphys, l_advect, l_diverge, l_pupdate &
       , l_fix_qv, l_nomphys_qv, l_noadv_qv, l_posadv_qv &
       , l_fix_theta, l_nomphys_theta, l_noadv_theta  &
       , l_noadv_hydrometeors, l_nodiv_hydrometeors, l_sediment &
       , isurface, l_noadv_aerosols, l_nodiv_aerosols, l_fix_aerosols &
       , l_diverge_advection, l_periodic_bound  &
       , l_force_positive, l_raut, l_constant_density &
       ! TAU bin switches
       , l_act, l_cond_evap, l_coll_coal, l_break, l_fix_supersat &
       , l_dist_activated_drops, mu_act, act_bin_ind

  logical :: iiwarm=.false.
  character(200) :: KiD_outdir=''
  character(200) :: KiD_outfile=''
  real :: set_Nc = 100.  ! for setting cloud number concentration (cm-3)
  real :: amp_fact= 1.0 ! factor applied to the amplitude of the Cu velocity
                        ! field to obtain different clouds. Default is
                        ! the standard case  	 	 
  real :: cltop=2700.   ! cloud top for the standard Cu case  	 	 
  real :: no_precip_time = 0.0 ! time that there is no precip process or sedimentation
  real :: smax = 999.0 ! maximum supersaturation that can be used in cond and act
                       ! set to 999 % as default so that this is not used unless 
                       ! smax is set in user namelist
  real :: smax_limit_time = 0.0
  namelist/addcontrol/iiwarm, KiD_outdir, KiD_outfile  &
#if SHIPWAY_MICRO == 1
       ! Shipway 4A ...
     , option, l_evap, l_sed_3mdiff &
     , max_step_length, max_sed_length, diag_mu, l_rsource, l_raut &
     , l_subseds_maxv &
     , cfl_vt_max &
     , l_evaporation, l_rain, l_sed, l_boussinesq, diag_mu_option   &
     , p1, p2, p3 &
     , sp1, sp2, sp3 &
     , l_abelshipway, l_cons &
     , l_coll_coal, l_cond, l_condensation, iopt_act &
     , aerosol_option, l_aaut, l_aacc, l_aevp, l_ased &
     , l_warm, l_inuc, iopt_rcrit, process_level   &
     , l_active_inarg2000, iopt_inuc, l_cu_cold &
     , l_separate_rain,         &
     ! Add in all switches so that all testing can be done through namelist
     l_g      , & ! switch for graupel
     l_sg     , & ! switch for snow and graupel
     l_pcond  , & ! Condensation
     l_praut  , & ! Autoconversion cloud -> rain
     l_pracw  , & ! Accretion  cloud -> rain
     l_pracr  , & ! aggregation of rain drops
     l_prevp  , & ! evaporation of rain
     l_psedl  , & ! sedimentation of cloud
     l_psedr  , & ! sedimentation of rain
     l_ptidy  , & ! tidying term 1
     l_ptidy2 , & ! tidying term 2
     l_pinuc  , & ! ice nucleation
     l_pidep  , & ! ice deposition
     l_piacw  , & ! ice accreting water
     l_psaut  , & ! ice autoconversion ice -> snow
     l_psdep  , & ! vapour deposition onto snow
     l_psacw  , & ! snow accreting water
     l_pgdep  , & ! vapour deposition onto graupel
     l_pseds  , & ! snow sedimentation
     l_psedi  , & ! ice sedimentation
     l_psedg  , & ! graupel sedimentation
     l_psaci  , & ! snow accreting ice
     l_praci  , & ! rain accreting ice
     l_psacr  , & ! snow accreting rain
     l_pgacr  , & ! graupel accreting rain
     l_pgacw  , & ! graupel accreting cloud water
     l_pgaci  , & ! graupel accreting ice
     l_pgacs  , & ! graupel accreting snow
     l_piagg  , & ! aggregation of ice particles
     l_psagg  , & ! aggregation of snow particles
     l_pgagg  , & ! aggregation of graupel particles
     l_psbrk  , & ! break up of snow flakes
     l_pgshd  , & ! shedding of liquid from graupel
     l_pihal  , & ! hallet mossop
     l_psmlt  , & ! snow melting
     l_pgmlt  , & ! graupel melting
     l_phomr  , & ! homogeneous freezing of rain
     l_phomc  , & ! homogeneous freezing of cloud droplets
     l_pssub  , & ! sublimation of snow
     l_pgsub  , & ! sublimation of graupel
     l_pisub  , & ! sublimation of ice
     l_pimlt    & ! ice melting
#endif
     ! Thompson 09...
     , l_reuse_thompson_lookup &
     ! Flags for KiD-A intercomparison
     , set_Nc, amp_fact, cltop, no_precip_time, smax, smax_limit_time  
  

  ! Namelist input...

  character(200) :: fileName=''
  character(200) :: fileNameIn=''
  character(200) :: fileNameOut=''
  character(200) :: namelistIn='namelists/input.nml'
  character(200) :: fileIn=''
  character(200) :: fileOut=''
  logical :: fexist

  namelist/namelistToUse/fileIn, fileOut

contains

  subroutine read_namelist
    !
    ! Read the namelists from file
    !

#if COMMANDLINE == 1
    ! This bit is F2003 - If your compiler doesnt support 
    ! it you need to comment the line out you can then specify 
    ! which namelist to use throught namelists/input.nml file 
    ! or else use command line processing that works with your 
    ! compiler (nearly all will do it but not in a portable way).
    write(*,*) 'Querying command line'
    CALL GET_COMMAND_ARGUMENT(1,fileNameIn)
    CALL GET_COMMAND_ARGUMENT(2,fileNameOut)
#endif

    if (trim(fileNameIn)=='')then  ! Not input at the command line 
                                   ! so use input.nml
#ifdef DEF_CASE
      write(namelistIn, '(A,A,A)') 'namelists/', DEF_CASE, '_input.nml'
#endif

      write(*,*) 'Unable to determine input file from command line, so querying ', trim(namelistIn), ' instead...'
      inquire(file=namelistIn, exist=fexist)
      if (fexist) then
        open(2, file=namelistIn)
        read(2, namelistToUse)
        close(2)
        write(*, namelistToUse)
      end if
      fileNameIn  = fileIn
      if (trim(fileOut)/='')fileNameOut = fileOut
    end if
    
    if (trim(fileNameIn)/='')fileName=fileNameIn

    write(6,*) 'Using namelist: ', trim(fileName)

    open(1, file=fileName)
!    rewind(1)
    read(1,mphys) 
!    rewind(1)
    read(1,case)
!    rewind(1)
    read(1,control) 
!    rewind(1)
    read(1,switch) 
!    rewind(1)
    read(1,addcontrol) 
    close(1)

    select case(mphys_scheme)
    case('lem2.4')
       imphys=imphys_lem2_4
       mphys_id='LEM2.4'
    case('tau_bin')
       imphys=imphys_tau_bin
       mphys_id='TAU_bin'
    case('thompson') ! NB same as thompson09
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson09')
       imphys=imphys_thompson09
       mphys_id='thompson09'
    case('thompson06')
       imphys=imphys_thompson06
       mphys_id='thompson06'
    case('thompson07')
       imphys=imphys_thompson07
       mphys_id='thompson07' 
    case('morr_two_moment')
       imphys=imphys_morr_two_moment
       mphys_id='morr_two_moment'
    case('um7_3')
       imphys=imphys_um7_3
       mphys_id='um7_3'
    case('wsm6')
       imphys=imphys_wsm6
       mphys_id='wsm6'
    case('wdm6')
       imphys=imphys_wdm6
       mphys_id='wdm6'
    case('casim')
       imphys=imphys_casim
       mphys_id='casim'
    case default
       print*, 'Mphys scheme not recognized: ', mphys_scheme
       print*, 'Did you mean:' 
       print*, '   lem2.4?'
       print*, '   um7_3?'
       print*, '   tau_bin?'
       print*, '   thompson09?'
       print*, '   thompson07?'
       print*, '   morr_two_moment?'
       print*, '   wsm6?'
       print*, '   4A?'
       print*, '(NB not all available in release version)'
       stop
    end select

    if (trim(input_file)=='')l_input_file=.False.

    if (.not. l_input_file)ifiletype=itest_case

  end subroutine read_namelist

end module namelists

