! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Driver for 1D Kinematic Driver (KiD) model 
!
! Author: Ben Shipway
!
! For version details see header_data.f90
!

Module main
 
  Use typeKind
  Use parameters, only : dt, dg_dt, nx
  Use namelists, only : read_namelist
  Use runtime, only : time, time_step, n_times
  Use switches
  Use set_profiles,  only : read_profiles
  Use interpolation, only : interpolate_input, interpolate_forcing
  Use diagnostics,   only : save_diagnostics_1d, save_diagnostics_2d, &
       write_diagnostics, query_dgstep, l_dgstep
  Use derived_fields, only : calc_derived_fields
  Use advection_interface, only : advect_column
  Use mphys_interface, only : mphys_column
  Use stepfields, only : step_column
  Use divergence, only : diverge_column

  use mphys_parameters, only : shipway_time, micro_time

#ifdef USE_MPI
  Use mpi_mod
#endif

  Implicit none

  real(wp) :: micro_time_max, micro_time_min, micro_time_avg
  real(wp) :: shipway_time_max, shipway_time_min, shipway_time_avg



contains

  subroutine main_loop
    
    integer :: itime      ! loop counter for time
    real :: stime,etime
    !
    ! Start by reading in namelists
    !

#ifdef USE_MPI
    call KiD_mpi_init
#endif
    if (l_namelists) call read_namelist

    ! Set up the initial fields and forcing
    if (l_input_file)then
       call read_profiles(input_file)
    else
       call read_profiles(icase)
    end if

!  call messages('interp')
    call interpolate_input(ifiletype)

!  call messages('forcing')
    call interpolate_forcing

!  call messages('derive')
    call calc_derived_fields

    ! Do we want to do diagnostics on this timestep?
    call query_dgstep

    if ( nx == 1 ) then 
       call save_diagnostics_1d
    else 
!  call messages('save')
       call save_diagnostics_2d
    endif

!   call cpu_time(stime)
    do itime=1,n_times

       time=time+dt
       time_step=time_step+1
  
       ! Do we want to do diagnostics on this timestep?
!  call messages('qstep')
       call query_dgstep

!  call messages('forcing')
       call interpolate_forcing

!  call messages('derive')
       call calc_derived_fields

       if (l_advect)then
!  call messages('advect')
          call advect_column(scheme_id=0)
       end if


       if (l_diverge)then
!  call messages('diverge')
          call diverge_column
       end if

       if (l_mphys)then
!  call messages('mphys')
          call mphys_column(scheme_id=imphys)
       end if

!  call messages('step')
       call step_column
       
       if ( nx == 1 ) then
          call save_diagnostics_1d
       else
!  call messages('save')
          call save_diagnostics_2d
       endif

    end do

#ifdef USE_MPI
   call mpi_allreduce(micro_time, micro_time_max,1,mpi_double, &
                      mpi_max,mpi_comm_world,mpi_ierr)
   call mpi_allreduce(micro_time, micro_time_min,1,mpi_double, &
                      mpi_min,mpi_comm_world,mpi_ierr)
   call mpi_allreduce(micro_time, micro_time_avg,1,mpi_double, &
                      mpi_sum,mpi_comm_world,mpi_ierr)

   call mpi_allreduce(shipway_time, shipway_time_max,1,mpi_double, &
                      mpi_max,mpi_comm_world,mpi_ierr)
   call mpi_allreduce(shipway_time, shipway_time_min,1,mpi_double, &
                      mpi_min,mpi_comm_world,mpi_ierr)
   call mpi_allreduce(shipway_time, shipway_time_avg,1,mpi_double, &
                      mpi_sum,mpi_comm_world,mpi_ierr)


   if ( rank == master ) then
      micro_time_avg = micro_time_avg / ncpu
      shipway_time_avg = shipway_time_avg / ncpu

      print*,'CASIM shipway_microphysics time MAX : ',shipway_time_max
      print*,'CASIM shipway_microphysics time MIN : ',shipway_time_min
      print*,'CASIM shipway_microphysics AVG : ',shipway_time_avg

      print*, '              '
      print*,'CASIM microphysics_common time MAX : ',micro_time_max
      print*,'CASIM microphysics_common time MIN : ',micro_time_min
      print*,'CASIM microphysics_common AVG : ',micro_time_avg



   endif

#endif

#ifdef USE_MPI
    call KiD_mpi_barrier
    if ( rank == master ) then
#endif
    if (l_write_dgs) call write_diagnostics
#ifdef USE_MPI
    endif
    call KiD_mpi_finalize
#endif


  end subroutine main_loop

End Module main
