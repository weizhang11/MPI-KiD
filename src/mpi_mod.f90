! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to interface with choice of microphysics schemes
!

module mpi_mod

  Use typeKind
  Use parameters, only : nz, pnx, nx, nspecies, naerosol, unset_real,max_nbins
  Use class_species, only :  species
  Use column_variables

  Implicit none

  integer,parameter :: master = 0, mpi_tag0=0, mpi_tag1=1
  integer,allocatable :: pidx(:),sloc(:),eloc(:),subnxall(:)
  integer,allocatable :: ksloc(:),keloc(:)
  integer :: sproc(2),rproc(2)
  integer :: isx,iex
  integer :: mpi_ierr,nworker,ncpu,rank
  integer :: mpi_4reqs(4),mpi_2reqs(2)
  real(wp) :: xs,xe
  real(sp) :: savescalar_sp
  real(dp) :: savescalar_dp
  real(sp),dimension(:),allocatable :: save1dfield_sp,field1d_sp
  real(dp),dimension(:),allocatable :: save1dfield_dp,field1d_dp
  real(sp),dimension(:,:),allocatable :: save2dfield_sp,field2d_sp
  real(dp),dimension(:,:),allocatable :: save2dfield_dp,field2d_dp
  real(sp),dimension(:,:,:),allocatable :: save3dfield_sp
  real(dp),dimension(:,:,:),allocatable :: save3dfield_dp
  real(wp),dimension(:,:),allocatable :: gfield_mask
  real(wp),dimension(:),allocatable :: gx

#ifdef USE_MPI

  include 'mpif.h'

contains

  subroutine KiD_mpi_init
    integer :: i,j,ip, itmp1,itmp2,itmp3
    real(wp) :: tmp1,tmp2,tmp3
    real(wp),dimension(:,:),allocatable :: pvar1,pvar2,gvar1

    ! MPI initializations
    call mpi_init(mpi_ierr)
    call mpi_comm_rank(mpi_comm_world,rank,mpi_ierr)
    call mpi_comm_size(mpi_comm_world,ncpu,mpi_ierr)

    nworker = ncpu - 1

    ! Define x-subdomain ---
    allocate(subnxall(0:nworker))
    itmp1 = int(pnx/ncpu)
    itmp2 = itmp1 * ncpu
    itmp3 = pnx - itmp2

    subnxall(:) = itmp1

    do i = 0,itmp3-1
      subnxall(i) = subnxall(i)+1 
    end do

    nx = subnxall(rank) ! nx for subdomain

    if ( rank == 0 ) then
       print*, 'Initializing KiD-MPI'
       print*, 'Total number of processors= ',ncpu
    endif

! Check Min. nx > 1
    if ( minval(subnxall) <= 1 ) then
       if ( rank == 0 ) then
          print*, 'nx for each processor is lower than 1.'
          print*, 'Try to decrease the number of processors'
          print*, 'ABORT: KiD-parallel'
       endif
       call KiD_mpi_barrier
       call KiD_mpi_finalize
    endif

    ! Subdomain grid indexing as a global index
    allocate(pidx(0:nx+1))
    itmp1 = 0 
    do i = 0,rank-1
      itmp1 = itmp1 + subnxall(i)
    end do
    
    do i = 1,nx
      pidx(i) = itmp1 + i
    end do 

    pidx(0)    = pidx(1)-1
    pidx(nx+1) = pidx(nx)+1


    if ( rank == 0 ) then
       allocate(sloc(0:nworker),eloc(0:nworker))

       sloc(0) = 1 ; eloc(0) = nx
       do i = 1,nworker
          sloc(i) = eloc(i-1) + 1
          eloc(i) = sloc(i)   + subnxall(i) - 1 
       end do
    endif ! master
   
    if ( rank /= master ) then
       deallocate(subnxall)
    endif

    ! Allocations of column_variables ---------------------------------------
    ! fields carried on each grid level
    allocate(theta_ref(nz))     ! Smoothed initial theta profile used as
                                ! a reference profile

    allocate(                   &
          z(nz)                 &   ! height (m)
         ,x(0:nx+1)            &   ! horizontal dimension (m)
         ,exner(nz, 0:nx+1)    &   ! exner pressure
         ,w(nz, 0:nx+1)        &   ! vertical velocity (m/s)
         ,v(nz, 0:nx+1)        &   ! horizontal velocity (m/s)
         ,rho(nz)               &   ! density (kg/m3)
         ,theta(nz, 0:nx+1)    &   ! potential temperature(K)
         ,qv(nz, 0:nx+1)       &   ! water vapour mixing ratio (kg/kg)
         ,ss(nz, 0:nx+1)       )   ! supersaturation (transported as a
                                    ! passive scalar)

          w(:,:) = 0
          v(:,:) = 0
          ss(:,:) = 0

    allocate(                  &
          z_half(nz)           &  ! height on half levels (m)
         ,w_half(nz,0:nx+1)   &  ! vertical velocity on half levels(m/)
         ,rho_half(nz)         &  ! density on half levels(kg/m3)
         ,x_half(0:nx+1)      &  ! horizontal half levels (m)
         ,v_half(nz,0:nx+1)   )  ! horizontal winds on half horiz. grid (m/s)

          w_half(:,:) = 0 
          rho_half(:) = 0

    !grid spacing
    allocate(               &
          dz(nz)            &   ! dz (m)
         ,dz_half(nz)       &   ! dz on half levels (m)
         ,dx(0:nx+1)       &   ! dx (m)
         ,dx_half(0:nx+1)  )

    allocate(                  &
          pmb(nz, 0:nx+1)     & ! Pressure in mb
         ,TdegK(nz, 0:nx+1)   ) ! Temperature in K

          pmb(:,:) = unset_real
          TdegK(:,:) = unset_real

    allocate(                                  &
          hydrometeors(nz, 0:nx+1, nspecies)  &  ! hydrometeors
         ,aerosol(nz, 0:nx+1, naerosol)       )  ! aerosols

    ! Initial profiles
    allocate(                   &
          thinit(nz, 0:nx+1)   &   ! initial potential temperature
         ,qvinit(nz, 0:nx+1)   )   ! initial water vapour

    ! Applied forcing tendencies (if used)
    allocate(                   &
          Tforce(nz,0:nx+1)    &   ! (Horizontal) Advective forcing:
                                    ! temperature(K)
         ,qforce(nz,0:nx+1)    )   ! (Horizontal) Advective forcing:
                                    ! vapour(kg/kg)

          Tforce(:,:) = 0
          qforce(:,:) = 0

    allocate(Trelax(nz))   ! relaxation timescale (if using relaxation)
          Trelax(:) = 1.

    ! Tendencies calculated by the model (advection/microphysics)
    allocate(                       &
          dtheta_adv(nz, 0:nx+1)   &  ! Advective tendency: Theta (K/s)
         ,dqv_adv(nz, 0:nx+1)      &  ! Advective tendency: qv (kg/kg/s)
         ,dss_adv(nz, 0:nx+1)      )  ! Advective super sat tendency for bin micro

    allocate(                                      &
          dhydrometeors_adv(nz, 0:nx+1, nspecies) & ! Advective tendency:
                                                     ! hydrometeors (?/s)
         ,daerosol_adv(nz, 0:nx+1, naerosol)      ) ! Advective tendency:
                                                     ! aerosol (?/s)

    allocate(                         &
          dTheta_mphys(nz, 0:nx+1)   &   ! Mphys tendency: Theta (K/s)
         ,dqv_mphys(nz, 0:nx+1)      &   ! Mphys tendency: qv (kg/kg/s)
         ,dss_mphys(nz, 0:nx+1)      )   ! Mphys tendency: ss

    allocate(                         &
          dhydrometeors_mphys(nz, 0:nx+1, nspecies) & ! Mphys tendency:
                                                       ! hydrometeors (?/s)
         ,daerosol_mphys(nz, 0:nx+1, naerosol)      ) ! Mphys tendency:
                                                       ! aerosol (?/s)

    allocate(                       &
          dTheta_div(nz, 0:nx+1)   &   ! Divergence tendency: Theta (K/s)
         ,dqv_div(nz, 0:nx+1)      &   ! Divergence tendency: qv (kg/kg/s)
         ,dss_div(nz, 0:nx+1)      )   ! Divergence tendency: ss

    allocate(                         &
          dhydrometeors_div(nz, 0:nx+1, nspecies) & ! Divergence tendency:
                                                     ! hydrometeors (?/s)
         ,daerosol_div(nz, 0:nx+1, naerosol)      ) ! Divergence tendency:
                                                     ! aerosol (?/s)

    ! Some test cases use initial values for hydrometeors as follows
    allocate(                         &
          hydrometeors_init(nz,nspecies) & ! Initial values for
                                           ! hydrometeors
         ,aerosol_init(nz,naerosol)      ) ! Initial values for
                                           ! aerosol (?/s)

! NB Forcing of variables not yet full functional
    allocate(                         &
          dhydrometeors_force(nz,0:nx+1,nspecies) & ! Imposed forcing:
                                                     ! hydrometeors (?/s)
         ,daerosol_force(nz,0:nx+1,naerosol)      ) ! Imposed forcing:
                                                     ! aerosol (?/s)
    ! scalars
    allocate(                  &
          wth_surf(0:nx+1)    &   ! surface potential temperature flux (Km/s)
         ,wqv_surf(0:nx+1)    )   ! surface vapour flux (kg/kg m/s)

          wth_surf(:) = 0
          wqv_surf(:) = 0

    ! Mask for applying increments. Increments are not applied where
    ! field is zero (and are multiplied by weight otherwise).

    allocate(field_mask(nz,0:nx+1))  ! Mask for increments

          field_mask(:,:) = 1.

    ! MPI SEND-RECV tag setup ---

    ! SEND & RECV tags --------------------------------------------------!
    !
    ! tag 1 = Left CPU ; tag 2 = Right CPU
    !
    !--------------------------------------------------------------------!
    !
    !            sproc(1)                              sproc(2)
    !
    !                       SEND              SEND
    !           [        ] <------ [      ] ------> [        ]
    !           [Proc n-1]         [Proc n]         [Proc n+1]
    !           [        ] ------> [      ] <------ [        ]
    !                       RECV              RECV
    !
    !            rproc(1)                              rproc(2)
    !
    !---------------------------------------------------------------------

    if ( rank == 0 ) then
      sproc(1) = nworker
      sproc(2) = 1
      rproc(1) = nworker
      rproc(2) = 1
    elseif ( rank == nworker ) then
      sproc(1) = nworker-1
      sproc(2) = 0
      rproc(1) = nworker-1
      rproc(2) = 0
    else
      sproc(1) = rank-1
      sproc(2) = rank+1
      rproc(1) = rank-1
      rproc(2) = rank+1
    endif

  end subroutine KiD_mpi_init


  subroutine KiD_thompson_k_decomp(ntot, ks, ke)
    integer,intent(in) :: ntot
    integer,intent(out) :: ks,ke
    integer,dimension(:),allocatable :: subnkall,pkidx
    integer :: i,j,k,ip,nk,itmp1,itmp2,itmp3

    if ( .not. allocated(ksloc)) then

       ! Define x-subdomain ---
       allocate(subnkall(0:nworker))
       itmp1 = int(ntot/ncpu)
       itmp2 = itmp1 * ncpu
       itmp3 = ntot - itmp2 
   
       subnkall(:) = itmp1
   
       do i = 1,itmp3
         subnkall(i) = subnkall(i)+1 
       end do
   
       nk = subnkall(rank) ! nx for subdomain
   
       allocate(ksloc(0:nworker),keloc(0:nworker))
   
       ksloc(0) = 1 ; keloc(0) = subnkall(0) 
       do i = 1,nworker
          ksloc(i) = keloc(i-1) + 1
          keloc(i) = ksloc(i)   + subnkall(i) -1
       end do
   
       do i = 0,nworker
          ksloc(i) = ksloc(i) - 1
          keloc(i) = keloc(i) - 1
       end do
       ks = ksloc(rank) ; ke = keloc(rank)
     
       deallocate(subnkall)
    else

       ks = ksloc(rank) ; ke = keloc(rank)
     
    endif

  end subroutine KiD_thompson_k_decomp

!----------------------------------------------------------------------!
! MPI COMM routine
!
! - Update buffer for a given 'signle-precision real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_real_sp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(sp),dimension(nz,0:nx+1) :: var
!----------------------------------------------------------------------!
! Communication: update buffer

      CALL MPI_SENDRECV(var(:,1),nz,                      &
                        MPI_REAL4,sproc(1),0,              &
                        var(:,nx+1),nz,                  &
                        MPI_REAL4,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
      CALL MPI_SENDRECV(var(:,nx),nz,                  &
                        MPI_REAL4,sproc(2),1,              &
                        var(:,0),nz,                      &
                        MPI_REAL4,rproc(1),1,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_real_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM routine
!
! - Update buffer for a given 'double-precision real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_real_dp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(wp),dimension(nz,0:nx+1) :: var
!----------------------------------------------------------------------!
! Communication: update buffer


!     CALL MPI_SENDRECV(var(:,1),nz,                      &
!                       MPI_REAL8,sproc(1),0,              &
!                       var(:,nx+1),nz,                  &
!                       MPI_REAL8,rproc(2),0,              &
!                       MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
!     CALL MPI_SENDRECV(var(:,nx),nz,                  &
!                       MPI_REAL8,sproc(2),1,              &
!                       var(:,0),nz,                      &
!                       MPI_REAL8,rproc(1),1,              &
!                       MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

      CALL MPI_IRECV(var(:,nx+1),nz,                             &
                     MPI_REAL8,rproc(2),mpi_tag0,MPI_COMM_WORLD, &
                     mpi_4reqs(1),mpi_ierr )
      CALL MPI_IRECV(var(:,0),nz,                                &
                     MPI_REAL8,rproc(1),mpi_tag1,MPI_COMM_WORLD, &
                     mpi_4reqs(2),mpi_ierr )

      CALL MPI_ISEND(var(:,nx),nz,                               &
                     MPI_REAL8,sproc(2),mpi_tag1,mpi_comm_world, &
                     mpi_4reqs(3),mpi_ierr)
      CALL MPI_ISEND(var(:,1),nz,                                &
                     MPI_REAL8,sproc(1),mpi_tag0,mpi_comm_world, &
                     mpi_4reqs(4),mpi_ierr)

      CALL MPI_WAITALL(4,mpi_4reqs,MPI_STATUS_IGNORE,mpi_ierr)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_real_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM routine
!
! - Update right buffer for a given 'sigle-precision real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_right_buffer_real_sp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(sp),dimension(nz,0:nx+1) :: var
!----------------------------------------------------------------------!
! Communication: update buffer

      CALL MPI_SENDRECV(var(:,1),nz,                      &
                        MPI_REAL4,sproc(1),0,              &
                        var(:,nx+1),nz,                  &
                        MPI_REAL4,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_right_buffer_real_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM routine
!
! - Update buffer for a given 'reversed dimension dp-real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_real_reverse_dp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(wp),dimension(0:nx+1,nz) :: var
!----------------------------------------------------------------------!
! Communication: update buffer

      CALL MPI_SENDRECV(var(1,:),nz,                      &
                        MPI_REAL8,sproc(1),0,              &
                        var(nx+1,:),nz,                  &
                        MPI_REAL8,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
      CALL MPI_SENDRECV(var(nx,:),nz,                  &
                        MPI_REAL8,sproc(2),1,              &
                        var(0,:),nz,                      &
                        MPI_REAL8,rproc(1),1,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_real_reverse_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM rouine
!
! - Update buffer for a given 'reversed dimension sp-real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_real_reverse_sp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(sp),dimension(0:nx+1,nz) :: var
!----------------------------------------------------------------------!
! Communication: update buffer

      CALL MPI_SENDRECV(var(1,:),nz,                      &
                        MPI_REAL4,sproc(1),0,              &
                        var(nx+1,:),nz,                  &
                        MPI_REAL4,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
      CALL MPI_SENDRECV(var(nx,:),nz,                  &
                        MPI_REAL4,sproc(2),1,              &
                        var(0,:),nz,                      &
                        MPI_REAL4,rproc(1),1,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!     CALL MPI_IRECV(var(nx+1,:),nz,                             &
!                    MPI_REAL4,rproc(2),mpi_tag0,MPI_COMM_WORLD, &
!                    mpi_4reqs(1),mpi_ierr )
!     CALL MPI_IRECV(var(0,:),nz,                                &
!                    MPI_REAL4,rproc(1),mpi_tag1,MPI_COMM_WORLD, &
!                    mpi_4reqs(2),mpi_ierr )

!     CALL MPI_ISEND(var(nx,:),nz,                               &
!                    MPI_REAL4,sproc(2),mpi_tag1,mpi_comm_world, &
!                    mpi_4reqs(3),mpi_ierr)
!     CALL MPI_ISEND(var(1,:),nz,                                &
!                    MPI_REAL4,sproc(1),mpi_tag0,mpi_comm_world, &
!                    mpi_4reqs(4),mpi_ierr)

!     CALL MPI_WAITALL(4,mpi_4reqs,MPI_STATUS_IGNORE,mpi_ierr)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_real_reverse_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM rouine
!
! - Update right buffer for a given 'reversed dimension sp-real' variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_right_buffer_real_reverse_sp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(sp),dimension(0:nx+1,nz) :: var
!----------------------------------------------------------------------!
! Communication: update buffer

      CALL MPI_SENDRECV(var(1,:),nz,                      &
                        MPI_REAL4,sproc(1),0,              &
                        var(nx+1,:),nz,                  &
                        MPI_REAL4,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!     CALL MPI_ISEND(var(1,:),nz,                                &
!                    MPI_REAL4,sproc(1),mpi_tag0,mpi_comm_world, &
!                    mpi_2reqs(2),mpi_ierr)
!     CALL MPI_IRECV(var(nx+1,:),nz,                             &
!                    MPI_REAL4,rproc(2),mpi_tag0,MPI_COMM_WORLD, &
!                    mpi_2reqs(1),mpi_ierr )
!     CALL MPI_WAITALL(2,mpi_2reqs,MPI_STATUS_IGNORE,mpi_ierr)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_right_buffer_real_reverse_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM rouine
!
! - Update buffer for a given 'reversed dimension sp-real' variable
! - No update for x-domain boundaries
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_nobdy_real_reverse_sp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(sp),dimension(0:nx+1,nz) :: var
      real(sp),dimension(nz) :: zval,svar1,rvar1,svar2,rvar2
!----------------------------------------------------------------------!
! Communication: update buffer

      if ( rank == 0 ) then
         zval(:) = var(0,:)
      elseif ( rank == nworker ) then
         zval(:) = var(nx+1,:)
      endif

      CALL MPI_SENDRECV(var(1,:),nz,                      &
                        MPI_REAL4,sproc(1),0,              &
                        var(nx+1,:),nz,                  &
                        MPI_REAL4,rproc(2),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
      CALL MPI_SENDRECV(var(nx,:),nz,                  &
                        MPI_REAL4,sproc(2),1,              &
                        var(0,:),nz,                      &
                        MPI_REAL4,rproc(1),1,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

!     CALL MPI_ISEND(var(nx,:),nz,                               &
!                    MPI_REAL4,sproc(2),mpi_tag1,mpi_comm_world, &
!                    mpi_4reqs(3),mpi_ierr)
!     CALL MPI_ISEND(var(1,:),nz,                                &
!                    MPI_REAL4,sproc(1),mpi_tag0,mpi_comm_world, &
!                    mpi_4reqs(4),mpi_ierr)
!     CALL MPI_IRECV(var(nx+1,:),nz,                             &
!                    MPI_REAL4,rproc(2),mpi_tag0,MPI_COMM_WORLD, &
!                    mpi_4reqs(1),mpi_ierr )
!     CALL MPI_IRECV(var(0,:),nz,                                &
!                    MPI_REAL4,rproc(1),mpi_tag1,MPI_COMM_WORLD, &
!                    mpi_4reqs(2),mpi_ierr )

!     CALL MPI_WAITALL(4,mpi_4reqs,MPI_STATUS_IGNORE,mpi_ierr)

      if ( rank == 0 ) then
         var(0,:) = zval(:)
      elseif ( rank == nworker ) then
         var(nx+1,:) = zval(:)
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_nobdy_real_reverse_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM routine
!
! - Update buffer for a given 'dp-real' variable
! - No update for x-domain boundaries
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_buffer_nobdy_real_dp (var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      real(wp),dimension(nz,0:nx+1) :: var
      real(wp),dimension(nz) :: zval
!----------------------------------------------------------------------!
! Communication: update buffer


      if ( rank == 0 ) then
         zval(:) = var(:,0)
      elseif ( rank == nworker ) then
         zval(:) = var(:,nx+1)
      endif


      CALL MPI_IRECV(var(:,nx+1),nz,                             &
                     MPI_REAL8,rproc(2),mpi_tag0,MPI_COMM_WORLD, &
                     mpi_4reqs(1),mpi_ierr )
      CALL MPI_IRECV(var(:,0),nz,                                &
                     MPI_REAL8,rproc(1),mpi_tag1,MPI_COMM_WORLD, &
                     mpi_4reqs(2),mpi_ierr )

      CALL MPI_ISEND(var(:,nx),nz,                               &
                     MPI_REAL8,sproc(2),mpi_tag1,mpi_comm_world, &
                     mpi_4reqs(3),mpi_ierr)
      CALL MPI_ISEND(var(:,1),nz,                                &
                     MPI_REAL8,sproc(1),mpi_tag0,mpi_comm_world, &
                     mpi_4reqs(4),mpi_ierr)

      CALL MPI_WAITALL(4,mpi_4reqs,MPI_STATUS_IGNORE,mpi_ierr)

!     CALL MPI_SENDRECV(var(1,:),nz,                      &
!                       MPI_REAL8,sproc(1),0,              &
!                       var(nx+1,:),nz,                  &
!                       MPI_REAL8,rproc(2),0,              &
!                       MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )
!     CALL MPI_SENDRECV(var(nx,:),nz,                  &
!                       MPI_REAL8,sproc(2),1,              &
!                       var(0,:),nz,                      &
!                       MPI_REAL8,rproc(1),1,              &
!                       MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

      if ( rank == 0 ) then
         var(:,0) = zval(:)
      elseif ( rank == nworker ) then
         var(:,nx+1) = zval(:)
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_buffer_nobdy_real_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------------------------------!
!
! MPI COMM routine
!
! - Update right buffer for advection
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_adv_right_buffer_sp (zf,v,rvar)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx
      Implicit none
      real(sp),dimension(0:nx+1,nz) :: zf,v
      real(sp),dimension(3,nz) :: rvar
      real(sp),dimension(2,nz) :: svar2
      real(sp),dimension(1,nz) :: svar1
!----------------------------------------------------------------------!
! Communication: update buffer

      ! svar2(1,:) = ZF - 1
      ! svar2(2,:) =  V - 1
      ! svar1(1,:) = ZF + 2

      svar2(1,:) = zf(nx-1,:)
      svar2(2,:) =  v(nx-1,:)

      CALL MPI_SENDRECV(svar2(:,:),nz*2,                      &
                        MPI_REAL4,sproc(2),0,              &
                        rvar(1:2,:),nz*2,                  &
                        MPI_REAL4,rproc(1),0,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )

      svar1(1,:) = zf(2,:)
      CALL MPI_SENDRECV(svar1,nz,                  &
                        MPI_REAL4,sproc(1),1,              &
                        rvar(3,:),nz,                      &
                        MPI_REAL4,rproc(2),1,              &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr )


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine update_adv_right_buffer_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_x_sp (isign,var_x_sp,p_var_x_sp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(pnx) :: var_x_sp
      real(sp),dimension(nx) :: p_var_x_sp
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var_x_sp(is:ie),                       &
                        subnxall(n),MPI_REAL4,            &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var_x_sp(:) = var_x_sp(is:ie)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var_x_sp(:),                         &
                        nx,MPI_REAL4,                     &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var_x_sp(:),                         &
                        nx,MPI_REAL4,                     &
                        master,1,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var_x_sp(is:ie),                       &
                        subnxall(n),MPI_REAL4,            &
                        n,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var_x_sp(is:ie) = p_var_x_sp(:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_x_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_x_dp (isign,var,p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(dp),dimension(pnx) :: var
      real(dp),dimension(nx) :: p_var
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var(is:ie),                       &
                        subnxall(n),MPI_REAL8,            &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var(:) = var(is:ie)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var(:),                         &
                        nx,MPI_REAL8,                     &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var(:),                         &
                        nx,MPI_REAL8,                     &
                        master,1,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var(is:ie),                       &
                        subnxall(n),MPI_REAL8,            &
                        n,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var(is:ie) = p_var(:)

      endif ! MASTER


      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_x_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_zx_sp (isign,var_zx_sp,p_var_zx_sp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use typeKind
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(nz,pnx) :: var_zx_sp
      real(sp),dimension(nz,nx) :: p_var_zx_sp
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var_zx_sp(:,is:ie),                       &
                        nz*subnxall(n),MPI_REAL4,          &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var_zx_sp(:,:) = var_zx_sp(:,is:ie)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var_zx_sp(:,:),                         &
                        nz*nx,MPI_REAL4,                    &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var_zx_sp(:,:),                         &
                        nz*nx,MPI_REAL4,                    &
                        master,1,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var_zx_sp(:,is:ie),                       &
                        nz*subnxall(n),MPI_REAL4,           &
                        n,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var_zx_sp(:,is:ie) = p_var_zx_sp(:,:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_zx_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_zx_dp (isign,var_zx_dp,p_var_zx_dp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(dp),dimension(nz,pnx) :: var_zx_dp
      real(dp),dimension(nz,nx) :: p_var_zx_dp
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var_zx_dp(:,is:ie),                       &
                        nz*subnxall(n),MPI_REAL8,          &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var_zx_dp(:,:) = var_zx_dp(:,is:ie)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var_zx_dp(:,:),                         &
                        nz*nx,MPI_REAL8,                    &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank /= master ) then ! WORKER
          call MPI_SEND(p_var_zx_dp(:,:),                         &
                        nz*nx,MPI_REAL8,                    &
                        master,1,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var_zx_dp(:,is:ie),                       &
                        nz*subnxall(n),MPI_REAL8,           &
                        n,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var_zx_dp(:,is:ie) = p_var_zx_dp(:,:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_zx_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_xz_sp (isign,var,p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use typeKind
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(pnx,nz) :: var
      real(sp),dimension(nx,nz) :: p_var
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var(is:ie,:),                       &
                        nz*subnxall(n),MPI_REAL4,          &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var(:,:) = var(is:ie,:)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var(:,:),                         &
                        nz*nx,MPI_REAL4,                    &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var(:,:),                         &
                        nz*nx,MPI_REAL4,                    &
                        master,1,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var(is:ie,:),                       &
                        nz*subnxall(n),MPI_REAL4,           &
                        n,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var(is:ie,:) = p_var(:,:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_xz_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_zxb_sp (isign,var_zxb_sp,p_var_zxb_sp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(nz,pnx,max_nbins) :: var_zxb_sp
      real(sp),dimension(nz,nx,max_nbins) :: p_var_zxb_sp
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var_zxb_sp(:,is:ie,:),                       &
                        nz*subnxall(n)*max_nbins,MPI_REAL8,         &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var_zxb_sp(:,:,:) = var_zxb_sp(:,is:ie,:)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var_zxb_sp(:,:,:),                         &
                        nz*nx*max_nbins,MPI_REAL8,                  &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var_zxb_sp(:,:,:),                         &
                        nz*nx*max_nbins,MPI_REAL8,                   &
                        master,0,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var_zxb_sp(:,is:ie,:),                       &
                        nz*subnxall(n)*max_nbins,MPI_REAL8,          &
                        n,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var_zxb_sp(:,is:ie,:) = p_var_zxb_sp(:,:,:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_zxb_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_master_worker_zxb_dp (isign,var_zxb_dp,p_var_zxb_dp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(dp),dimension(nz,pnx,max_nbins) :: var_zxb_dp
      real(dp),dimension(nz,nx,max_nbins) :: p_var_zxb_dp
!----------------------------------------------------------------------!

! isign (+1) : Scattering (MASTER -> WORKER) --------------------------!
      if ( isign == +1 ) then

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_SEND(var_zxb_dp(:,is:ie,:),                       &
                        nz*subnxall(n)*max_nbins,MPI_REAL8,         &
                        n,0,mpi_comm_world,mpi_ierr)
        end do
          is = sloc(rank) ; ie = eloc(rank)
          p_var_zxb_dp(:,:,:) = var_zxb_dp(:,is:ie,:)

      endif ! MASTER

      if ( rank > master ) then ! WORKER
          call MPI_RECV(p_var_zxb_dp(:,:,:),                         &
                        nz*nx*max_nbins,MPI_REAL8,                  &
                        master,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
      endif ! WORKER

! isign (-1) : Gathering (WORKER -> MASTER) ---------------------------!
      elseif ( isign == -1 ) then

      if ( rank > master ) then ! WORKER
          call MPI_SEND(p_var_zxb_dp(:,:,:),                         &
                        nz*nx*max_nbins,MPI_REAL8,                   &
                        master,0,MPI_COMM_WORLD,mpi_ierr)
      endif ! WORKER

      if ( rank == master ) then ! MASTER
        do n = 1,nworker
          is = sloc(n) ; ie = eloc(n)
          call MPI_RECV(var_zxb_dp(:,is:ie,:),                       &
                        nz*subnxall(n)*max_nbins,MPI_REAL8,          &
                        n,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierr)
        end do

          is = sloc(rank) ; ie = eloc(rank)
          var_zxb_dp(:,is:ie,:) = p_var_zxb_dp(:,:,:)

      endif ! MASTER

      endif ! isign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine comm_master_worker_zxb_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_averaging_sp (g_scalar_sp,p_scalar_sp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),intent(in) :: p_scalar_sp
      real(sp) :: g_scalar_sp,scalar_sp

      scalar_sp = p_scalar_sp * nx

      call mpi_reduce(scalar_sp,g_scalar_sp,1,mpi_real4,mpi_sum, &
                      master,mpi_comm_world,mpi_ierr)
      
      if ( rank == master ) then
         g_scalar_sp = g_scalar_sp / float(pnx)
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END subroutine comm_averaging_sp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_averaging_dp (g_scalar_dp,p_scalar_dp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(wp),intent(in) :: p_scalar_dp
      real(wp) :: g_scalar_dp,scalar_dp

      scalar_dp = p_scalar_dp * nx

      call mpi_reduce(scalar_dp,g_scalar_dp,1,mpi_real8,mpi_sum, &
                      master,mpi_comm_world,mpi_ierr)
      
      if ( rank == master ) then
         g_scalar_dp = g_scalar_dp / dfloat(pnx)
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END subroutine comm_averaging_dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comm_gatherv_dp (field, elemsize, sidx, eidx)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: elemsize, sidx, eidx,i
      integer, dimension(:), allocatable :: recvcounts, displs
      real(wp) :: field(0:*)
      real(wp) :: field_local((eidx-sidx+1)*elemsize)
 
      allocate(recvcounts(ncpu), displs(ncpu))

      i = (eidx-sidx+1) * elemsize
      call mpi_allgather(i,1,mpi_integer,recvcounts,1,mpi_integer,  &
                         mpi_comm_world,mpi_ierr)
      i = sidx*elemsize
      call mpi_allgather(i,1,mpi_integer,displs,1,mpi_integer,  &
                         mpi_comm_world,mpi_ierr)
     

      do i = 1,elemsize*(eidx-sidx+1)
         field_local(i) = field(i+elemsize*sidx-1)
      end do

      call mpi_allgatherv (field_local,             &
                           (eidx-sidx+1)*elemsize,  &
                           mpi_real8,               &
                           field,                   &
                           recvcounts,              &
                           displs,                  &
                           mpi_real8,               &
                           mpi_comm_world,          &
                           mpi_ierr )
 
      deallocate(recvcounts,displs)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END subroutine comm_gatherv_dp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine messages (msgs)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit none
      character(*),intent(in) :: msgs 
   
      if ( rank == 0 ) then
         print*,msgs
      end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine messages 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_1d_sp (p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(pnx) :: var
      real(sp),dimension(nx) :: p_var
!----------------------------------------------------------------------!

      call comm_master_worker_x_sp (-1,var,p_var)

      if ( rank == master ) then
!        do j = 1,pnx
!           print*, var(j)
!        end do
         print*,minval(var), maxval(var)
      endif

!     call kid_mpi_barrier

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine check_1d_sp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_1d_dp (p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(dp),dimension(pnx) :: var
      real(dp),dimension(nx) :: p_var
!----------------------------------------------------------------------!

      call comm_master_worker_x_dp (-1,var,p_var)

      if ( rank == master ) then
         do j = 1,pnx
            print*, j,var(j)
         end do
!        print*,minval(var), maxval(var)
      endif

!     call kid_mpi_barrier

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine check_1d_dp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_2d_sp (p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(nz,pnx) :: var
      real(sp),dimension(nz,nx) :: p_var
!----------------------------------------------------------------------!

      call comm_master_worker_zx_sp (-1,var,p_var)

      if ( rank == master ) then
         do k = 1,nz
         do j = 1,pnx
            print*, k,j,var(k,j)
         end do
         end do
!        print*,minval(var), maxval(var)
      endif

      call kid_mpi_barrier

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine check_2d_sp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_2d_dp (p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(dp),dimension(nz,pnx) :: var
      real(dp),dimension(nz,nx) :: p_var
!----------------------------------------------------------------------!

      call comm_master_worker_zx_dp (-1,var,p_var)

      if ( rank == master ) then
         do k = 1,nz
         do j = 1,pnx
            print*, k,j,var(k,j)
         end do
         end do
!        print*,minval(var), maxval(var)
      endif

      call kid_mpi_barrier

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine check_2d_dp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_2d_reverse_sp (p_var)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Use parameters, only: nz,nx,pnx
      Implicit none
      integer :: i,j,k,n,isign,is,ie
      real(sp),dimension(pnx,nz) :: var
      real(sp),dimension(nx,nz) :: p_var
!----------------------------------------------------------------------!

      call comm_master_worker_xz_sp (-1,var,p_var)

      if ( rank == master ) then
         do k = 1,nz
         do j = 1,pnx
            print*, k,j,var(j,k)
         end do
         end do
!        print*,minval(var), maxval(var)
      endif

!     call kid_mpi_barrier

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine check_2d_reverse_sp 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine KiD_mpi_barrier 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
      call mpi_barrier(mpi_comm_world,mpi_ierr)
!     stop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine KiD_mpi_barrier 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine KiD_mpi_finalize 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
      call mpi_barrier(mpi_comm_world,mpi_ierr)
      call mpi_finalize(mpi_ierr) 
  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine KiD_mpi_finalize 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#endif

end module mpi_mod
