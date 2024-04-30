Module activated_bin_dist

  Use typeKind
  Use special_kid, only: gamma_kid
  Use class_species, only: species

  Implicit None

  real, parameter :: pi=3.14159


contains

  subroutine set_activated_bins_dist(mass_bin, nd_bin, q, n, mu)
    ! Code to set the size distribution of the activated mass
    ! and number in the TAU 2M  bin microphysics. Code is based
    ! set_init_bin routine
    !
    ! NB: Currently only set up to work for TAU warm bin
    !

    integer, parameter :: nbins=34
    real(wp), intent(inout) :: mass_bin(nbins), nd_bin(nbins)
    real(wp), intent(in) ::             &
       q                            & ! bulk mass mixing ratio
       ,n                            ! bulk number concentration
    real(wp), intent(in) :: mu                            ! shape parameter
    
    real(wp) :: lambda_0, n_0
    integer :: k
    real(wp) :: dk, dkm1

    ! TAU values
    real(wp) :: D0=3.125e-6       ! Diameter of smallest bin upper bound
    real(wp) :: m0
    real(wp) :: mk(nbins)
    real(wp) :: dbar, mbar
    real(wp), parameter :: rhow=1000.
    real(wp) :: c=pi*rhow/6.     ! M=cD^d
    real(wp) :: d=3               ! M=cD^d
    real(wp) :: sumq, sumn

    lambda_0 = (n*c*gamma_kid(1.+mu+d)/(q*gamma_kid(1.+mu)))**(1./d)
    n_0 = n*lambda_0**(1.+mu)/gamma_kid(1.+mu)

    print *, n, q
    
    m0=c*D0**d
    mk(1)=m0
    do k=2,nbins
      mk(k)=mk(k-1)*2.
    end do

    dk=(mk(1)/c)**(1./d)
    mbar=.5*(c*dk**d)
    nd_bin(1)=n_0*lambda_0**(-(1.+mu))* &
         (gamma_kid(mu+1.)- gamma_kid(mu+1., lambda_0*dk))
    mass_bin(1)=n_0*c*lambda_0**(-(1.+mu+d))* &
         (gamma_kid(d+mu+1.) - gamma_kid(d+mu+1., lambda_0*dk))
    !mass_bin(1)=nd_bin(1)*mbar
    do k=2,nbins
      dk=(mk(k)/c)**(1./d)
      dkm1=(mk(k-1)/c)**(1./d)
      mbar=.5*(c*dk**d+c*dkm1**d)
      nd_bin(k)=n_0*lambda_0**(-(1.+mu))* &
         (gamma_kid(mu+1., lambda_0*dkm1) - gamma_kid(mu+1., lambda_0*dk))
      mass_bin(k)=n_0*c*lambda_0**(-(1.+mu+d))* &
         (gamma_kid(d+mu+1., lambda_0*dkm1) - gamma_kid(d+mu+1., lambda_0*dk)) 
      !mass_bin(k)=nd_bin(k)*mbar
      !nd_bin(k)=mass_bin(k)/mbar
      if (nd_bin(k)<2.)then
        nd_bin(k)=0.
        mass_bin(k)=0.
      end if
    end do

    sumq=0.
    sumn=0.
    do k=1,nbins
      sumq=sumq+mass_bin(k)
      sumn=sumn+nd_bin(k)
    !  print *, 'mass_bin, k', mass_bin(k), k
    end do
    print *, 'mass', sumq, q
    print *, 'number', sumn, n
    
  end subroutine set_activated_bins_dist

end Module activated_bin_dist




  
