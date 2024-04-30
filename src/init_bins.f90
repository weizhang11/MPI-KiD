Module init_bins

  Use typeKind
  Use special_kid, only: gamma_kid
  Use class_species, only: species

  Implicit None

  real, parameter :: pi=3.14159


contains

  subroutine set_init_bins(hydrometeor, q, n, mu)
    ! Code to initialize a distribution for bin 
    ! microphysics.  
    !
    ! NB: Currently only set up to work for TAU warm bin
    !

    type(species), intent(inout) :: &
       hydrometeor                    ! hydrometeor to initialise
    real(wp), intent(in) ::             &
       q                            & ! bulk mass mixing ratio
       ,n                           & ! bulk number concentration
       ,mu                            ! shape parameter

    real(wp) :: lambda_0, n_0
    integer :: k
    real(wp) :: dk, dkm1

    ! TAU values
    real(wp) :: D0 = 3.125e-6       ! Diameter of smallest bin lower bound
    real(wp) :: m0
    integer, parameter :: nbins=34
    real(wp) :: mk(nbins)
    real(wp) :: dbar, mbar
    real(wp), parameter :: rhow=1000.
    real(wp) :: c=pi*rhow/6.     ! M=cD^d
    real(wp) :: d=3               ! M=cD^d

    real(wp) :: qk(nbins), nk(nbins)
    real(wp) :: sumq, sumn

    lambda_0 = (n*c*gamma_kid(1.+mu+d)/(q*gamma_kid(1.+mu)))**(1./d)
    n_0 = n*lambda_0**(1.+mu)/gamma_kid(1.+mu)

    m0=(c*(D0)**d)*2.0
    mk(1)=m0
    do k=2,nbins
      mk(k)=mk(k-1)*2.
    end do

    dk=(mk(1)/c)**(1./d)
    mbar=.5*(c*dk**d)
    nk(1)=n_0*lambda_0**(-(1.+mu))* &
         (gamma_kid(mu+1.)- gamma_kid(mu+1., lambda_0*dk))
    qk(1)=n_0*c*lambda_0**(-(1.+mu+d))* &
         (gamma_kid(d+mu+1.) - gamma_kid(d+mu+1., lambda_0*dk))
    !qk(1)=nk(1)*(1.000001*(c*D0**d))
    do k=2,nbins
      dk=(mk(k)/c)**(1./d)
      dkm1=(mk(k-1)/c)**(1./d)
      
      mbar=.5*(c*dk**d+c*dkm1**d)
      nk(k)=n_0*lambda_0**(-(1.+mu))* &
         (gamma_kid(mu+1., lambda_0*dkm1) - gamma_kid(mu+1., lambda_0*dk))
      qk(k)=n_0*c*lambda_0**(-(1.+mu+d))* &
         (gamma_kid(d+mu+1., lambda_0*dkm1) - gamma_kid(d+mu+1., lambda_0*dk)) 
      !qk(k)=nk(k)*mbar
      !nk(k)=qk(k)/mbar
      if (nk(k)<2.)then
        nk(k)=0.
        qk(k)=0.
      end if
    end do

    sumq=0.
    sumn=0.
    do k=1,nbins
      hydrometeor%moments(k,1)=qk(k)
      hydrometeor%moments(k,2)=nk(k)
      print *, hydrometeor%moments(k,1), hydrometeor%moments(k,2)
      sumq=sumq+hydrometeor%moments(k,1)
      sumn=sumn+hydrometeor%moments(k,2)
    end do
    print *, sumq, sumn

  end subroutine set_init_bins

end Module init_bins



  
