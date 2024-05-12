module Scaterring_source
    
    use Quadratures,                   only : C_o_n, C_o_W
    use Legendre,                        only : All_poly
    use Integration_with_weights,        only : Integral
    use, intrinsic :: iso_fortran_env,   only : dp=>real64
    
    implicit none
    
    private
    public  :: Ssource
    
    contains
    
    subroutine Ssource(Flux, Cross_sec_scat, order, energy, res)
    
    real(dp),              intent(in   ) :: Flux(:,:,:), Cross_sec_scat(:,:,:,:)
    real(dp), allocatable                :: Poly(:,:), nodes(:), weights(:), Integr(:,:,:)   !Integr is auxiliary integral of flux * Pn
    real(dp), allocatable, intent(inout) :: res(:,:,:)                                       !res(result) is final value of Ssource
    integer,               intent(in   ) :: order, energy                                    !max order of Legendre's polynomial and energy group to which neutrons transfer
    integer                              :: i, j, k                                          !loop counters
    
    !Flux(ncord, nangles, nenergy)
    !Cross_sec_scat(ncord,nenergy (out),nenergy (in),order)
    
    allocate( res( SIZE(Flux, 1), SIZE(Flux, 2), SIZE(Flux, 3) ), source = 0.0_dp )
    allocate( Poly( order + 1 , SIZE(Flux, 2) ) )
    call C_o_n( SIZE(Flux, 2), nodes)
    call C_o_w( nodes, weights)
    call All_poly( Poly)
    allocate( Integr( SIZE(Flux, 1), order + 1, SIZE(Flux, 3)) )
    do i = 1, SIZE(Flux, 1)
        do j = 1, order + 1
            do k = 1, SIZE(Flux, 3)
                Integr(i,j,k) = Integral( Poly(j,:) * Flux(i, :, k ), weights)
                res(i,:,k) = res(i,:,k) +  Poly(j,:) * Cross_sec_scat(i, k, energy, j) * real( 2 * ( j - 1 ) + 1, dp) / 2_dp * Integr(i,j,k)
            end do
        end do
    end do
    
    end subroutine Ssource
    
end module Scaterring_source

program test

    use Scaterring_source, only : Ssource
    use Quadratures, only     : C_o_n
    use Legendre, only     : All_poly
    use, intrinsic        :: iso_fortran_env, only: dp=>real64
    
    implicit none
    
    real(dp), allocatable :: nodes(:), Flux(:,:,:), Cross_sec_scat(:,:,:,:), SS(:,:,:), Poly(:,:), SS_an(:,:,:)
    real(dp)              :: integr_analytic(6)
    integer               :: order,ncord,nenergy,nangles, i, j, k,energy
    
    !values for variables used to check module
    order = 5
    ncord = 3
    nenergy = 5
    energy = 1
    nangles = 5
    call C_o_n(nangles,nodes)
    allocate(Flux(ncord, nangles, nenergy))
    do i = 1,SIZE(Flux,1)
        do k = 1,SIZE(Flux,3)
            Flux(i,:,k) = COS(nodes)
        end do
    end do
    allocate(Cross_sec_scat(ncord,nenergy,nenergy,order + 1),source = 1.0_dp)
    !analyticaly calculated integral of Pn*cos(mu)
    integr_analytic(1) = 2_dp * SIN(1.0_dp)
    integr_analytic(2) = 0_dp
    integr_analytic(3) = 6_dp * COS(1.0_dp) - 4_dp * SIN(1.0_dp)
    integr_analytic(4) = 0_dp
    integr_analytic(5) = 122_dp * SIN(1.0_dp) - 190_dp * COS(1.0_dp)
    integr_analytic(6) = 0_dp
    allocate( Poly(order + 1, nangles) )
    call All_poly( Poly)
    call Ssource( Flux, Cross_sec_scat, order, energy, SS)
    allocate( SS_an(ncord, nangles, nenergy ), source = 0.0_dp)
    do i = 1, SIZE(Flux, 1)
        do j = 1, order + 1
           do k = 1, SIZE(Flux, 3)
               SS_an(i,:,k) = SS_an(i,:,k) + Poly(j,:) * Cross_sec_scat(i, k, energy, j) * ( 2_dp * ( j - 1_dp ) + 1_dp) / 2_dp * integr_analytic(j)
           end do
      end do
    end do
    !absolute difference
    write(*,*) MAXVAL(ABS( SS_an - SS ) )
    !relative difference
    write(*,*) MAXVAL(ABS( (SS_an - SS) / SS ) )
    deallocate(Flux)
    deallocate(Cross_sec_scat)
    deallocate(SS_an)
    deallocate(nodes)
    deallocate(Poly)
    deallocate(SS)
    
end program test
