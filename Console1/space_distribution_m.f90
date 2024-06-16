module space
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    implicit none

    private
    public  :: Space_distrib
    
    contains
    
    subroutine Space_distrib(Flux, SS, FS, Cross_sec_tot, nodes,h)
    real(dp),              intent(in   ) :: Cross_sec_tot(:,:,:), SS(:,:,:), FS(:,:), nodes(:), h ! h - step of space grid
    real(dp), allocatable, intent(inout) :: Flux(:,:,:)                                           ! final result
    integer                              :: i, j, k                                               ! loop counters
    real(dp), parameter                  :: PI = ACOS( -1.D0 )                                    !the number of pi
    
    !Flux(ncord, nangles, nenergy)
    !only vacuum in boundary conditions
    
    do j = 1,SIZE(Flux,2)
        do k = 1,SIZE(Flux,3)
            if(j > SIZE(Flux,2)/2) then ! positive angles
                Flux(1,j,k) = 0.0_dp
                do i = 2,SIZE(Flux,1)
                    Flux(i,j,k) = ( SS(i,j,k) + SS(i - 1,j,k) + ( FS(i,k) + FS(i - 1,k) ) / (2.0_dp * PI * SUM(FS)) ) &
                        * h / (2.0_dp * nodes(j) + h * ( Cross_sec_tot(i,j,k) + Cross_sec_tot(i - 1,j,k) ) / 2.0_dp) &
                        + Flux(i - 1,j,k) * ( 4.0_dp * nodes(j) - h * ( Cross_sec_tot(i,j,k) + Cross_sec_tot(i - 1,j,k) ) ) &
                        / ( 4.0_dp * nodes(j) + h * (Cross_sec_tot(i,j,k) + Cross_sec_tot(i - 1,j,k) ) )
                end do
            else ! negative angles
            Flux(SIZE(Flux,1),j,k) = 0.0_dp
                do i = SIZE(Flux,1) - 1, 1, -1
                    Flux(i,j,k) = ( SS(i,j,k) + SS(i + 1,j,k) + ( FS(i,k) + FS(i + 1,k) ) / (2.0_dp * PI * SUM(FS)) ) &
                        * h / ( h * ( Cross_sec_tot(i,j,k) + Cross_sec_tot(i + 1,j,k) ) / 2.0_dp - 2.0_dp * nodes(j)) &
                        + Flux(i + 1,j,k) * ( 4.0_dp * nodes(j) + h * ( Cross_sec_tot(i,j,k) + Cross_sec_tot(i + 1,j,k) ) ) &
                        / ( 4.0_dp * nodes(j) - h * ( Cross_sec_tot(i,j,k) + Cross_sec_tot(i + 1,j,k) ) )
                end do
            end if
        end do
    end do
    
    end subroutine Space_distrib
    
end module space

program test

    use Space, only             : Space_distrib
    use Scaterring_source, only : Ssource
    use Fission_source, only    : Fsource
    use Legendre, only       : nodes_values
    
    use, intrinsic      :: iso_fortran_env, only: dp=>real64
    
    implicit none
    
    real(dp), allocatable :: nodes(:), Flux(:,:,:), Flux_s(:,:,:), Cross_sec_scat(:,:,:,:), SS(:,:,:), Poly(:,:), FS(:,:), Nuf(:,:), Cross_sec_fis(:,:), Hi(:),Cross_sec_tot(:,:,:)
    integer               :: ncord, nenergy, nangles
    real(dp)              :: h = 1e-2_dp, gamma = 6.02_dp / 235 * 19.05_dp * 10 ! h - step of space grid, gamma - nuclear density in 10**-24 cm**-3
    real(dp), parameter   :: eps = 1e-2_dp                                     ! difference between iterations
    integer, parameter    :: order = 5                                          ! order of Legendre polynomial
    
    ncord = 10000
    nenergy = 2
    nangles = 8
    allocate(Cross_sec_tot(ncord, nangles, nenergy), Source = 3.09_dp * gamma )
    Cross_sec_tot(:,:,2) = 1.89_dp * gamma
    allocate(Cross_sec_scat(ncord,nenergy,nenergy,order + 1),Source = 0.0_dp)
    Cross_sec_scat(:,1,2,:) = 1.8_dp * gamma
    allocate(Nuf(ncord, nenergy ), source = 2.8_dp )
    Nuf(:,2) = 2.52_dp
    allocate(Cross_sec_fis(ncord, nenergy ), Source = 1.24_dp * gamma )
    Cross_sec_fis(:, 2 ) = 1.54_dp * gamma
    allocate(Hi(nenergy ), source = 0.53_dp )
    Hi(2) = 0.47_dp
    allocate(Flux(ncord, nangles, nenergy), Source = 1.0_dp)
    allocate(Flux_s(ncord, nangles, nenergy), Source = 0.0_dp)
    call nodes_values(nangles, nodes)
    
    do while ( MAXVAL( ABS(Flux-Flux_s) ) >= eps )
        Flux_s = Flux
        call Ssource( Flux, Cross_sec_scat, order, SS)
        call Fsource( Flux, Nuf, Cross_sec_fis, Hi, FS )
        call Space_distrib(Flux, SS, FS, Cross_sec_tot, nodes,h)
        deallocate (SS)
        deallocate (FS)
    end do
    
    call Fsource( Flux, Nuf, Cross_sec_fis, Hi, FS )
    write(*,*) 'K'
    write(*,*) SUM(FS)
    write(*,*) 'thermal'
    write(*,*) SUM(Flux(:,:,1),2)
    write(*,*) 'fast'
    write(*,*) SUM(Flux(:,:,2),2)
    
end program test