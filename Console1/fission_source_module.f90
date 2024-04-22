module Fission_source
    
    use Integration_with_weights, only: Integral
    use Calc_of_nodes,            only: C_o_n
    use Calc_of_weights,          only: C_o_W
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    
    implicit none
    private
    public :: Fsource
    
    contains
    
    subroutine Fsource(Flux, Nuf, cross_sec_fis, Hi, Keff, res )
    real(dp),              intent(in   ) :: Flux(:,:,:), Nuf(:,:), cross_sec_fis(:,:), Hi(:)    ! Hi is a fission cpectrum, Nuf is a the average number of neutrons per act of fission 
    real(dp), allocatable, intent(inout) :: res(:,:)                                            ! final value
    real(dp), allocatable                :: nodes(:), weights(:), Integr(:,:)
    real(dp)                             :: Keff
    integer                              :: i, j                                                ! loop counters
    
    call C_o_n( SIZE(Flux, 2), nodes)
    call C_o_w( nodes, weights)
    allocate( Integr( SIZE(Flux, 1), SIZE(Flux, 3) ) )
    allocate( res( SIZE(Flux, 1), SIZE(Flux, 3) ) )
    ! integration and dividing by 2
    do i = 1, SIZE(Flux, 1)
        do j = 1, SIZE(Flux, 3)
            Integr(i,j) = Integral( Flux( i,:, j ) , weights) / 2_dp
        end do
    end do
    do i = 1, SIZE(Flux, 1)
        do j = 1, SIZE(Flux, 3)
            res(i,j) = Hi(j) / Keff * SUM( Integr(i, :) * Nuf(i, :) * cross_sec_fis(i, :) )
        end do
    end do
    
    end subroutine Fsource
    
end module Fission_source
    
program FS

    use Fission_source, only : Fsource
    use, intrinsic        :: iso_fortran_env, only: dp=>real64
    
    implicit none 
    
    real(dp), allocatable :: Flux(:,:,:), Nuf(:,:), cross_sec_fis(:,:), Hi(:), Fiss_S(:,:)
    real(dp)              :: Keff
    integer               :: ncord, nenergy, nangles
    
    ncord = 3
    nenergy = 5
    nangles = 5
    Keff = 1.0_dp
    allocate( Flux(ncord, nangles, nenergy ), source = 1.0_dp )
    allocate( Nuf(ncord, nenergy ), source = 1.0_dp )
    allocate( cross_sec_fis(ncord, nenergy ), source = 1.0_dp )
    allocate( Hi( nenergy ), source = 1.0_dp )
    call Fsource(Flux, Nuf, cross_sec_fis, Hi, Keff, Fiss_S )
    deallocate(Flux)
    deallocate(Nuf)
    deallocate(cross_sec_fis)
    deallocate(Hi)
    
end program FS