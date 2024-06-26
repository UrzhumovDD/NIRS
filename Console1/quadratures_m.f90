module Quadratures

    use legendre, only : Pn
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none
    
    private
    public  :: weights_values
    
    contains
    
    subroutine weights_values(nodes,weights)
    
        real(dp), intent(  out), allocatable    :: weights(:) !values of weights
        real(dp), intent(in   )                 :: nodes(:)   !values of nodes
        integer                                 :: i, n       !i is a loop counter, n is number of nodes
        
        n = SIZE(nodes)
        !weights of the Gauss-Legendre quadrature
        !the derivative of the polynomial calculated by the recurrent formula
        allocate( weights(n), source = [ ( 2_dp / ( 1_dp - nodes(i) ** 2_dp ) / &
            ( ( ( Pn(n - 1, nodes(i)) - nodes(i) * Pn(n, nodes(i)) ) * real(n, dp) / &
            ( 1_dp - nodes(i) ** 2 ) ) ** 2 ), i = 1, n) ] )
            
    end subroutine weights_values
        
end module Quadratures