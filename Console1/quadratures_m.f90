module Quadratures

    use legendre, only : Pn
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none
    
    private
    public  :: nodes_values, weights_values
    
    contains
    
    subroutine nodes_values( n, nodes)
    
        integer,  intent(in   )              :: n                   !number of nodes
        real(dp), intent(  out), allocatable :: nodes(:)            !values of nodes
        real(dp), allocatable                :: nodes1(:)           !auxiliary values for the loop
        real(dp)                             :: eps = 1e-12         !accuracy
        integer                              :: i                   !loop counter
        real(dp), parameter                  :: PI = ACOS( -1.D0 )  !the number of pi
        
        !first approximation
        nodes = [ ( -COS( PI * ( 4.0_dp * real(i, dp) - 1.0_dp ) / ( 4.0_dp * real(n, dp) + 2.0_dp ) ), i = 1, n ) ]
        !first value for auxiliary massive is 10 which is greater then any value of first approximation modulo
        allocate( nodes1( n ), source = 10.0_dp )
        !loop with iterative Newton method for nodes
        do while ( MAXVAL(ABS( nodes - nodes1 ) ) >= eps)
            nodes1 = nodes
            do i = 1 ,n        
                nodes(i) = nodes1(i) - Pn(n, nodes1(i) ) * ( 1.0_dp - nodes1(i) ** 2.0_dp ) / &
                    (( Pn(n - 1, nodes1(i) ) - nodes1(i) * Pn(n, nodes1(i) ) ) * real(n, dp))       
            end do
        end do
        
    end subroutine nodes_values
    
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