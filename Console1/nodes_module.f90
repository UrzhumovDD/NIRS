module Calc_of_nodes

    use legendre, only : Pn
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none
    
    private
    public  :: C_o_n
    
    contains
    
    subroutine C_o_n( n, nodes)
    
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
        do while ( ABS(MAXVAL( nodes - nodes1 ) ) >= eps)
            nodes1 = nodes
            do i = 1 ,n        
                nodes(i) = nodes1(i) - Pn(n, nodes1(i) ) * ( 1.0_dp - nodes1(i) ** 2.0_dp ) / &
                    (( Pn(n - 1, nodes1(i) ) - nodes1(i) * Pn(n, nodes1(i) ) ) * real(n, dp))       
            end do
        end do
        
    end subroutine C_o_n
    
end module Calc_of_nodes