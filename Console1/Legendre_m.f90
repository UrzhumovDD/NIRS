module legendre

    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    implicit none

    private
    public  :: nodes_values, Pn, All_poly
    
    contains

    !> Calculate Legendre polynomial of the n-th order at the point x
    pure elemental real(dp) function Pn(n, x) 
        
        real(dp), intent(in)    :: x        !! argument
        integer,  intent(in)    :: n        !! polynomial order
        real(dp)                :: a, b     !  variables for recurrence relation: a = P_{n-2}, b = P_{n-1}
        integer                 :: i        !  loop counter

        ! Bonnet's recursion formula
        if (n < 0) then
            error stop "Negative order of Legendre polynomial"
        else if (n == 0) then
            Pn = 1.0_dp
        else if (n == 1) then
            Pn = x
        else
            a = 1.0_dp
            b = x
            do i = 2, n
                Pn = b * x * real(2 * i - 1, dp) / real(i, dp) - a * real(i - 1, dp) / real(i, dp)
                a = b
                b = Pn
            end do
        end if

    end function Pn
    
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
        
    !calculating polynomials of all orders from 0 to SIZE(Poly, 1) - 1 in nodes
    subroutine All_poly(Poly)
    
        real(dp), intent(inout)     :: Poly(:,:)    ! values of legendre polynomials in nodes
        real(dp), allocatable       :: nodes(:)     ! values of nodes
        integer                     :: i,j          ! loop counters
        
        call nodes_values( SIZE(Poly, 2), nodes)
        do i = 1,SIZE(Poly, 1) ! orders of polynomials
            do j = 1,SIZE(Poly, 2) 
                Poly(i,j) = Pn(i - 1, nodes(j) )
            end do
        end do
        
    end subroutine All_poly
    
end module legendre