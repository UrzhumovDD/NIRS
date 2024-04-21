module legendre

    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    implicit none

    private
    public  :: Pn
    
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

end module legendre