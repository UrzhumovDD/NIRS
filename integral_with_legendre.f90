module legendre
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains

! calculating polynomial 

real(dp) function Lezh(n,x) result(c)
real(dp), intent(in   )     :: x
integer,  intent(in   )     :: n
real(dp)                    :: a, b
integer                     :: i

if ( n == 0 ) then
    
c = 1.0_dp

else if ( n == 1 ) then
    
c = x

else
    
a = 1.0_dp
b = x

do i = 2, n
    
c = b * x * ( 2.0_dp*i - 1.0_dp ) / i - a * ( i - 1.0_dp ) / i
a = b
b = c

end do

end if

end function Lezh


!calculating nods and their weights


subroutine nodes_weights( n, a, b, nodes, weights )
integer,  intent(in   )              :: n
real(dp), intent(in   )              :: a,b
real(dp), intent(  out), allocatable :: nodes(:), weights(:)
real(dp), allocatable                :: nodes1(:)
real(dp)                             :: eps = 1e-12
integer                              :: i 
real(dp), parameter                  :: PI = ACOS( -1.D0 )

!first approximation

nodes = [ ( -cos( PI * ( 4.0_dp * i - 1.0_dp ) / ( 4.0_dp * n + 2.0_dp ) ), i = 1, n ) ]
allocate( nodes1( n ) )
nodes1 = 10

! cycle with Newton method for nods

do while ( abs( maxval( nodes - nodes1 ) ) >= eps)
nodes1 = nodes
do i = 1 ,n
nodes( i ) = nodes1( i ) - lezh( n,nodes1( i ) ) * ( 1.0_dp - nodes1( i ) ** 2.0_dp ) / &
    ( ( lezh( n - 1, nodes1( i ) ) - nodes1( i ) * lezh( n, nodes1( i ) ) )* n )
end do 
end do

!weights

allocate( weights( n ),source = [ ( 2_dp / ( 1_dp-nodes( i ) ** 2_dp ) / &
    ( ( ( lezh( n - 1, nodes( i ) ) - nodes( i ) * lezh( n, nodes( i ) ) ) * n / &
    ( 1_dp - nodes( i ) ** 2 ) ) **2 ), i = 1, n ) ] )

!affine transformations

nodes = ( nodes * ( b - a ) + a + b) * 0.5_dp
end subroutine nodes_weights


end module legendre


module integral
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
contains


subroutine F( nodes, Func )
real(dp),intent(in   ),allocatable :: nodes(:)
real(dp),intent(  out),allocatable :: Func(:)

Func = tanh( nodes )

end subroutine F


real(dp) function analytically( a, b )
real(dp),intent(in   ) :: a, b

analytically = Log( cosh( b ) )- Log( cosh( a ) )

end function analytically


real(dp) function numerically( Func, weights, a, b )
real(dp),intent(in   ),allocatable  :: Func(:) , weights(:)
real(dp),intent(in   )              :: a, b

numerically = sum( weights * Func ) * ( b - a ) * 0.5_dp

end function numerically 


subroutine poly( nodes, Func )
real(dp),intent(in   ),allocatable :: nodes(:)
real(dp),intent(  out),allocatable :: Func(:)

Func = nodes ** 3 / 4_dp + 3_dp * nodes ** 2 / 4_dp - 3_dp * nodes / 2_dp - 2_dp

end subroutine poly


real(dp) function an_poly( a, b )
real(dp),intent(in   ) :: a, b

an_poly = ( b ** 4 - a ** 4 ) / 16_dp + ( b ** 3 - a ** 3 ) / 4_dp - 3_dp * ( b ** 2 - a ** 2 ) / 4_dp &
    - 2_dp * ( b - a)  

end function an_poly


real(dp) function num_poly( Func, weights, a, b )
real(dp),intent(in   ),allocatable  :: Func(:) , weights(:)
real(dp),intent(in   )              :: a, b

num_poly = sum( weights * Func ) * ( b - a ) * 0.5_dp

end function num_poly


subroutine Fhard( nodes, Func )
real(dp),intent(in   ),allocatable :: nodes(:)
real(dp),intent(  out),allocatable :: Func(:)

Func = nodes * atan( nodes ) + nodes * cos( nodes )

end subroutine Fhard


real(dp) function an_hard( a, b )
real(dp),intent(in   ) :: a, b

an_hard = b * sin( b ) - a * sin( a ) + cos( b ) - cos( a ) &
    + ( ( b ** 2 + 1_dp) * atan( b ) - b ) / 2_dp &
    - ( ( a ** 2 + 1_dp) * atan( a ) - a ) / 2_dp

end function an_hard


real(dp) function num_hard( Func, weights, a, b )
real(dp),intent(in   ),allocatable  :: Func(:) , weights(:)
real(dp),intent(in   )              :: a, b

num_hard = sum( weights * Func ) * ( b - a ) * 0.5_dp

end function num_hard

end module integral


program test
use legendre, only : nodes_weights
use integral, only : F, analytically, numerically, poly, num_poly, an_poly , num_hard, an_hard, Fhard
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nodes(:),weights(:),Func(:)
real(dp)              :: a, b
integer               :: n

write(*,*) 'order'
read(*,*) n
write(*,*) 'left border'
read(*,*) a
write(*,*) 'right border'
read(*,*) b
call nodes_weights (n, a, b, nodes, weights )

call Fhard( nodes, Func )
!print *, nodes
!print *, weights
print *, 'x*arctan(x) + x*cos(x)'
print *, 'analytically'
print *, an_hard( a, b )
print *, 'numerically'
print *, num_hard( Func, weights, a, b )
print *, 'absolute difference'
print *, an_hard( a, b ) - num_hard (Func, weights, a, b )
print *, 'relative difference'
print *, ( an_hard( a, b ) - num_hard( Func, weights, a, b ) ) / an_hard( a, b )
deallocate(Func)

print *, 'x^3 / 4 + 3 * x^2 / 4 - 3 * x / 2 - 2'
call poly( nodes, Func )
print *, 'analytically'
print *, an_poly( a, b )
print *, 'numerically'
print *, num_poly( Func, weights, a, b )
print *, 'absolute difference'
print *, an_poly( a, b ) - num_poly (Func, weights, a, b )
print *, 'relative difference'
print *, ( an_poly( a, b ) - num_poly( Func, weights, a, b ) ) / an_poly( a, b )
deallocate(Func)

print *, 'tanh(x)'
call F( nodes, Func )
print *, 'analytically'
print *, analytically( a, b )
print *, 'numerically'
print *, numerically( Func, weights, a, b )
print *, 'absolute difference'
print *, analytically( a, b ) - numerically (Func, weights, a, b )
print *, 'relative difference'
print *, ( analytically( a, b ) - numerically( Func, weights, a, b ) ) / analytically( a, b )
deallocate(Func)

end program test
