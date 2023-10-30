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
integer                            :: i

Func = tanh(nodes)

end subroutine F

real(dp) function analytically( a, b )
real(dp),intent(in   ) :: a, b

analytically = Log( cosh( b ) )- Log( cosh( a ) )

end function analytically

real(dp) function numerically( Func, weights, a, b )
real(dp),intent(in   ),allocatable  :: Func(:) , weights(:)
real(dp),intent(in   )              :: a, b
integer                             :: i

numerically = sum( weights * func ) * ( b - a ) * 0.5_dp

end function numerically 
end module integral    
    
program test
use legendre, only : nodes_weights
use integral, only : F, analytically, numerically
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
call F( nodes, Func )
print *, nodes
print *, weights
print *, Func
print *, analytically( a, b )
print *, numerically( Func, weights, a, b )
print *, 'absolute difference'
print *, analytically( a, b ) - numerically (Func, weights, a, b )
print *, 'relative difference'
print *, ( analytically( a, b ) - numerically( Func, weights, a, b ) ) / analytically( a, b )

end program test
