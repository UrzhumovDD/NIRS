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


subroutine nodes_weights( n, nodes, weights)
integer,  intent(in   )              :: n
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


end subroutine nodes_weights

subroutine Poly_in_cord(Func,cord)
real(dp),allocatable, intent(inout) :: Func(:,:)
real(dp),allocatable, intent(in   ) :: Cord(:)
integer              :: i,j

do i = 1,SIZE(Func,1)

do j = 1,SIZE(Func,2)

Func(i,j) = Lezh(i - 1,Cord(j))

end do

end do
end subroutine Poly_in_cord


end module legendre



program test
use legendre, only : nodes_weights,Poly_in_cord
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nodes(:),weights(:),Func(:,:),Cord(:)
integer               :: m,ncord
integer              :: i,j

write(*,*) 'order'
read(*,*) m
write(*,*) 'number of coordinates'
read(*,*) ncord
allocate( Cord( ncord ) )
read(*,*) Cord
allocate( Func( m + 1 , ncord ) )
call nodes_weights (m, nodes, weights)
call Poly_in_cord(Func,cord)
do i = 1,SIZE(Func,1)

do j = 1,SIZE(Func,2)

write(*,*) Func(i,j)

end do

end do

deallocate(Cord)
deallocate(Func)
end program test
