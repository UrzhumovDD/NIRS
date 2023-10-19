module legendre
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none
contains
! calculating polynomial        
real(dp) function Lezh(n,x) result(c)
real(dp) , intent(in)    :: x
integer , intent(in)     :: n
real(dp)                 :: a, b
integer                  :: i
if( n == 0 ) then
c = 1.0
else if( n == 1 ) then
c = x
else
a = 1.0
b = x
do i = 2, n
c = b*x*( 2*i - 1 )/i - a*( i - 1 )/i
a = b
b = c
end do
end if
return
end function Lezh
!calculating nods and their weights
subroutine nodes_weights(n,a,b,nodes,weights) 
integer, intent(in), optional     :: n
real(dp), intent(in), optional        :: a,b
real(dp)                          :: eps = 1e-8_dp
integer                           :: i 
real, parameter                   :: PI = DACOS(-1.D0)
real(dp),allocatable,intent(out)  :: nodes(:) , weights(:)
real(dp),allocatable              :: nodes1(:)
!first approximation
allocate( nodes(n) )
nodes = [( -cos( PI*(4*i - 1)/(4*n + 2) ), i = 1 , n )]
allocate( nodes1(n) )
nodes1 = 10
! cycle with Newton method for nods
do while (abs(maxval(nodes - nodes1)) >= eps)   
nodes1 = nodes
do i = 1 , n
nodes(i) = nodes1(i) - lezh(n,nodes1(i)) * (1 - nodes1(i)**2) / &
    ( ( lezh(n - 1,nodes1(i)) - nodes1(i) * lezh(n,nodes1(i)) )* n )
end do 
end do
!weights
allocate( weights(n),source = [(2/(1-nodes(i)**2)/ &
    ( ( (lezh(n - 1,nodes(i))-nodes(i) * lezh(n,nodes(i)) )* n / (1-nodes(i)**2))**2),i = 1, n)] )
!affine transformations 
nodes = (nodes * (b-a) + a + b)/2
end subroutine nodes_weights  
end module legendre
  
module integral  
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
contains 
    
subroutine F(nodes,Func)
real(dp),intent(in),allocatable    :: nodes(:)
integer                :: i
real(dp),allocatable,intent(out) :: Func(:)
allocate(Func(size(nodes)))
do i =1,size(nodes)
Func(i) = TANH(nodes(i))    
end do
end subroutine F

real function analytically(a,b)
real(dp),intent(in) :: a, b
analytically = Log(cosh(b))- Log(cosh(a))
end function analytically

real function numerically(Func,weights,a,b)
real(dp),allocatable,intent(in)  :: Func(:) , weights(:)
real(dp),intent(in)  :: a, b
integer :: i
numerically = 0
do i = 1 ,size(weights)
numerically = numerically + weights(i)*Func(i)*(b - a)/2
end do
end function numerically 
end module integral    
    
program test
use legendre, only : nodes_weights
use integral, only : F, analytically , numerically
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nodes(:),weights(:),Func(:)
integer               :: n
real(dp)              :: a, b
write(*,*) 'order'
read(*,*) n
write(*,*) 'left border'
read(*,*) a
write(*,*) 'right border'
read(*,*) b
call nodes_weights(n,a,b,nodes,weights)
call F(nodes,Func)
print *, nodes
print *, weights
print *, Func
print *, analytically(a,b)
print *, numerically(Func,weights,a,b)
end program test
