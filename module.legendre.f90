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
subroutine nods_weights(n,a,b,nods,weights) 
integer, intent(in), optional     :: n,a,b
real(dp)                          :: eps = 1e-8
integer                           :: i 
real, parameter                   :: PI = DACOS(-1.D0)
real(dp),allocatable,intent(out)  :: nods(:) , weights(:)
real(dp),allocatable              :: nods1(:)
!first approximation
allocate( nods(n) )
nods = [( -cos( PI*(4*i - 1)/(4*n + 2) ), i = 1 , n )]
allocate( nods1(n) )
nods1 = 10
! cycle with Newton method for nods
do while (abs(maxval(nods - nods1)) >= eps)   
nods1 = nods
do i = 1 , n
nods(i) = nods1(i) - lezh(n,nods1(i)) * (1 - nods1(i)**2) / &
    ( ( lezh(n - 1,nods1(i)) - nods1(i) * lezh(n,nods1(i)) )* n )
end do 
end do
!weights
allocate( weights(n),source = [(2/(1-nods(i)**2)/ &
    ( ( (lezh(n - 1,nods(i))-nods(i) * lezh(n,nods(i)) )* n / (1-nods(i)**2))**2),i = 1, n)] )
!affine transformations
nods = (nods * 2 -(b + a) )/(b - a)
end subroutine nods_weights  
end module legendre
    
program test
use legendre, only : nods_weights
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nods(:),weights(:)
integer               :: n, a, b
write(*,*) 'order'
read(*,*) n
write(*,*) 'left border'
read(*,*) a
write(*,*) 'right border'
read(*,*) b
call nods_weights(n,a,b,nods,weights)

print *, nods
print *, weights
end program test
