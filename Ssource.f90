module legendre
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
contains


real(dp) function Legen(n,x) result(c)
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

end function Legen


subroutine nodes_and_weights( n, nodes, weights)
integer,  intent(in   )              :: n
real(dp), intent(  out), allocatable :: nodes(:),weights(:)
real(dp), allocatable                :: nodes1(:)
real(dp)                             :: eps = 1e-12
integer                              :: i 
real(dp), parameter                  :: PI = ACOS( -1.D0 )

!first approximation

nodes = [ ( -cos( PI * ( 4.0_dp * i - 1.0_dp ) / ( 4.0_dp * n + 2.0_dp ) ), i = 1, n ) ]

allocate( nodes1( n ), source = 10.0_dp )


! cycle with Newton method for nods

do while ( abs( maxval( nodes - nodes1 ) ) >= eps )

nodes1 = nodes

do i = 1 ,n

nodes( i ) = nodes1( i ) - Legen( n,nodes1( i ) ) * ( 1.0_dp - nodes1( i ) ** 2.0_dp ) / &
    ( (Legen( n - 1, nodes1( i )) - nodes1( i ) * Legen( n, nodes1( i ) ) )* n )

end do

end do


allocate( weights( n ),source = [ ( 2_dp / ( 1_dp-nodes( i ) ** 2_dp ) / &
    ( ( ( Legen( n - 1, nodes( i ) ) - nodes( i ) * Legen( n, nodes( i ) ) ) * n / &
    ( 1_dp - nodes( i ) ** 2 ) ) **2 ), i = 1, n ) ] )

end subroutine nodes_and_weights


subroutine Poly_in_cord( Func, cord)
real(dp), intent(inout) :: Func(:,:)
real(dp), intent(in   ) :: Cord(:)
integer                 :: i,j

do i = 1,SIZE(Func, 1)

do j = 1,SIZE(Func, 2)

Func(i,j) = Legen(i - 1, Cord(j) )

end do

end do

end subroutine Poly_in_cord


real(dp) function Integral( Func, weights)
real(dp),intent(in   ) :: weights(:)
real(dp),intent(in   ) :: Func(:)

Integral = sum( weights * Func )

end function Integral


subroutine Ssource(Flux, Cross_sec_scat, order, energy, res) 
real(dp),              intent(in   ) :: Flux(:,:,:), Cross_sec_scat(:,:,:,:)
real(dp), allocatable                :: Poly(:,:), nodes(:), weights(:), Integr(:,:,:)
real(dp), allocatable, intent(inout) :: res(:,:,:)
integer,               intent(in   ) :: order, energy
integer                              :: i, j, k 
!Flux(ncord, nangles, nenergy)
!Cross_sec_scat(ncord,nenergy (out),nenergy (in),order)

allocate( res( SIZE( Flux, 1), SIZE( Flux, 2), SIZE( Flux, 3) ) )
res = 0

allocate( Poly( order + 1 , SIZE( Flux, 2) ) )

call nodes_and_weights( SIZE( Flux, 2), nodes, weights)

call Poly_in_cord( Poly, nodes )

allocate( Integr( SIZE( Flux, 1), order + 1, SIZE( Flux, 3)) )

do i = 1, SIZE( Flux, 1)

do j = 1, order + 1

do k = 1, SIZE( Flux, 3)

Integr(i,j,k) = Integral( Poly(j,:) * Flux(i, :, k ), weights)

res(i,:,k) = res(i,:,k) +  Poly(j,:) * Cross_sec_scat(i, k, energy, j) * ( 2_dp * ( j - 1_dp ) + 1_dp) / 2_dp * Integr(i,j,k)

end do

end do

end do

end subroutine Ssource

end module legendre



program test
use legendre, only : Ssource, nodes_and_weights, Poly_in_cord
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nodes(:), weights(:), Flux(:,:,:), Cross_sec_scat(:,:,:,:), SS(:,:,:), Poly(:,:), SS_an(:,:,:)
real(dp)              :: integr_analytic(6)
integer               :: order,ncord,nenergy,nangles, i, j, k,energy

order = 5
ncord = 3
nenergy = 5
energy = 1
nangles = 5

call nodes_and_weights(nangles,nodes,weights)

allocate(Flux(ncord, nangles, nenergy))

do i = 1, SIZE( Flux, 1)

do k = 1, SIZE( Flux, 3)

Flux(i,:,k) = cos( nodes )

end do

end do



allocate( Cross_sec_scat(ncord, nenergy, nenergy, order + 1), source = 1.0_dp )



integr_analytic(1) = 2_dp * sin(1.0_dp)
integr_analytic(2) = 0_dp
integr_analytic(3) = 6_dp * cos(1.0_dp) - 4_dp * sin(1.0_dp)
integr_analytic(4) = 0_dp
integr_analytic(5) = 122_dp * sin(1.0_dp) - 190_dp * cos(1.0_dp)
integr_analytic(6) = 0_dp


allocate( Poly(order + 1, nangles) )

call Poly_in_cord( Poly, nodes)

call Ssource( Flux, Cross_sec_scat, order, energy, SS)

allocate( SS_an( ncord, nangles, nenergy ), source = 0.0_dp)

do i = 1, SIZE( Flux, 1)

do j = 1, order + 1

do k = 1, SIZE( Flux, 3)


SS_an(i,:,k) = SS_an(i,:,k) + Poly(j,:) * Cross_sec_scat(i, k, energy, j) * ( 2_dp * ( j - 1_dp ) + 1_dp) / 2_dp * integr_analytic(j)

end do

end do

end do

write(*,*) maxval( abs(SS_an - SS) )
write(*,*) maxval( abs( (SS_an - SS) / SS) )


deallocate(Flux)
deallocate(Cross_sec_scat)
deallocate(SS_an)
deallocate(nodes)
deallocate(Poly)
deallocate(SS)
deallocate(weights)

end program test
