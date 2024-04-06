module Legendre
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
!Calculating legendre's polynomial of the n order in x
real(dp) function Pn(n,x) result(c)
real(dp), intent(in   )     :: x
integer,  intent(in   )     :: n
real(dp)                    :: a, b
integer                     :: i
!
if ( n == 0 ) then
    c = 1.0_dp
!
else if ( n == 1 ) then
    c = x
!
else
    a = 1.0_dp
    b = x
!
do i = 2, n
!
    c = b * x * ( 2.0_dp*i - 1.0_dp ) / i - a * ( i - 1.0_dp ) / i
    a = b
    b = c
!
end do
end if
end function Pn
end module Legendre
!
module Calc_of_nodes
use legendre, only : Pn
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
subroutine C_o_n( n, nodes)
integer,  intent(in   )              :: n
real(dp), intent(  out), allocatable :: nodes(:)
real(dp), allocatable                :: nodes1(:)
real(dp)                             :: eps = 1e-12
integer                              :: i 
real(dp), parameter                  :: PI = ACOS( -1.D0 )
!first approximation
nodes = [ ( -COS( PI * ( 4.0_dp * i - 1.0_dp ) / ( 4.0_dp * n + 2.0_dp ) ), i = 1, n ) ]
allocate( nodes1( n ), source = 10.0_dp )
! cycle with Newton method for nodes
do while ( ABS(MAXVAL( nodes - nodes1 ) ) >= eps)
    nodes1 = nodes
    do i = 1 ,n        
        nodes( i ) = nodes1( i ) - Pn( n,nodes1( i ) ) * ( 1.0_dp - nodes1( i ) ** 2.0_dp ) / &
            ( (Pn( n - 1, nodes1( i )) - nodes1( i ) * Pn( n, nodes1( i ) ) )* n )       
    end do
end do
end subroutine C_o_n
end module Calc_of_nodes
!
module Calc_of_weights
use legendre, only : Pn
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
subroutine C_o_w(nodes,weights)
real(dp), intent(  out), allocatable :: weights(:)
real(dp), intent(in   )              :: nodes(:)
integer                              :: i,n
!
n = SIZE( nodes )
allocate( weights( n ),source = [ ( 2_dp / ( 1_dp-nodes( i ) ** 2_dp ) / &
    ( ( ( Pn( n - 1, nodes( i ) ) - nodes( i ) * Pn( n, nodes( i ) ) ) * n / &
    ( 1_dp - nodes( i ) ** 2 ) ) **2 ), i = 1, n ) ] )
end subroutine C_o_w
end module Calc_of_weights
!
module Poly_in_nodes
use legendre, only : Pn
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
!calculating polynomials of all orders from 0 to SIZE(Poly, 1) - 1 in nodes
subroutine All_poly( Poly, nodes)
real(dp), intent(inout) :: Poly(:,:)
real(dp), intent(in   ) :: nodes(:)
integer                 :: i,j
!
do i = 1,SIZE(Poly, 1)
    do j = 1,SIZE(Poly, 2)

        Poly(i,j) = Pn(i - 1, nodes(j) )

    end do
end do
!
end subroutine All_poly
end module Poly_in_nodes
!
module Integration_with_weights
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
!
real(dp) function Integral( Func, weights)
real(dp),intent(in   ) :: weights(:)
real(dp),intent(in   ) :: Func(:)
!
Integral = SUM( weights * Func )
!
end function Integral
end module Integration_with_weights
!
module Scaterring_source
use Calc_of_nodes,            only : C_o_n
use Calc_of_weights,            only : C_o_w
use Poly_in_nodes,            only : All_poly
use Integration_with_weights, only : Integral
use, intrinsic :: iso_fortran_env, only: dp=>real64
implicit none
    contains
subroutine Ssource(Flux, Cross_sec_scat, order, energy, res) 
real(dp),              intent(in   ) :: Flux(:,:,:), Cross_sec_scat(:,:,:,:)
real(dp), allocatable                :: Poly(:,:), nodes(:), weights(:), Integr(:,:,:)
real(dp), allocatable, intent(inout) :: res(:,:,:)
integer,               intent(in   ) :: order, energy
integer                              :: i, j, k 
!Flux(ncord, nangles, nenergy)
!Cross_sec_scat(ncord,nenergy (out),nenergy (in),order)
!res(result) is final value of Ssource
allocate( res( SIZE(Flux, 1), SIZE(Flux, 2), SIZE(Flux, 3) ), source = 0.0_dp )
allocate( Poly( order + 1 , SIZE(Flux, 2) ) )
call C_o_n( SIZE(Flux, 2), nodes)
call C_o_w( nodes, weights)
call All_poly( Poly, nodes )
allocate( Integr( SIZE(Flux, 1), order + 1, SIZE(Flux, 3)) )
!
do i = 1, SIZE(Flux, 1)
    do j = 1, order + 1
        do k = 1, SIZE(Flux, 3)
            Integr(i,j,k) = Integral( Poly(j,:) * Flux(i, :, k ), weights)
            res(i,:,k) = res(i,:,k) +  Poly(j,:) * Cross_sec_scat(i, k, energy, j) * ( 2_dp * ( j - 1_dp ) + 1_dp) / 2_dp * Integr(i,j,k)
        end do
    end do
end do
!
end subroutine Ssource
end module Scaterring_source
!
program test
use Scaterring_source, only : Ssource
use Calc_of_nodes, only     : C_o_n
use Poly_in_nodes, only     : All_poly
use, intrinsic        :: iso_fortran_env, only: dp=>real64
implicit none 
real(dp), allocatable :: nodes(:), Flux(:,:,:), Cross_sec_scat(:,:,:,:), SS(:,:,:), Poly(:,:), SS_an(:,:,:)
real(dp)              :: integr_analytic(6)
integer               :: order,ncord,nenergy,nangles, i, j, k,energy
!
order = 5
ncord = 3
nenergy = 5
energy = 1
nangles = 5
!
call C_o_n(nangles,nodes)
allocate(Flux(ncord, nangles, nenergy))
!
do i = 1,SIZE(Flux,1)
    do k = 1,SIZE(Flux,3)

        Flux(i,:,k) = COS(nodes)

    end do
end do
!
allocate(Cross_sec_scat(ncord,nenergy,nenergy,order + 1),source = 1.0_dp)
!analyticaly calculated integral of Pn*cos(mu)
integr_analytic(1) = 2_dp * SIN(1.0_dp)
integr_analytic(2) = 0_dp
integr_analytic(3) = 6_dp * COS(1.0_dp) - 4_dp * SIN(1.0_dp)
integr_analytic(4) = 0_dp
integr_analytic(5) = 122_dp * SIN(1.0_dp) - 190_dp * COS(1.0_dp)
integr_analytic(6) = 0_dp
!
allocate( Poly(order + 1, nangles) )
call All_poly( Poly, nodes)
call Ssource( Flux, Cross_sec_scat, order, energy, SS)
allocate( SS_an(ncord, nangles, nenergy ), source = 0.0_dp)
!
do i = 1, SIZE(Flux, 1)
    do j = 1, order + 1
        do k = 1, SIZE(Flux, 3)
            SS_an(i,:,k) = SS_an(i,:,k) + Poly(j,:) * Cross_sec_scat(i, k, energy, j) * ( 2_dp * ( j - 1_dp ) + 1_dp) / 2_dp * integr_analytic(j)
        end do
    end do
end do
!
!absolute difference
write(*,*) MAXVAL(ABS( SS_an - SS ) )
!relative difference
write(*,*) MAXVAL(ABS( (SS_an - SS) / SS ) )
!
deallocate(Flux)
deallocate(Cross_sec_scat)
deallocate(SS_an)
deallocate(nodes)
deallocate(Poly)
deallocate(SS)
end program test
