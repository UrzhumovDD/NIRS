module Poly_in_nodes
    
    use legendre, only : Pn
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    
    implicit none
    
    private
    public  :: All_poly
    
    contains
    
    !calculating polynomials of all orders from 0 to SIZE(Poly, 1) - 1 in nodes
    subroutine All_poly( Poly, nodes)
    
        real(dp), intent(inout)     :: Poly(:,:)    ! values of legendre polynomials in nodes
        real(dp), intent(in)        :: nodes(:)     ! values of nodes
        integer                     :: i,j          ! loop counters
    
        do i = 1,SIZE(Poly, 1) ! orders of polynomials
            do j = 1,SIZE(Poly, 2) 
                Poly(i,j) = Pn(i - 1, nodes(j) )
            end do
        end do
        
    end subroutine All_poly
    
end module Poly_in_nodes