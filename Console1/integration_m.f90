module Integration_with_weights

    use, intrinsic :: iso_fortran_env, only: dp=>real64
    
    implicit none
    
    private
    public  :: Integral
    
    contains
    !Integration of Function with weights of the Gauss-Legendre quadrature
    pure real(dp) function Integral( Func, weights)
    
        real(dp),intent(in)    :: Func(:), weights(:)    !function values in nodes and weights for these nodes
        
        Integral = SUM(weights * Func)
        
    end function Integral
    
end module Integration_with_weights