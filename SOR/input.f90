module input_mod
implicit none

integer::iter_option=3! 1:Jacobi, 2:Gauss-Seidel 3:SOR(optimal) 
integer,parameter::N=64
real,parameter::tol=1.e-15
real,parameter::pi=3.141592653589793238462643383


end module input_mod
