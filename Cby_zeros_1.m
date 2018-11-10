function roots = Cby_zeros_1(N)
%
%overview
%Cbr_zeros_1( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Chebyshev of the first kind polynomials.
%
%inputs: 
%---------------------------------------------------------   
%| N     : integer          : Chebyshev (1st) sentence N |  
%---------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [Nx1] double : Chebyshev (1st) zeros       |       
%------------------------------------------------------    
%
%----------------------------------------------------------------       

roots=jacobi_zeros( N, -0.5, -0.5 );
    
end