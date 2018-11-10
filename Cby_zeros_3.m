function roots = Cby_zeros_3(N)
%
%overview
%Cbr_zeros_3( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Chebyshev of the 3rd kind polynomials.
%
%inputs: 
%---------------------------------------------------------   
%| N     : integer          : Chebyshev (3rd) sentence N |  
%---------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [Nx1] double : Chebyshev (3rd) zeros       |       
%------------------------------------------------------    
%


roots=jacobi_zeros( N, 0.5, -0.5 );
    
end