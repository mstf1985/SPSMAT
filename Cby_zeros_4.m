function roots = Cby_zeros_4(N)
%
%overview
%Cbr_zeros_4( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Chebyshev of the forth kind polynomials.
%
%inputs: 
%-----------------------------------------------------   
%| N     : integer      : Chebyshev (4th) sentence N |  
%-----------------------------------------------------    
%
%Output:
%-----------------------------------------------------    
%| out   : [Nx1] double : Chebyshev (4th) zeros      |       
%-----------------------------------------------------    
%

roots=jacobi_zeros( N, -0.5, 0.5 );
    
end