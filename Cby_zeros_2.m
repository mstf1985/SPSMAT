function roots = Cby_zeros_2(N)

%overview
%Cbr_zeros_2( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Chebyshev of the second kind polynomials.
%
%inputs: 
%---------------------------------------------------------   
%| N     : integer          : Chebyshev (2nd) sentence N |  
%---------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [Nx1] double : Chebyshev (2nd) zeros       |       
%------------------------------------------------------    
%

 
%----------------------------------------------------------------       

roots=jacobi_zeros( N, 0.5, 0.5 );
    
end