function roots = Lgdr_zeros(N)
%
%
%overview
%Lgdr_zeros( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Legendre polynomials.
%
%inputs: 
%---------------------------------------------------   
%| N     : integer          : Legendre  sentence N |  
%---------------------------------------------------    
%
%Output:
%---------------------------------------------------    
%| out   : [Nx1] double     : Legendre zeros       |       
%---------------------------------------------------    
%

roots=jacobi_zeros( N, 0, 0 );
    
end