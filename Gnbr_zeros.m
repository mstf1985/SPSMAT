function roots = Gnbr_zeros( N, tetha_ )
%
%
%overview
%Gnbr_zeros( N ) is a function returing N roots of the N-th sentence of 
% Standard orthogonal Gegenbauer polynomials.
%
%inputs: 
%-----------------------------------------------------   
%| N     : integer          : Gegenbauer  sentence N |  
%| theta_: double           : Gegenbauer parameter   |
%-----------------------------------------------------    
%
%Output:
%-----------------------------------------------------    
%| out   : [Nx1] double     : Gegenbauer zeros       |       
%-----------------------------------------------------    
%

roots=jacobi_zeros( N, tetha_-0.5, tetha_-0.5 );
    
end