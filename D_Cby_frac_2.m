function out=D_Cby_frac_2(N,derv,u_)
%
%
%
% Overview
% This function returns an operational matrix of derivative of Chebyshev(2nd) functions. 
%     
%out = D_Cby_frac_2(N,r,u_)
%
%inputs: 
%------------------------------------------------------------ 
%| N   : integer          : From  Chebyshev(2nd) sentence 0 |  
%|                            to Chebyshev(2nd) sentence N  |
%| r   : double           : Order of derivative             |
%| u_  : symbolic function: Shifting parameter              |          
%------------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational derivative Chebyshev |       
%-----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%

theta_=1;
out=D_Gnbr_frac(N,theta_,derv,u_);

end