function out=I_Cby_frac_2(N,r,u_)
%
%
% Overview
% This function returns an operational matrix of Integral of Chebyshev(2nd) functions. 
%     
%out = I_Cby_frac_2(N,r,u_)
%
%inputs: 
%------------------------------------------------------------ 
%| N   : integer          : From  Chebyshev(2nd) sentence 0 |  
%|                            to Chebyshev(2nd) sentence N  |
%| r   : double           : Order of integral               |
%| u_  : symbolic function: Shifting parameter              |          
%------------------------------------------------------------    
%
%Output:
%---------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational Integral Chebyshev |       
%---------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%


theta_=1;
out=I_Gnbr_frac(N,theta_,r,u_);

end