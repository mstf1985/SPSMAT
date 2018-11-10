function out=D_Cby_2(N,u_)
%
%
%
% Overview
% This function returns a operational matrix of derivative of Chebyshev(2nd) functions. 
%     
%out = D_Cby_2(N,u_) 
%
%inputs: 
%-------------------------------------------------------------- 
%| N     : integer          : From  Chebyshev(2nd) sentence 0 |  
%|                            to Chebyshev(2nd) sentence N    |
%|                            Chebyshev(2nd) functions        |
%| u_    : symbolic function: Shifting parameter              |          
%--------------------------------------------------------------    
%
%Output:
%---------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : derivative Operational matrix|       
%---------------------------------------------------------------    
% 
%
-% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
%
theta_=1;
out=D_Gnbr(N,theta_,u_);
end