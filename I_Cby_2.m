function out=I_Cby_2(N,u_)
%
% Overview
% This function returns a operational matrix of Integral of Chebyshev(2nd) functions. 
%     
%out = I_Cby_2(N,u_) 
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
%-------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Integral Operational matrix|       
%-------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%
%
%
theta_=1;
out=I_Gnbr(N,theta_,u_);
end