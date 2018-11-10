function out=D_Lgdr(N,u_)
%
%
%
% Overview
% This function returns a operational matrix of derivative of Legendre functions. 
%     
%out = D_Lgdr(N,u_)
%
%inputs: 
%-------------------------------------------------------- 
%| N     : integer          : From  Legendre sentence 0 |  
%|                            to Legendre sentence N    |
%|                            Legendre functions        | 
%| u_    : symbolic function: Shifting parameter        |          
%--------------------------------------------------------    
%
%Output:
%---------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Derivative Operational matrix|       
%---------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
%
theta_=0.5;
out=D_Gnbr(N,theta_,u_);

end