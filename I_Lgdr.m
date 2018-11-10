function out=I_Lgdr(N,u_)
%
% Overview
% This function returns a operational matrix of Integral of Legendre functions. 
%     
%out = I_Lgdr(N,u_)
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
%| out   : [(N+1)x(N+1)] double : Integral Operational matrix|       
%---------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%
%
%
theta_=0.5;
out=I_Gnbr(N,theta_,u_);

end