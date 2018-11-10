function out=I_Lgdr_frac(N,r,u_)
%
%
%
% Overview
% This function returns an operational matrix of Integral of Legendre functions. 
%     
%out = Lgdr_frac(N,r,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N   : integer          : From  Legendre sentence 0 |  
%|                            to Legendre sentence N  |
%| r   : double           : Order of integral         |
%| u_  : symbolic function: Shifting parameter        |          
%------------------------------------------------------    
%
%Output:
%---------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational Integral  Legendre |       
%---------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%

theta_=0.5;
out=I_Gnbr_frac(N,theta_,r,u_);
end