function out=D_Lgdr_frac(N,derv,u_)
%
%
%
% Overview
% This function returns an operational matrix of derivative of Legendre functions. 
%     
%out = D_Lgdr_frac(N,derv,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N   : integer          : From  Legendre sentence 0 |  
%|                            to Legendre sentence N  |
%| r   : double           : Order of derivative       |
%| u_  : symbolic function: Shifting parameter        |          
%------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational derivative  Legendre |       
%-----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
theta_=0.5;
out=D_Gnbr_frac(N,theta_,derv,u_);
end