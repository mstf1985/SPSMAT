function out=Lgdr_(N,derv,X,u_)


% Overview
% This function returns a Shifted Legendre matrix functions. 
%     
%out = Lgdr_(N,derv,X,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From Legendre sentence 0|
%|                             to Legendre sentence  N|
%| derv  : integer          : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Legendre functions      |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%--------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Legendre functions|           
%--------------------------------------------------------    
% 
%

theta_=0.5;
out=Gnbr_(N,theta_,derv,X,u_);

end