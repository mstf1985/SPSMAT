function out=Lgdr_frac(N,derv,X,u_)
%
% 
% Overview
% This function returns a Shifted fractional Legendre matrix functions. 
%     
%out = Lgdr_frac(N,theta_,derv,X,u_)
%
%inputs: 
%---------------------------------------------------------------- 
%| N     : integer          : From fractional Legendre sentence |
%|                               0 to fractional Legendre       |
%| derv  : double           : derivative order                  |
%|   X   : [1xm] double     : Inputs of u_(x) in  fractional    | 
%|                            Legendre functions                |
%|   u_  : symbolic function: Shifting parameter                |           
%----------------------------------------------------------------    
%
%Output:
%------------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted fractional Legendre   |
%|                             functions                    |       
%------------------------------------------------------------    
% 
%
theta_=0.5;
out=Gnbr_frac(N,theta_,derv,X,u_);

end