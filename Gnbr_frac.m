function out=Gnbr_frac(N,theta_,derv,X,u_)

% 
% Overview
% This function returns a Shifted fractional Gegenbauer matrix functions. 
%     
%out = Gnbr_frac(N,theta_,derv,X,u_)
%
%inputs: 
%------------------------------------------------------------------ 
%| N     : integer          : From fractional Gegenbauer sentence |
%|                               0 to fractional Gegenbauer       |
%|
%| theta_: double           : Gegenbauer parameter                |
%| derv  : double           : derivative order                    |
%|   X   : [1xm] double     : Inputs of u_(x) in                  | 
%|                            fractional Gegenbauer functions     |
%|                            functions                           |
%|   u_  : symbolic function: Shifting parameter                  |           
%------------------------------------------------------------------    
%
%Output:
%------------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted fractional Gegenbauer |
%|                             functions                    |       
%------------------------------------------------------------    
% 

alpha_=theta_-0.5;
beta_=theta_-0.5;

multiplier=zeros(1,length(X));
for i=0:N
  multiplier(i+1)=factorial(i) *gamma(theta_+0.5)/gamma(i+theta_+0.5);
 end% for

 out=multiplier.*jacobi_frac(N,alpha_,beta_,derv,X,u_);
 
 end
 
 