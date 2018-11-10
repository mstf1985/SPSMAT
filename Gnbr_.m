function out=Gnbr_(N,theta_,derv,X,u_)
% 
% Overview
% This function returns a Shifted Gegenbauer matrix functions. 
%     
%out = Gnbr_(N,theta_,derv,X,u_) 
%
%inputs: 
%---------------------------------------------------------- 
%| N     : integer          : From  Gegenbauer sentence 0 |  
%|                            to Gegenbauer sentence N    |
%| theta_: double           : Gegenbauer parameter        |
%| derv  : integer          : derivative order            |
%|   X   : [1xm] double     : Inputs of u_(x) in          | 
%|                            Gegenbauer functions        |
%|   u_  : symbolic function: Shifting parameter          |          
%----------------------------------------------------------    
%
%Output:
%----------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Gegenbauer functions|       
%----------------------------------------------------------    
% 

alpha_=theta_-0.5;
beta_= theta_-0.5;

multiplier=zeros(1,length(X));
for i=0:N
  multiplier(i+1)=factorial(i) *gamma(theta_+0.5)/gamma(i+theta_+0.5);
 end% for

 out=multiplier.*jacobi_(N,alpha_,beta_,derv,X,u_);
end