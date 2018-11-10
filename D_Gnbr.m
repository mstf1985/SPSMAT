function out=D_Gnbr(N,theta_,u_)
%
%
%
%
% Overview
% This function returns a operational matrix of derivative of Gegenbauer functions. 
%     
%out = D_Gnbr(N,theta_,u_) 
%
%inputs: 
%---------------------------------------------------------- 
%| N     : integer          : From  Gegenbauer sentence 0 |  
%|                            to Gegenbauer sentence N    |
%|                            Gegenbauer functions        |
%| theta_ : double          : Gegenbauer parameter        |  
%| u_    : symbolic function: Shifting parameter          |          
%----------------------------------------------------------    
%
%Output:
%---------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Derivative Operational matrix|       
%---------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [a,b].
%
%

%



alpha_=theta_-0.5;
beta_= theta_-0.5;

multiplier=zeros(1,N+1);
for i=0:N
  multiplier(i+1)= factorial(i) *gamma(theta_+0.5)/gamma(i+theta_+0.5);
end% for

out=diag(1./multiplier)*D_jacobi(N,alpha_,beta_,u_)*diag(multiplier);
end