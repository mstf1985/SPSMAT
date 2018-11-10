function out=D_Gnbr_frac(N,theta_,derv,u_)

%
% Overview
% This function returns an operational matrix of derivative of Gegenbauer functions. 
%     
%out = D_Gnbr_frac(N,theta_,r,u_)
%
%inputs: 
%-------------------------------------------------------- 
%| N   : integer          : From  Gegenbauer sentence 0 |  
%|                            to Gegenbauer sentence N  |
%| theta_ : double          : Gegenbauer parameter      |  
%| r   : double           : Order of derivative         |
%| u_  : symbolic function: Shifting parameter          |          
%--------------------------------------------------------    
%
%Output:
%------------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational derivative Gegenbauer |       
%------------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
alpha_=theta_-0.5;
beta_= theta_-0.5;

multiplier=zeros(1,N+1);
for i=0:N
  multiplier(i+1)= factorial(i) *gamma(theta_+0.5)/gamma(i+theta_+0.5);
 end% for

 
 out=diag(1./multiplier)*D_jacobi_frac(N,alpha_,beta_,derv,u_)*diag(multiplier);



end