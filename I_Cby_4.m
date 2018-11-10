function out=I_Cby_4(N,u_)
%
% Overview
% This function returns a operational matrix of Integral of Chebyshev(4th) functions. 
%     
%out = I_Cby_4(N,u_) 
%
%inputs: 
%-------------------------------------------------------------- 
%| N     : integer          : From  Chebyshev(4th) sentence 0 |  
%|                            to Chebyshev(4th) sentence N    |
%|                            Chebyshev(4th) functions        |
%| u_    : symbolic function: Shifting parameter              |          
%--------------------------------------------------------------    
%
%Output:
%--------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Integral Operational matrix |       
%--------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%
%
%

alpha_=0.5;
beta_=-0.5;

multiplier=zeros(1,N+1);
for i=0:N
  multiplier(i+1)=2^(2*i)/nchoosek(2*i,i);
 end% for

 out=diag(1./multiplier)*I_jacobi(N,alpha_,beta_,u_)*diag(multiplier);

 end