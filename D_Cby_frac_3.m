function out=D_Cby_frac_3(N,derv,u_)
%
%
% Overview
% This function returns an operational matrix of derivative of Chebyshev(3rd) functions. 
%     
%out = D_Cby_frac_3(N,r,u_)
%
%inputs: 
%------------------------------------------------------------ 
%| N   : integer          : From  Chebyshev(3rd) sentence 0 |  
%|                            to Chebyshev(3rd) sentence N  |
%| r   : double           : Order of derivative             |
%| u_  : symbolic function: Shifting parameter              |          
%------------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational derivative Chebyshev |       
%-----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%

alpha_=-0.5;
beta_=0.5;

multiplier=zeros(1,N+1);
for i=0:N
  multiplier(i+1)=2^(2*i)/nchoosek(2*i,i);
 end% for

 out=diag(1./multiplier)*D_jacobi_frac(N,alpha_,beta_,derv,u_)*diag(multiplier);

end