function out=D_Cby_3(N,u_)
%
% Overview
% This function returns a operational matrix of derivative of Chebyshev(3rd) functions. 
%     
%out = I_Cby_3(N,u_) 
%
%inputs: 
%-------------------------------------------------------------- 
%| N     : integer          : From  Chebyshev(3rd) sentence 0 |  
%|                            to Chebyshev(3rd) sentence N    |
%|                            Chebyshev(3rd) functions        |
%| u_    : symbolic function: Shifting parameter              |          
%--------------------------------------------------------------    
%
%Output:
%-------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Derivative Operational matrix|       
%-------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
alpha_=-0.5;
beta_=0.5;

multiplier=zeros(1,N+1);
for i=0:N
  multiplier(i+1)=2^(2*i)/nchoosek(2*i,i);
 end% for

 out=diag(1./multiplier)*D_jacobi(N,alpha_,beta_,u_)*diag(multiplier);

end