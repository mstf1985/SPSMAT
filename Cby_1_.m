function out=Cby_1_(N,derv,X,u_)
% 
% Overview
% This function returns a Shifted Chebyshev (1st kind) matrix functions. 
%     
%out = Cby_1_(N,derv,X,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From Chebyshev(1st kind)|
%|                            sentence 0 to Chebyshev |
%|                            sentence number N       |
%| derv  : integer          : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Chebyshev (1st kind)    |
%|                            functions               |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Chebyshev       |
%|                            (1st kind) functions    |       
%------------------------------------------------------    
% 

alpha_=-0.5;
beta_=-0.5;

multiplier=zeros(1,length(X));
for i=0:N
  multiplier(i+1)=factorial(i) *gamma(0.5)/gamma(i+0.5);
 end% for

 out=multiplier.*jacobi_(N,alpha_,beta_,derv,X,u_);
end