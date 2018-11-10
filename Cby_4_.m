function out=Cby_4_(N,derv,X,u_)

% Overview
% This function returns a Shifted Chebyshev (4th kind) matrix functions. 
%     
%out = Cby_4_(N,derv,X,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From Chebyshev(4th kind)|
%|                            sentence 0 to Chebyshev |
%|                            sentence number N       |
%| derv  : integer          : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Chebyshev (4th kind)    |
%|                            functions               |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Chebyshev       |
%|                            (4th kind) functions    |       
%------------------------------------------------------    
% 
%
%reference:[ORIGINAL ARTICLE On using third and fourth kinds Chebyshev
%polynomials for solving the integrated forms of high odd-order linear boundary value problems
% By E.H. Doha, W.M. Abd-Elhameed, and M.M. Alsuyuti]

alpha_=0.5;
beta_=-0.5;

multiplier=zeros(1,length(X));
for i=0:N
  multiplier(i+1)=2^(2*i)/nchoosek(2*i,i);
 end% for

 out=multiplier.*jacobi_(N,alpha_,beta_,derv,X,u_);
end