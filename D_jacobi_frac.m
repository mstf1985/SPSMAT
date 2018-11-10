function out=D_jacobi_frac(n,alpha_,beta_,r,u_)
%
% Overview
% This function returns a operational matrix of derivative of Jacobi functions. 
%     
%out = D_jacobi_frac(N,alpha_,beta_,r,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%| r     : double           : Order of derivative     |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%----------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Operational derivative Jacobi |       
%----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);

% now calculate the derivative of the function

g=diff(u, x);
g=function_handle(g);

%u_ == Ax+B
A=g(0);
B=u(0);

%[a,b]
a=(-B-1)/A;
b=(1-B)/A;


out=zeros(n+1,n+1);

for i=0:n
 for j=0:n
 sum_=0;
 for k=ceil(r):i
 for l=0:j 
 
  sum_=sum_+ (A**(alpha_+beta_+1)/2**(alpha_+beta_+1))*((-1)**(i-k)*(b-a)**(1+alpha_+beta_-r)...
  *gamma(j+beta_+1)*...
  gamma(i+beta_+1)*gamma(k+alpha_+beta_+1+i)*(-1)**(j-l)*...
  gamma(j+alpha_+beta_+1+l)*gamma(alpha_+1)*gamma(k-r+beta_+1+l))/...
  (v(j,alpha_,beta_)*gamma(k+beta_+1)*gamma(j+alpha_+beta_+1)...
  *gamma(alpha_+beta_+1+i)*gamma(k+1-r)*gamma(i-k+1)*gamma(l+1)*...
   gamma(j-l+1)*gamma(l+beta_+1)*gamma(l+k-r+alpha_+beta_+2));

  end%for l
  end%for k
  out(j+1,i+1)=sum_; %this j+1 and i+1 are exchanged intentionally.
 end%for j
end%for i




end%function


%------------------------------------------------------------

function out=v(N, alpha_,beta_)

out=(gamma(N+alpha_+1)*gamma(N+beta_+1))/...
(gamma(N+1)*gamma(N+alpha_+beta_+1)*(2*N+alpha_+beta_+1));


end%func