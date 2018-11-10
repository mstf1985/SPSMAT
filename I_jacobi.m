function out=I_jacobi(N,alpha_,beta_,u_)
%
% Overview
% This function returns a operational matrix of Integral of Jacobi functions. 
%     
%out = I_jacobi(N,alpha_,beta_,u_) 
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%----------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Operational Integral Jacobi |       
%----------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [0,b].
%
%Reference: A generalized operational method for solving 
%integro–partial differential equations based on Jacobi 
%polynomials Abdollah Borhanifar, and Khadijeh Sadri

f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);


g=diff(u, x);
g=function_handle(g);

%u_ == Ax+B
A=g(0);
B=u(0);

%[a,b]
a=(-B-1)/A;
% a must be 0
b=(1-B)/A;


I=zeros(N+1,N+1);
for i=0:N
for j=0:N
su=0;
for l=0:i
su=su+w(i,j,l,alpha_,beta_);
end% for l
I(i+1,j+1)=b*su;
end%for j

end%for i
out=I;
end%func

%--------------------------------------
function out=v(N, alpha_,beta_)

out=(gamma(N+alpha_+1)*gamma(N+beta_+1))/...
(gamma(N+1)*gamma(N+alpha_+beta_+1)*(2*N+alpha_+beta_+1));


end%func
%--------------------------------------

function out=w(i,j,l,alpha_,beta_)
  
  
  mul=((-1)^(i-l)*gamma(i+beta_+1)*gamma(i+l+alpha_+beta_+1))...
  /(gamma(l+beta_+1)*gamma(i+alpha_+beta_+1)*gamma(l+2)*gamma(i-l+1));

  su=0;
for k=0:j
su=su+((-1)^(j-k)*gamma(j+k+alpha_+beta_+1)*gamma(j+beta_+1)*gamma(k+l+beta_+2)*gamma(alpha_+1))...
/(v(j,alpha_,beta_)*gamma(k+beta_+1)*gamma(j+alpha_+beta_+1)*gamma(k+l+alpha_+beta_+3)*gamma(k+1)*gamma(j-k+1));
end%for
out=su*mul;
end%func


