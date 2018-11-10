function out=I_jacobi_frac(N,alpha_,beta_,r,u_)

%
% Overview
% This function returns a operational matrix of Integral of Jacobi functions. 
%     
%out = I_jacobi_frac(N,alpha_,beta_,r,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%| r     : double           : Order of integral       |
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


% Reference: Kazem S. An integral operational matrix based on
% Jacobi polynomials for solving fractional-order differential
% equations Applied Mathematical Modelling.2013. 37 (3), 1126-1136 

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
  for k=0:i
   for l=0:j
     su=su+p(i,k,alpha_,beta_)*p(j,l,alpha_,beta_)*...
     (gamma(k+1)*B_(k+l+r+beta_+1,alpha_+1))/(gamma(k+r+1)*v(j,alpha_,beta_));
  

   end% for l
 end% for k 
 I(i+1,j+1)=su;
end%for j

end%for i
out=b*I;


end%func

%--------------------------------------
function out=v(N, alpha_,beta_)

out=(gamma(N+alpha_+1)*gamma(N+beta_+1))/...
(gamma(N+1)*gamma(N+alpha_+beta_+1)*(2*N+alpha_+beta_+1));


end%func
%--------------------------------------

function out=B_(a,b)

out=(gamma(a)*gamma(b))/(gamma(a+b));

end
%--------------------------------------
function out=p(n,i,alpha_,beta_)

out=((-1)^(n-i))*(gamma(n+alpha_+beta_+i+1)./(gamma(n+alpha_+beta_+1)*gamma(i+1)))...
*(gamma(n+alpha_+1)./(gamma(alpha_+i+1)*gamma(n-i+1)));


end































