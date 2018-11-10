function out=G_jacobi(n,alpha_,beta_,u_)
% 
% Overview
% This function returns a the G matrix disscussed in the article. 
%
% out = G_jacobi(n,alpha_,beta_,u_)
%
%inputs: 
%------------------------------------------------------ 
%|   N   : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%--------------------------------------------    
%| out   : [(N+1)x(N+1)] double : G matrix  |       
%--------------------------------------------    
% 
%
a_=alpha_;
b_=beta_;

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

G=zeros(n+1,n+1);

for i=1:1:n+1
 for j=1:1:n+1
  if i==j
    G(i,j)=(b_^2-a_^2)/((a_+b_+2*i+1)*(a_+b_+2*i+2));
   elseif i==j-1
    G(i,j)=2*(i+a_)*(i+b_)/((a_+b_+2*i)*(a_+b_+2*i+1));
   elseif i==j+1
    G(i,j)=2*(i)*(i+b_+a_)/((a_+b_+2*i)*(a_+b_+2*i-1));
  endif
 end%for
end%for


out=2*G/(b-a);
end%function