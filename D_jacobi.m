function out=D_jacobi(n,alpha_,beta_,u_)

%
% Overview
% This function returns a operational matrix of derivative of Jacobi functions. 
%     
%out = D_jacobi(n,alpha_,beta_,u_)
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

ffd = diff(u, x);
fffd=diff(ffd,x);

u_derv = function_handle(ffd);

a_=alpha_;
b_=beta_;


r=diff(u, x);
r=function_handle(r);

%u_ = Ax+B
A=r(0);
B=u(0);

%[a,b]
a=(-B-1)/A;
b=(1-B)/A;


E=zeros(n,n);
for i=1:1:n
 for j=1:1:n
  if i==j
    E(i,j)=2*(a_+b_+i)/((a_+b_+2*i-1)*(a_+b_+2*i));
   elseif i==j-1
    E(i,j)=2*(a_-b_)/((a_+b_+2*i)*(a_+b_+2*i+2));
   elseif i==j-2
    E(i,j)=-2*(i+a_+1)*(i+b_+1)/((a_+b_+i+1)*(a_+b_+2*i+2)*(a_+b_+2*i+3));
  endif
 end%for
end%for




out=[zeros(n+1,1),[inv(E);zeros(1,n)]];
out=2*out/(b-a);
end%function

