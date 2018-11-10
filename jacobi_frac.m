function out=jacobi_frac(N,alpha_,beta_,derv,X,u_)

%

% 
% Overview
% This function returns a Shifted fractional Jacobi matrix functions. 
%     
%out = jacobi_(N,alpha_,beta_,derv,X,u_) 
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%| derv  : double           : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Jacobi functions        |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%--------------------------------------------------  
%| out   : [mx(N+1)] double : shifted fractional  |
%                             Jacobi functions    |       
%--------------------------------------------------    
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

alpha_=alpha_;
beta_=beta_;


r=diff(u, x);
r=function_handle(r);
%u_ == ax+b
a=r(0);
b=u(0);

%if the interval be over[A,B]
A=(-b-1)/a;
B=(1-b)/a;

out=zeros(length(X),N+1);
for l=1:length(X)
 out(l,1)=1;
  for i=0:N
    su=0;
    for k=0:i
    temp=((-1)^(i-k)*gamma(i+beta_+1)*gamma(i+k+alpha_+beta_+1))/...
        (gamma(k+beta_+1)*gamma(i+alpha_+beta_+1)*factorial(i-k)*factorial(k)*(B-A)^k);
      for j=0:k
        if j>= derv
        su=su+temp...
        *nchoosek(k,j)*(-1)^(k-j)*(A)^(k-j)*X(l)^(j-derv)*(gamma(j+1))/(gamma(j+1-derv));
        end%if
        
      end%for j

    end%for k
        
    out(l,i+1)=su;
  end%for i
end%l
 end%function
 