function out=jacobi_(N,alpha_,beta_,derv,X,u_)

% 
% Overview
% This function returns a Shifted Jacobi matrix functions. 
%     
%out = jacobi_(N,alpha_,beta_,derv,X,u_) 
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From  Jacobi sentence 0 |  
%|                            to Jacobi sentence N    |
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%| derv  : integer          : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Jacobi functions        |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Jacobi functions|       
%------------------------------------------------------    
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

alpha_=alpha_+derv;
beta_=beta_+derv;

for i=1:length(X)
 if N<0 
   out =0 ;
   return
 elseif N==0 
   out(i,:)=1;
   return   
 else 
   Ja(1)=1;
   Ja(2)=(0.5)*(alpha_+beta_+2)*u(X(i))+(0.5)*(alpha_-beta_);
   for n=1:1:N-1
     a=((2*n+alpha_+beta_+1)*(2*n+alpha_+beta_+2))/((2*(n+1)*(n+alpha_+beta_+1)));
     b=((alpha_^2-beta_^2)*(2*n+alpha_+beta_+1))/(2*(n+1)*(n+alpha_+beta_+1)*(2*n+alpha_+beta_));
     c=((n+alpha_)*(n+beta_)*(2*n+alpha_+beta_+2))/((n+1)*(n+alpha_+beta_+1)*(2*n+alpha_+beta_));
     Ja(n+2)=((a)*u(X(i))+(b))*Ja(n+1)-(c)*Ja(n);
   end%for

 end%else
 
 out(i,:)=Ja;
 if derv~=0 
 for n=1:1:N+1
   if (n-derv>0)
    Jaa(n)=(u_derv(X(i)))^derv*((gamma((alpha_-derv)+(beta_-derv)+(n-1)+1+derv))/(2^derv*gamma((alpha_-derv)+(beta_-derv)+(n-1)+1)))*Ja(n-derv);
   else 
    Jaa(n)=0;
   end %if
 end%for 
 out(i,:)=Jaa; 
 end %if
end %for X
 
end%function
  
 