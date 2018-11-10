function out=lgr_Jacobi_(N,X,points, u_)
%
%% 
% Overview
% This function returns a Shifted Lagrange Jacobi matrix functions. 
%     
%out = lgr_Jacobi_(N,X,points, u_)
%
%inputs: 
%------------------------------------------------------ 
%|   N   : integer          : N+1 sentences are       |  
%|                            considered              |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Lagrange functions      |
%|  pints: [1xN] double     : Inputs for making       |
%|                            Lagrange polynomilas    |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------    
%| out   : [mxN] double : shifted Lagrange Jacobi functions|       
%-----------------------------------------------------------    
                         


%The difference between this lgr_Jacobi_(N,X,points, u_) and lgr_Jacobi(N,X,points, u_) is that 
% this function recieve an array and return a matrix but the other one is recieveing a scalar
% and returns  a vector.
%
%
                         
warning off;
f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);
                         
 
   %initialization
   L=ones(length(X),N); 
 
    for j=1:1:N
   
     for i=1:1:N
      if i!=j
        for k=1:1:length(X)
         L(k,j)=L(k,j)*(u(X(k))-u(points(i)))/(u(points(j))-u(points(i)));
        end
      end%if
     end%for 
   end%for

   out=L;
end