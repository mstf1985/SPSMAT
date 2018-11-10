function out=lgr_Jacobi(N,X,points, u_)
%
% Overview
% This function returns a Shifted Lagrange Jacobi matrix functions. 
%     
%out = lgr_Jacobi(N,z,points, u_)
%
%inputs: 
%------------------------------------------------------ 
%|   N   : integer          : N+1 sentences are       |  
%|                            considered              |
%|   z   : double           : Inputs of u_(x) in      | 
%|                            Lagrange functions      |
%|  pints: [1xN] double     : Inputs for making       |
%|                            Lagrange polynomilas    |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------   
%Output:
%-----------------------------------------------------------    
%| out   : [1xN] double : shifted Lagrange Jacobi functions|       
%-----------------------------------------------------------   
                         
%Rule #1=The length of points must be N
%Rule #2= X is only one scalar number (not vecotr)
%Rule #3=There is no difference between guass/gaussRadau/gaussLobatto points
#Rule #4= points  must be N roots of jacobi polynomials(either guass/gaussRadau/gaussLobatto)
                         

                         
warning off;
f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);
                         
   L=zeros(1,N); 
 
    for j=1:1:N
   L(j)=1;
     for i=1:1:N
      if i!=j
        L(j)=L(j)*(u(X)-u(points(i)))/(u(points(j))-u(points(i)));
      end
     end%for 
   end%for
%
   out=L;
end