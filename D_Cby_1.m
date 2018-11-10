function out=D_Cby_1(N,u_)
%
% Overview
% This function returns an operational matrix of derivative of Chebyshev(1st) functions. 
%     
%out = D_Cby_1(N,u_) 
%
%inputs: 
%-------------------------------------------------------------- 
%| N     : integer          : From  Chebyshev(1st) sentence 0 |  
%|                            to Chebyshev(1st) sentence N    |
%|                            Chebyshev(1st) functions        |
%| u_    : symbolic function: Shifting parameter              |          
%--------------------------------------------------------------    
%
%Output:
%----------------------------------------------------------------    
%| out   : [(N+1)x(N+1)] double : Derivative Operational matrix |       
%----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
%Reference: A Chebyshev spectral method based on operational matrix for initial and
%boundary value problems of fractional order
%E.H. Doha a, A.H. Bhrawy b,∗, S.S. Ezz-Eldien
%
%
%
%
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
b=(1-B)/A;

out=zeros(N+1,N+1);

for i=1:N+1
for j=1:N+1
 if (j>=i) out(i,j)=0;
 elseif (mod(i,2)==0) && (j==1) out(i,j)=i-1;  
 elseif (mod(i,2)==1) && (j~=1) && (mod(j,2)==0) out(i,j)=2*(i-1);  
 elseif (mod(i,2)==0) && (mod(j,2)==1) && j!=1   out(i,j)=2*(i-1);
 endif
end
end

out=(2/(b-a))*out';
 end%func
