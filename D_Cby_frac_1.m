function out=D_Cby_frac_1(N,derv,u_)
%
%
% Overview
% This function returns an operational matrix of derivative of Chebyshev(1st) functions. 
%     
%out = D_Cby_frac_1(N,r,u_)
%
%inputs: 
%------------------------------------------------------------ 
%| N   : integer          : From  Chebyshev(1st) sentence 0 |  
%|                          to Chebyshev(1st) sentence N    |
%| r   : double           : Order of derivative             |
%| u_  : symbolic function: Shifting parameter              |          
%------------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational derivative Chebyshev |       
%-----------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the derivative
%  operational matrix for interval [a,b].
%
%
%
%%Reference: A Chebyshev spectral method based on operational matrix for initial and
%boundary value problems of fractional order
%E.H. Doha a, A.H. Bhrawy b,∗, S.S. Ezz-Eldien
%
%
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

out=zeros(N+1,N+1);

 
 
 for i=0:N 
  for j=0:N 
   su=0;
 
   for k=ceil(derv):i 
    su=su+((-1)**(i-k)*2*i*gamma(i+k)*gamma(k-derv+0.5))/...
    (epp(j)*(b-a)**derv*gamma(k+0.5)*gamma(k-derv-j+1)*...
    gamma(k+j+1-derv)*gamma(i-k+1));
   end%k
out(i+1,j+1)=su; 
 end%j 
 end%i
 
 out=out';
 end%func
 
function out=epp(j)
 if (j==0) out=2;
  else out=1;
  endif
end 
