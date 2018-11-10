function out=I_Cby_frac_1(N,r,u_)
%
%
% Overview
% This function returns an operational matrix of Integral of Chebyshev(1st) functions. 
%     
%out = I_Cby_frac_1(N,r,u_)
%
%inputs: 
%------------------------------------------------------------ 
%| N   : integer          : From  Chebyshev(1st) sentence 0 |  
%|                          to Chebyshev(1st) sentence N    |
%| r   : double           : Order of integral               |
%| u_  : symbolic function: Shifting parameter              |          
%------------------------------------------------------------    
%
%Output:
%---------------------------------------------------------------    
%| out : [(N+1)x(N+1)] double : Operational Integral Chebyshev |       
%---------------------------------------------------------------    
% 
%
% Caution!:This function works for returning the integral
%  operational matrix for interval [a,b].
%

%
%%Reference: The Construction of Operational Matrix of Fractional
%Integration Using the Fractional Chebyshev Polynomials
%E. Fathizadeh1 · R. Ezzati1 · K. Maleknejad
%
%
f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);

%% now calculate the derivtive of the function

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
   if (i==0) 
    out(1,j+1)=c(0,j,a,b,r)/gamma(1+r);
    else 
     su=0;
 
     for k=0:i 
       su=su+c(k,j,a,b,r)*((-1)**(i-k)*i*gamma(i+k)*2**(2*k)*gamma(k+1))/...
       ((b-a)**(k)*gamma(2*k+1)*gamma(k+r+1)*...
       gamma(i-k+1));
      end%k
    out(i+1,j+1)=su;
    endif   
 
 end%j 
 end%i
 
 out=out';
 end%func
 %---------------------------------------------
function out=c(k,j,a,b,r)

if (j==0) 
 out=((b-a)**(k+r)*gamma(k+r+0.5))/(sqrt(pi)*gamma(k+r+1));
else
su=0;
for l=0:j 
su=su+(j*(b-a)**(k+r)/sqrt(pi))*...
((-1)**(j-l)*gamma(j+l)*2**(2*l+1)*...
gamma(l+k+r+0.5))/(gamma(j-l+1)*gamma(2*l+1)*gamma(l+k+r+1));
end%for 
out=su;
 endif

 
end 
 