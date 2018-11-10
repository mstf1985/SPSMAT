function w=jacobi_w(n,alpha_,beta_,type)
%
% Overview
% This function returns Gauss/Gauss-Radau/Gauss-Lobatto 
% weights for Jacobi polynomials over [-1,1].
%    
%out = jacobi_w(n,alpha_,beta_,type) 
%
%inputs: 
%------------------------------------------------------ 
%| n     : integer          : Jacobi sentence n       |  
%| alpha_: double           : Jacobi parameter        |
%| beta_ : double           : Jacobi parameter        |
%| type  : string           : Gauss type              |
%------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------    
%| out   : [1x(N+1)] double :  Weight of Jacobi polynomials|       
%-----------------------------------------------------------    
% 
 if strcmp(type,'gauss')
    X=jacobi_zeros(n+1,alpha_,beta_)';
    J=jacobi_(n+1,alpha_,beta_,1,X,@(x)x);
    for i=1:n+1
       w(i)=1./(J(i,n+2)**2*(1-X(i).^2))*g(n,alpha_,beta_);
    end %for
 
 elseif strcmp(type,'gauss_rdu')
     X=[-1,jacobi_zeros(n,alpha_,beta_+1)'];
     for i=1:n+1
       if (i==1)
          w(i)=(gamma(n+alpha_+1)*(beta_+1)*gamma(beta_+1)**2*gamma(n+1))...
          /(gamma(n+beta_+2)*gamma(n+alpha_+beta_+2));
       else
          
          J=jacobi_(n,alpha_,beta_+1,1,X,@(x)x); 
          w(i)=1./(J(i,n+1)**2*(1-X(i))*(1+X(i)).^2)*g(n-1,alpha_,beta_+1);
       endif
     end %for
 
 elseif type='gauss_lbt'
     X=[-1,jacobi_zeros(n-1,alpha_+1,beta_+1)',1];
     for i=1:n+1
       if (i==1)
          w(i)=(gamma(n+alpha_+1)*(beta_+1)*gamma(beta_+1)**2*gamma(n)*2**(alpha_+beta_+1))...
          /(gamma(n+beta_+1)*gamma(n+alpha_+beta_+2));
       elseif (i==n+1)
          w(i)=(gamma(n+beta_+1)*(alpha_+1)*gamma(alpha_+1)**2*gamma(n)*2**(alpha_+beta_+1))...
          /(gamma(n+alpha_+1)*gamma(n+alpha_+beta_+2));
       else
          
          J=jacobi_(n-1,alpha_+1,beta_+1,1,X,@(x)x); 
          w(i)=1./(J(i,n)**2*(1-X(i).^2).^2)*g(n-2,alpha_+1,beta_+1);
       endif
     end %for
 endif
 
end%function


%............................
function out=g(n,a_,b_)


out= (gamma(a_+n+2)*gamma(b_+n+2)*2**(a_+b_+1))/(gamma(a_+b_+n+2)*gamma(n+2));

end





