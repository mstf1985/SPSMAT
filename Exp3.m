% y(1.5)(x)+y(2)(x)+y=1+x. the exact solution is y=x+1
% This is the "Example 3"  in documentation file. 
function tee

N=3;
alpha_=0;
beta_=1;


X=(1+jacobi_zeros(N+1,alpha_,beta_))/2;


H=F(X)*inv(jacobi_(N,alpha_,beta_,0,X,@(x)2*x-1)');

D_3_2=D_jacobi_frac(N,alpha_,beta_,1.5,@(x)2*x-1);
D_2=  D_jacobi_frac(N,alpha_,beta_,2,  @(x)2*x-1);

A=D_3_2+D_2+eye(N+1,N+1);

% Boundary condition y(0)=1, y'(0)=1
A(1,:)=jacobi_(N,alpha_,beta_,0,0,@(x)2*x-1);
H(1,1)=1;

A(2,:)=jacobi_(N,alpha_,beta_,1,0,@(x)2*x-1);
H(1,2)=1;

% Ac=H'-->H'=?
c=A\H';

poi=linspace(0,1,20);

%Error ploting 
 for i=1:1:20
  y(i)=jacobi_(N,alpha_,beta_,0,poi(i),@(x)2*x-1)*c-(1+poi(i));  
 end



figure (3);
plot(poi,abs(y),"-s", "markersize", 10)
title("Error illustration");
xlabel ("x");
saveas (3, "Pic_Exp3.eps");
end
%-----------------------------------

function out=F(X)
% F(X) return [f(X0),...,f(Xn)]
  for i=1:length(X)
    out(1,i)=1+X(i);
  end%for

end






