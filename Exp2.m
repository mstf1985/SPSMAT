% 2*y''-y=4-x**2. the exact solution is y=x**2
% This is the "Example 2"  in documentation file. 
function Exp2

N=2;
alpha_=1;
beta_=1;

%collocation points
X=(1+jacobi_zeros(N+1,alpha_,beta_))/2;


H=F(X)*inv(jacobi_(N,alpha_,beta_,0,X,@(x)2*x-1)');

D_=D_jacobi(N,alpha_,beta_,@(x)2*x-1);
A=2*(D_*D_)-eye(N+1,N+1);

% Boundary condition y(0)=0
A(1,:)=jacobi_(N,alpha_,beta_,0,0,@(x)2*x-1);
H(1,1)=0;

% Ac=H'-->c=A\H'
c=A\H';

poi=linspace(0,1,20);

%Error ploting 
 for i=1:1:20
  y(i)=jacobi_(N,alpha_,beta_,0,poi(i),@(x)2*x-1)*c-(poi(i)**2);  
 end


figure (2);
plot(poi,abs(y),"-s", "markersize", 10)
title("Error illustration");
xlabel ("x");
saveas (2, "Pic_Exp2.eps");
end
%-----------------------------------

function out=F(X)
% F(X) return [f(X0),...,f(Xn)]
  for i=1:length(X)
    out(1,i)=4-X(i)**2;
  end%for

end






