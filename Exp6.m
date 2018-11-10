% y(1.5)(x)+y(2)(x)+y=1+x. the exact solution is y=x+1 with gengbaure
function tee

N=10;
theta_=1;
u=@(x)2*x-1;

X=(1+Gnbr_zeros(N+1,theta_))/2;


H=F(X);

D_0  =Gnbr_(N,theta_,0  ,X,u);
D_3_2=Gnbr_frac(N,theta_,1.5,X,u);
D_2  =Gnbr_(N,theta_,2  ,X,u);


A=D_3_2+D_2+D_0;

% Boundary condition y(0)=1, y'(0)=1
A(1,:)=Gnbr_(N,theta_,0,0,u);
H(1,1)=1;
A(2,:)=Gnbr_(N,theta_,1,0,u);
H(1,2)=1;

% Ac=H'-->H'=?
c=A\H';

poi=linspace(0,1,20);

%Error ploting 
 for i=1:1:20
  y(i)=Gnbr_(N,theta_,0,poi(i),@(x)2*x-1)*c-(1+poi(i));  
 end


plot(poi,abs(y),"s-", "markersize", 10)
title("Error illustration")
hgsave ("f")
end
%-----------------------------------

function out=F(X)
% F(X) return [f(X0),...,f(Xn)]
  for i=1:length(X)
    out(1,i)=1+X(i);
  end%for

end






