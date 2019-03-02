
function Fokker_nonlinear


global delta_t
global N
global point
global teta
global A1
global A2
global A3
global A4
global A5
global C

alpha_=0.5;


N=7
teta=0.5;
delta_t=0.001
Maxstep=1/delta_t
 
output_precision(20);
 


point=Gnbr_zeros(N-1,alpha_+1)';
point=[0,(point+1)/2,1];

A1=eye(1,N+1);% for the Bounday Condition y(0,t)=0

A2=eye(N+1,N+1);% I
A2=A2(2:N,:);


A3=[zeros(1,N),1];%for the Bounday Condition  y(1,t)=exp(t)


D_1=D_lgr_Gnbr(N+1,alpha_,point', @(x) 2*x-1,'gauss_lbt');%Q_x
D_2=(0+diag(ones(N+1,1)*2)*D_1)*(diag(ones(N+1,1)*2)^(-1))*D_1;
A4=D_2(2:N,:);


A5=D_1(2:N,:);


C(1:(N+1),1)=point'.^2;
 
for step=2:Maxstep

tem1=ones(1,N-1)'-teta*(delta_t)*(ones(1,N-1)'.*(1./3)+...
     C(2:N,step-1).*(4./(point(2:N)'.^2))-...
     8./(point(2:N)').*A5*C(1:N+1,step-1)+...
     (2*A4*C(1:N+1,step-1)));

tem2=teta*(delta_t)*(point(2:N)'/3+2*A5*C(1:N+1,step-1));


A=      [A1;  diag(tem1)*A2-diag(tem2)*(A5);A3];


C(1:(N+1),step) = A\b(step);

end%for step

%------------------test----------

 
max_poi=25; 
poi=linspace(0,1,max_poi);


for ti=1:1:Maxstep/40
for xi=1:1:max_poi
  y(xi,ti)=lgr_Gnbr_(N+1,poi(xi),point, @(x) 2*x-1)*C(1:(N+1),ti*40)...
  -(poi(xi)^2*exp((ti*40)*delta_t));
  end
end

[px py]=ndgrid(poi,poi);

figure(1);
surf(py,px,abs(y))

end%function

%------------------------------------------------------
function out=b(step)

global delta_t;
global N
global point
global teta
global C
global A1
global A2
global A3
global A4
global A5

t=delta_t*step;

out(1,1)=0;
out(2:N,1)=(C(2:(N),step-1))+(delta_t)*(1-teta)*(   
              (ones(1,N-1)'.*(1./3)+...
              C(2:N,step-1).*(4./(point(2:N)'.^2))-...
              8./(point(2:N)').*A5*C(1:N+1,step-1)+...
              (2*A4*C(1:N+1,step-1))).*C(2:N,step-1)+...
              (point(2:N)'/3+2*A5*C(1:N+1,step-1)).*A5*C(1:N+1,step-1));
   
out(N+1,1)=exp(t);
 

end %function b
