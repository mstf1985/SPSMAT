% This is the "Example 5"  in documentation file. 
 function Exp5

 alpha_=1;
 beta_=0;
 N=5;
 
 
 u_=@(x)x/2;
 X=[-2,2*jacobi_zeros(N-1,alpha_+1,beta_+1)',2]; 
 
 
 H=F(X);
 
 D_1=D_lgr_jacobi(N+1,alpha_,beta_,X,u_,'gauss_lbt');
 D_2=   (diag(ones(N+1,1)*(1/2.0))*D_1)...
       *(diag(ones(N+1,1)*(1/2.0))^(-1))*D_1;%
      
 A=2*(D_2)-eye(N+1,N+1);

 % Boundary condition y(-2)=0
 A(1,:)=[1,zeros(1,N)];
 H(1,1)=4;
  
 A(N+1,:)=[zeros(1,N),1];
 H(1,N+1)=4;
 
 % Ac=H'-->c=A\H'
 c=A\H';

 poi=linspace(0,1,20);

 %Error ploting 
 for i=1:1:20
  y(i)=lgr_Jacobi_(N+1,poi(i),X, u_)*c-(poi(i)**2);  
 end


 figure (5);
 plot(poi,abs(y),"-s", "markersize", 10)
 title("Error illustration");
 xlabel ("x");
 saveas (5, "Pic_Exp5.eps");
end
%-----------------------------------

function out=F(X)
% F(X) return [f(X0),...,f(Xn)]
  for i=1:length(X)
    out(1,i)=4-X(i)**2;
  end%for

end






