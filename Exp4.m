% This function solves a differntial equation by a jacobi integral operational matrix.
% the problem is intoduced in documentation as Example 4.
function Exp4

 N=1;
 alpha_=1;
 beta_=1;
 v=0.2;

 u_=@(x)2*x-1;
 X=(1+jacobi_zeros(N+1,alpha_,beta_))/2;


 H=f(X,v)*inv(jacobi_(N,alpha_,beta_,0,X,u_)');


 A=I_jacobi_frac(N,alpha_,beta_,1  ,u_)+...
   I_jacobi_frac(N,alpha_,beta_,1-v,u_);

 %Ac=H-->c=?
 c=H*inv(A);

poi=linspace(0,1,20);

%Error ploting 
 for i=1:1:20
  y(i)=c*jacobi_(N,alpha_,beta_,0,poi(i),u_)'-(-1+2*poi(i));  
 end

figure (4);
plot(poi,abs(y),"-s", "markersize", 10)
title("Error illustration");
xlabel ("x");
saveas (4, "Pic_Exp4.eps");
end
%-----------------------------------

function out=f(X,v)

  for i=1:length(X)
    out(1,i)=2*X(i)^(2-v)/gamma(3-v)-X(i)^(1-v)/gamma(2-v)+X(i)^2-X(i);
  end%for

end






