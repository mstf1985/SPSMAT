
 function Exp0

 alpha_=1;
 beta_=3;
 N=4;
 a=0;
 b=1;

 X=[-1,jacobi_zeros(N-1,alpha_+1,beta_+1)',1]; 
 r=f(((b-a)/2).*X+(b+a)/2);
 w=((b-a)/2)*jacobi_w(N,alpha_,beta_,'gauss_lbt');
 r*w'
 end


 function out=f(X)
   for i=1:length(X)
     out(i)=3*X(i)**3+X(i)**2;
   end%for
 end
 
 