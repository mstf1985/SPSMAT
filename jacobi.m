
function out=jacobi(N,alpha_,beta_,derv,shify,X,shifx)

%
%
% Just for simplicity in this guidance, imagine a_:=alpha_ , b_:=beta_ , s_:=shify,b_:=shifx , d:=derv
%
%x=[-1,1]-->t=[0,1] i.e the t=2x-1 then jacobi(1,a_,b_,d,s_*X+b_)=jacobi(1,a_,b_,d,2*X-1) and it 
%is going to be jacobi(1,a_,b_,d,2,X,-1) in function calling.
%
%if X contains only one element, jacobi(N,a_,b_,d,s_*X,b_) return a sequence(Vector)  like
%[ jacobi(1,a_,b_,d,s_*X+b_), jacobi(1,a_,b_,d,s_*X+b_),..., jacobi(N,a_,b_,d,s_*X+b_)] i.e N+1 sentences.

%  if X contains more than one element i.e a vector, jacobi(N,alpha_,beta_,d,shify*X) return a Matrix as like as
%[[ jacobi(1,a_,b_,d,s_*X(1)+b_), jacobi(1,a_,b_,d,s_*X(1)+b_),..., jacobi(N,a_,b_,d,s_*X(1)+b_)],
% [ jacobi(1,a_,b_,d,s_*X(2)+b_), jacobi(1,a_,b_,d,s_*X(2)+b_),..., jacobi(N,a_,b_,d,s_*X(2)+b_)],
% .
% .
% .
% [ jacobi(1,a_,b_,d,s_*X(N)+b_), jacobi(1,a_,b_,d,s_*X(N)+b_),..., jacobi(N,a_,b_,d,s_*X(N)+b_)]]
%i.e a (N+1) by (N+1) Matrix .
%
%Clearly for d=m where m>=0 that means the derivative of standard jacobi polynomials of order m
%
%



alpha_=alpha_+derv;
beta_=beta_+derv;

for i=1:length(X)
if N<0 
   out =0 ;
   return
 elseif N==0 
     out(i,:)=1;
     return   
else 
   Ja(1)=1;
   Ja(2)=(0.5)*(alpha_+beta_+2)*((X(i)*shify)+shifx)+(0.5)*(alpha_-beta_);
   for n=1:1:N-1
     a=((2*n+alpha_+beta_+1)*(2*n+alpha_+beta_+2))/((2*(n+1)*(n+alpha_+beta_+1)));
     b=((alpha_^2-beta_^2)*(2*n+alpha_+beta_+1))/(2*(n+1)*(n+alpha_+beta_+1)*(2*n+alpha_+beta_));
     c=((n+alpha_)*(n+beta_)*(2*n+alpha_+beta_+2))/((n+1)*(n+alpha_+beta_+1)*(2*n+alpha_+beta_));
     Ja(n+2)=((a)*(X(i)*shify+shifx)+(b))*Ja(n+1)-(c)*Ja(n);
   end%for

 end%else
out(i,:)=Ja;
 if derv~=0 
 for n=1:1:N+1
   if (n-derv>0)
    Jaa(n)=(shify)^derv*((gamma((alpha_-derv)+(beta_-derv)+(n-1)+1+derv))/(2^derv*gamma((alpha_-derv)+(beta_-derv)+(n-1)+1)))*Ja(n-derv);
   else 
    Jaa(n)=0;
   end %if
 end%for 
 out(i,:)=Jaa; 
 end %if
end %for x
 
 end%function
 
 
 
 
 
 