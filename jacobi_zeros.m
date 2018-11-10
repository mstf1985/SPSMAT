function roots = jacobi_zeros( N, alpha_, beta_ )
%
%
%Overview
%jacobi_zeros( N, alpha_, beta_ ) is a function returing N roots 
%of the N-th sentence of  Standard orthogonal Jacobi polynomials
%over the interval [-1,1].
%
%inputs: 
%-------------------------------------------------------   
%| N     : integer          : Jacobi sentence N        |  
%| alpha_: double           : Jacobi parameter         |
%| beta_ : double           : Jacobi parameter         | 
%-------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [Nx1] double : Jacobi zeros                |       
%------------------------------------------------------    
%
%%Reference: Stoer J and Bulirsch R. Introduction to numerical analysis.
%Springer Science and Business Media; 2013

output_precision(7);
for j=0:1:N-1
 a(j+1)=(beta_^2-alpha_^2)/((2*j+alpha_+beta_+eps(0))*(2*j+alpha_+beta_+2+eps(0)));
 b(j+1)=(4*j*(j+alpha_)*(j+beta_)*(j+alpha_+beta_))/...
 ((2*j+alpha_+beta_-1+eps(0))*((2*j+alpha_+beta_+eps(0))^2)*...
 (2*j+alpha_+beta_+1+eps(0)));
end

A=zeros(N,N);

for i=1:N
 A(i,i)=a(i);
   if i<N
    A(i,i+1)=sqrt(b(i+1));
    A(i+1,i)=sqrt(b(i+1));
   end
end

roots=sort(eig(A));

end