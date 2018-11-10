%
%D_lgr_jacobi(5,1,0,0.1, @(x) 2*x-1,'gauss')
%D_lgr_jacobi(5,1,0,[0.1, 0.4], @(x) 2*x-1,'gauss')
%--------------------------------------------
%
% u= @(x) 2*x-1
%D_lgr_jacobi(5,1,0,[0.1, 0.4], u,'gauss')
%
%
% u= @(x) x
%D_lgr_jacobi(5,1,0,jacobi_zeros(5,1,0), u,'gauss')
%
%

% u= @(x) 2*x-1
% D_lgr_jacobi(7,0,1,(jacobi_zeros(7,0,1)+1)/2, @(x) 2*x-1,'gauss')
%
%
%D_lgr_jacobi(6,0,1,[-1;jacobi_zeros(5,0,2)], @(x) x,'gauss_rdu_f')
%
%D_lgr_jacobi(6,0,1,[0;(jacobi_zeros(5,0,2)+1)/2], @(x) 2*x-1,'gauss_rdu_f')
%
%
%D_lgr_jacobi(6,0,1,[(jacobi_zeros(5,0,2)+1)/2;1], @(x) 2*x-1,'gauss_rdu_l')
%
%D_lgr_jacobi(5,0,1,[jacobi_zeros(4,0,2);1], @(x) x,'gauss_rdu_l')
%
%
%%D_lgr_jacobi(6,0,1,[-1;jacobi_zeros(4,1,2);1], @(x) x,'gauss_lbt')
%
%D_lgr_jacobi(6,0,1,[0;(1+jacobi_zeros(4,1,2))/2;1], @(x) 2*x-1,'gauss_lbt')
%--------------------------------------------
%Rule 1#: X must include N elements
%

function out=D_lgr_jacobi(N,a_,b_,X, u_,gu)
output_precision(5);

 
f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);

% now calculate the derivative of the function

ffd = diff(u, x);
fffd=diff(ffd,x);
% and convert it back to an Anonymous function

%u=function_handle(fd);
u_derv = function_handle(ffd);
u_derv_2 = function_handle(fffd);
u_derv(X);

u(X);
D=zeros(N,N);

%--------------------------------------Gauss------------------------------------
if isequal(gu ,'gauss')

n=N-1;
% P=J_(n+1)-->P'=u'J_n
J_n=jacobi(n,a_+1,b_+1,0,1,u(X),0);

J_n(:,n+1);
u_derv(X(2));

   for j=1:1:n+1;
    for k=1:1:n+1;
     if k!=j 
      D(k,j)=(u_derv(X(k))/(u(X(k))-u(X(j))))*(J_n(k,n+1)/J_n(j,n+1));
     else
      D(k,k)=(u_derv(X(k))) *(a_-b_+(2+a_+b_)*u(X(k)))/(2*(1-u(X(k))^2)) ;
     endif

    end 
  end; 

%----------------------------------Gauss Radau----------------------------------
elseif isequal(gu ,'gauss_rdu_f')

n=N-1;
% P=(u-u0)J_(n)-->P'=u'J_n+(u-u0)J_n_1
J_n=jacobi(n,a_,b_+1,0,1,u(X),0);
J_n_1=jacobi(n-1,a_+1,b_+2,0,1,u(X),0);

J_n_2=jacobi(n-2,a_+2,b_+3,0,1,u(X),0);
%J_n_3=jacobi(n-3,a_+3,b_+3,0,1,u(X),0);

J_n_1(:,n);
u_derv(X(2));

   for j=1:1:n+1;
    for k=1:1:n+1;
    
     if k==j &&  j==1
       D(k,j)=(n*u_derv(X(j))*(a_+b_+n+2))/(-2*(b_+2)); 
     elseif j==1 && 2<=k<=n+1 
      D(k,j)=((gamma(n+1)*gamma(b_+2))/(gamma(b_+n+2)*2*(-1)^n))*...
      (u_derv(X(k))*(a_+b_+n+2)*J_n_1(k,n));
     elseif k==1 && 2<=j<=n+1
       D(k,j)=(2*u_derv(X(k))*gamma(b_+n+2)*(-1)^n)/((u(X(j))-u(X(1)))*...
       (a_+b_+n+2)*J_n_1(j,n)*gamma(n+1)*gamma(b_+2)*(u(X(1))-u(X(j))));
     elseif k!=j &&  2<=k<=n+1 && 2<=j<=n+1
       D(k,j)=(u_derv(X(k))*(u(X(k))-u(X(1)))*J_n_1(k,n))/((u(X(j))-u(X(1)))*...
       J_n_1(j,n)*(u(X(k))-u(X(j))));
     elseif   2<=k<=n+1 && 2<=j<=n+1 && k==j 
       D(k,j)=u_derv(X(j))/(u(X(k))-u(X(1)))+u_derv(X(j))*(a_+b_+n+3)*J_n_2(k,n-1)/...
      (4*J_n_1(j,n));
    
     endif

    end 
  end; 

  %------------------------------------------------------------------
  elseif isequal(gu ,'gauss_rdu_l')

n=N-1;
% P=(u-un)J_(n)-->P'=u'J_n+(u-un)J_n_1
J_n=jacobi(n,a_,b_+1,0,1,u(X),0);
J_n_1=jacobi(n-1,a_+1,b_+2,0,1,u(X),0);

J_n_2=jacobi(n-2,a_+2,b_+3,0,1,u(X),0);


J_n_1(:,n);
u_derv(X(2));

   for j=1:1:n+1;
    for k=1:1:n+1;
     
     if k==j &&  j==n+1  
       D(k,j)=(n*u_derv(X(j))*(a_+b_+n+2))/(-2*(a_+1)); 
     elseif j==n+1 && 1<=k<=n 
      D(k,j)=((gamma(n+1)*gamma(a_+1))/(gamma(a_+n+1)*2*(-1)^n))*...
      (u_derv(X(k))*(a_+b_+n+2)*J_n_1(k,n));
     elseif k==n+1 && 1<=j<=n
       D(k,j)=(2*u_derv(X(k))*gamma(a_+n+1))/((u(X(j))-u(X(n+1)))*...
       (a_+b_+n+2)*J_n_1(j,n)*gamma(n+1)*gamma(a_+1)*(u(X(n+1))-u(X(j))));
     elseif k!=j &&  1<=k<=n && 1<=j<=n
       D(k,j)=(u_derv(X(k))*(u(X(k))-u(X(n+1)))*J_n_1(k,n))/((u(X(j))-u(X(n+1)))*...
       J_n_1(j,n)*(u(X(k))-u(X(j))));
     elseif   1<=k<=n && 1<=j<=n && k==j 
       D(k,j)=u_derv(X(j))/(u(X(k))-u(X(n+1)))+u_derv(X(j))*(a_+b_+n+3)*J_n_2(k,n-1)/...
      (4*J_n_1(j,n));
    
     endif
     
    end 
  end; 
%----------------------------------Gauss Lobatto----------------------------------
elseif isequal(gu ,'gauss_lbt')
n=N-1
% P=(u-u0)J_(n-1)(u-un)-->P'=u'J_n_1)(u-un)+(u-u0)J_n_1)+c*u'(u-u0)(u-un)J_n_2
J_n_1=jacobi(n-1,a_+1,b_+1,0,1,u(X),0);
J_n_2=jacobi(n-2,a_+2,b_+2,0,1,u(X),0);
J_n_3=jacobi(n-3,a_+3,b_+3,0,1,u(X),0);

u_derv(X(2));

   for j=1:1:n+1;
    for k=1:1:n+1;
     
     if (j==1 && k>=2 && k<=n) 
       D(k,j)=(u_derv(X(k))*J_n_2(k,n-1)*(u(X(k))-u(X(n+1)))*gamma(n)*(a_+b_+n+2)*gamma(b_+2))...
       /(2*(u(X(1))-u(X(n+1)))*(-1)^(n-1)*gamma(b_+n+1));
      
     elseif ( k == n+1 && j==1   )
 
       D(k,j)=(u_derv(X(n+1))*gamma(a_+n+1)*gamma(b_+2))...
        /((u(X(1))-u(X(n+1)))*(-1)^(n-1)*gamma(b_+n+1)*gamma(a_+2));
     elseif k==1 && j>=2 && j<=n 
       D(k,j)=(-2*(-1)^(n-1)*u_derv(X(1))*(u(X(1))-u(X(n+1)))*gamma(b_+n+1))...
       /(gamma(b_+2)*gamma(n)*(a_+b_+n+2)*J_n_2(j,n-1)*(u(X(1))-u(X(j)))^2*(u(X(j))-u(X(n+1))));
     elseif k==1 && j==(n+1) 
      D(k,j)=(u_derv(X(1))*gamma(b_+n+1)*gamma(a_+2)*(-1)^(n-1))/((u(X(n+1))-u(X(1)))*gamma(a_+n+1)*gamma(b_+2));
     elseif j==n+1 && k>=2 && k<=n 
       D(k,j)=(u_derv(X(k))*J_n_2(k,n-1)*(u(X(k))-u(X(1)))*gamma(n)*(a_+b_+n+2)*gamma(a_+2))...
       /(2*(u(X(n+1))-u(X(1)))*gamma(a_+n+1));
     elseif j>=2 && j<=n && k>=2 && k<=n && j!=k
       D(k,j)=((u(X(k))-u(X(1)))*(u(X(k))-u(X(n+1)))*u_derv(X(k))*J_n_2(k,n-1))...
       /((u(X(j))-u(X(1)))*(u(X(j))-u(X(n+1)))*(u(X(k))-u(X(j)))*J_n_2(j,n-1));
     
     elseif j>=2 && j<=n &&  k==n+1;
       D(k,j)=(2*u_derv(X(n+1))*(u(X(n+1))-u(X(1)))*gamma(a_+n+1))...
       /(gamma(a_+2)*gamma(n)*(a_+b_+n+2)*J_n_2(j,n-1)*(u(X(1))-u(X(j)))*(u(X(j))-u(X(n+1)))^2);
     elseif j>=2 && j<=n && k>=2 && k<=n && j==k
       D(k,j)=u_derv(X(j))/((u(X(j))-u(X(1))))+u_derv(X(j))/((u(X(j))-u(X(n+1))))...
     +(u_derv(X(j))*(a_+b_+n+3)*J_n_3(j,n-2))/(4*J_n_2(j,n-1));
     elseif (j==1  ) && j==k
       D(k,j)=u_derv_2(X(j))/(2*u_derv(X(j)))-u_derv(X(j))/(u(X(n+1))-u(X(1)))+...
       (u_derv(X(j))*(a_+b_+n+2)*J_n_2(j,n-1))/(2*J_n_1(j,n));
      
      elseif ( j==n+1 ) && j==k
      D(k,j)=u_derv_2(X(j))/(2*u_derv(X(j)))+u_derv(X(j))/(u(X(n+1))-u(X(1)))+...
      (u_derv(X(j))*(a_+b_+n+2)*J_n_2(j,n-1))/(2*J_n_1(j,n));
     endif

    end 
  end 
endif % if guass/gauss_rdu/gauss_lbt
%-------------------------------------------------------------------------------  
out=D;
end%func