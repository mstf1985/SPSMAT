function out=D_lgr_Cby_3(N,theta,X, u_,gu)

a_=-0.5; ;
b_=0.5; 
out=D_lgr_jacobi(N,a_,b_,X, u_,gu);

end