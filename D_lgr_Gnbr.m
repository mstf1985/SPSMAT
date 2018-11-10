
function out=D_lgr_Gnbr(N,theta,X, u_,gu)

a_=theta-0.5; ;
b_=theta-0.5; 
out=D_lgr_jacobi(N,a_,b_,X, u_,gu);

end