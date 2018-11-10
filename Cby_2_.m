function out=Cby_2_(N,derv,X,u_)

% 
% Overview
% This function returns a Shifted Chebyshev (2nd kind) matrix functions. 
%     
%out = Cby_2_(N,derv,X,u_)
%
%inputs: 
%------------------------------------------------------ 
%| N     : integer          : From Chebyshev(2nd kind)|
%|                            sentence 0 to Chebyshev |
%|                            sentence number N       |
%| derv  : integer          : derivative order        |
%|   X   : [1xm] double     : Inputs of u_(x) in      | 
%|                            Chebyshev (2nd kind)    |
%|                            functions               |
%|   u_  : symbolic function: Shifting parameter      |          
%------------------------------------------------------    
%
%Output:
%------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted Chebyshev       |
%|                            (2nd kind) functions    |       
%------------------------------------------------------    
% 
theta_=1;
out=Gnbr_(N,theta_,derv,X,u_);

end