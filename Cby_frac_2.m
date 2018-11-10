function out=Cby_frac_2(N,derv,X,u_)

%
%

% 
% Overview
% This function returns a Shifted fractional Chebyshev (2nd kind) matrix functions. 
%     
%%out = Cby_frac_3(N,derv,X,u_)
%
%inputs: 
%----------------------------------------------------------------- 
%| N     : integer          : From fractional Chebyshev(2nd kind)|
%|                            sentence 0 to fractional Chebyshev |
%|                            sentence number N                  |
%| derv  : double           : derivative order                   |
%|   X   : [1xm] double     : Inputs of u_(x) in                 | 
%|                            fractional Chebyshev (2nd kind)    |
%|                            functions                          |
%|   u_  : symbolic function: Shifting parameter                 |          
%-----------------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted fractional Chebyshev       |
%|                            (2nd kind) functions               |       
%-----------------------------------------------------------------    
% 
theta_=1;
out=Gnbr_frac(N,theta_,derv,X,u_);

end