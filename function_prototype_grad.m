% find the gradient of one of the "Fair " (used for robust estimation according to the paper by hero blatta and gauchman
%functions which make up the total
% objective function, each function is different from the other one by its
% measurement , one of the values in Y_l passed as y here
% f_l(x) = (1/L) g(x -y_l) 
% g(x) = c^2 [ ( |x|/c ) - log ( 1 + ( |x|/c ) ) ]

function [grad] = function_prototype_grad( x, y, L ,c )
g_arg = x-y ;
abs_val = abs(g_arg);
grad = (1/L) * (   g_arg / ( 1 + (abs_val/c) ) ) ;

 %grad = 2* ( x-y)/L ; % this is the gradient for the the quadratic
 %functions used for experimenting
end