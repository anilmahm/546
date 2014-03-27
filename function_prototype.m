% the sub functions, the functions whose sum make up the total objective
% function, the y_l value or the measurement used differentiates one function from another,
% the value is passed as y 
% f_l(x) = (1/L) g(x -y_l) 
% g(x) = c^2 [ ( |x|/c ) - log ( 1 + ( |x|/c ) ) ]
function [ val ] = function_prototype( x, y, L ,c )
g_arg = x-y ;
abs_val = abs(g_arg);
val = (1/L) * c^2 *  (  ( abs_val /c )  -  log ( 1 + ( abs_val/c ) )  );
end
